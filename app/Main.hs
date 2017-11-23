{- |
Module      :  BMM simulator
Description :  Simulate populations using the boundary mutation model
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Enables simulation of sequence data for a different phylogenies using the
boundary mutation model. The output is a counts file.

* Changelog

-}

import qualified ArgParse                         as Args
import qualified BndModel                         as BMM
import qualified BndState                         as BMS
import qualified CFWriter                         as CF
import Control.Monad.Logger
import           Control.Monad.Random.Strict
import           Control.Monad.Trans.State.Strict
import qualified Data.Distribution                as D
import           Data.Maybe
import           Data.Random
-- import Data.Text (pack)
import qualified DNAModel                         as DNA
import qualified GammaRate                        as G
import           Numeric.LinearAlgebra
import qualified RateMatrix                       as RM
import qualified RTree                            as Tree
import qualified System.Environment               as Sys
import           System.IO
import qualified Transition                       as Trans

data BMMParams = BMMParams
  { progName     :: String
  , args         :: [String]
  , bmmArgs      :: Args.BMMArgs
  , gen          :: StdGen
  , logHandle    :: Handle
  , treeFilePath :: FilePath }

-- | The logging monad is wrapped around the actual state monad of the
-- simulator. However, this logging feature is not used at the moment (for
-- historical reasons).
type Simulation = LoggingT (StateT BMMParams IO)

-- | Initialize the parameters for the simulator (command line arguments etc.).
getParams :: IO BMMParams
getParams = do
  p    <- Sys.getProgName
  a    <- Sys.getArgs
  bmmA <- Args.parseBMMArgs
  let seedArg = Args.seed bmmA
  unless (seedArg == "random") $ setStdGen (read $ seedArg ++ " 1")
  g    <- getStdGen
  let l  = Args.outFileName bmmA ++ ".log"
      t  = Args.outFileName bmmA ++ ".tree"
  lh <- openFile l WriteMode
  return $ BMMParams p a bmmA g lh t

-- Use something like this for printing info. Problem: the nice variables
-- defined in simulate are not accessible. Maybe include them in the data
-- structure?
-- printInfo :: BMMParams -> IO ()
-- printInfo = do ...

getCommandLineStr :: String -> [String] -> String
getCommandLineStr n as = unlines
  [ "Boundary mutation model simulator version 0.1.0.1."
  , "Command line: " ++ n ++ " " ++ unwords as ]

getGeneratorStr :: String -> StdGen -> String
getGeneratorStr s g = unlines
  [ "Seed: " ++ s ++ " (Generator: " ++ show g ++ ")" ]

getHeadlineStr :: String -> String
getHeadlineStr h = unlines
  [ ""
  , "--"
  , h ]

getFileNamesStr :: String -> String -> String
getFileNamesStr d t = unlines
  [ "Counts file name: " ++ d
  , "Tree file name: " ++ t
  , "Branch lengths are measured in mutations and frequency shifts." ]

printTreeToFile :: Tree.RTree String Tree.BranchLn -> FilePath -> IO ()
printTreeToFile tree fn = do
  fh <- openFile fn WriteMode
  hPutStrLn fh $ Tree.toNewick tree
  hClose fh

logStr :: String -> Simulation ()
logStr s = do
  params <- lift get
  liftIO $ hPutStr (logHandle params) s
  liftIO $ putStr s

simulate :: Simulation ()
simulate = do
  params <- lift get
  -- Mutation model.
  let bmmA       = bmmArgs params
      dnaModelSpec = Args.dnaModelSpec bmmA
      dnaModel   = DNA.rateMatrix dnaModelSpec
  -- Gamma rate heterogeneity, handled with the Maybe monad.
  let gammaNCat  = Args.gammaNCat bmmA
      gammaShape = Args.gammaShape bmmA
      gammaMeans = liftM2 G.getMeans gammaNCat gammaShape
  -- Boundary mutation model.
  let popSize          = Args.popSize bmmA
      heterozygosity   = Args.heterozygosity bmmA
      bmm              = BMM.createBMM dnaModel popSize heterozygosity
      bmmStationaryGen  = Trans.stationaryDistToGenerator (BMM.bmmStationaryDist bmm)
  -- Tree.
  let treeHeight             = Args.treeHeight bmmA
      treeType               = Args.treeType bmmA
      maybeTreeYuleRate = Args.treeYuleRate bmmA
      (treeSubs, scenario)   = case treeType of
                   "ILS"  -> Tree.ils treeHeight
                   "Yule" -> fst $ sampleState (Tree.yule treeHeight recipRate) (gen params)
                     where recipRate = fromMaybe (error "No Yule reciprocal speciation rate specified.")
                                       maybeTreeYuleRate
                   _      -> error $ "Tree type not recognized: " ++ treeType
      treeBMM     = BMM.scaleTreeToBMM bmm treeSubs
      popNames   = Tree.getLeaves treeBMM
      treePrb    = Trans.branchLengthsToTransitionProbs (BMM.bmmRateMatrix bmm) treeBMM
      treeGen    = Trans.treeProbMatrixToTreeGenerator treePrb
  -- Other options.
  let nSites       = Args.nSites bmmA
      fileName     = Args.outFileName bmmA

  -- $(logError) $ pack "An error ocurred."

  -- Output.
  logStr $ getCommandLineStr (progName params) (args params)

  logStr $ getHeadlineStr "General options."
  logStr $ getGeneratorStr (Args.seed bmmA) (gen params)
  logStr $ getFileNamesStr fileName (treeFilePath params)
  logStr $ "Number of simulated sites: " ++ show nSites ++ "\n"

  logStr $ getHeadlineStr "Boundary mutation model options."
  logStr $ BMM.getBMMInfoStr bmm gammaShape gammaMeans

  logStr $ getHeadlineStr "Tree options."
  logStr $ Tree.getTreeStr scenario treeSubs treeBMM

  -- Also output tree to special file.
  liftIO $ printTreeToFile treeBMM (treeFilePath params)

  logStr $ getHeadlineStr "Performing simulation."
  -- Prepare output file.
  treeHandle <- liftIO $ CF.open fileName
  liftIO $ CF.writeHeader treeHandle nSites popNames
  -- Simulation helpers.
  let toStates = map (BMS.rmStateToBMState popSize . snd) :: [(a, RM.State)] -> [BMS.State]
      writer   = CF.writeLine treeHandle "SIM"    :: CF.Pos -> CF.DataOneSite -> IO ()
  -- Simulation; loop over positions.
  let simAndPrintOneSite pos =
        if isNothing gammaShape then
          evalRandIO (Trans.simulateSite bmmStationaryGen treeGen)
          >>= writer pos . toStates
        else
          -- Gamma rate heterogeneity is activated.
          evalRandIO (Trans.simulateSiteGen uniformGen bmmStationaryGens treesGen)
          >>= writer pos . toStates
        where nCat               = fromMaybe (error "The number of gamma rate categories was not given.") gammaNCat
              uniformGen         = D.fromDistribution $ D.uniform [0 .. nCat - 1]
              means              = fromMaybe (error "No gamma shape parameter given.") gammaMeans
              mutationModel      = BMM.bmmMutModel bmm
              mutationModels     = [ DNA.dnaModelScale s mutationModel | s <- means ]
              bmms               = [ BMM.createBMMNormalized m popSize | m <- mutationModels ]
              bmmStationaryDists = map BMM.bmmStationaryDist bmms
              bmmStationaryGens  = map Trans.stationaryDistToGenerator bmmStationaryDists
              bmmRateMatrices    = map BMM.bmmRateMatrix bmms
              treesPrb           = [ Trans.branchLengthsToTransitionProbs b treeBMM | b <- bmmRateMatrices ]
              treesGen           = map Trans.treeProbMatrixToTreeGenerator treesPrb
  liftIO $ forM_ [1..nSites] simAndPrintOneSite

  -- Done.
  logStr "Done.\n"
  liftIO $ hClose treeHandle
  liftIO $ hClose (logHandle params)

main :: IO BMMParams
main = getParams >>= execStateT (runStderrLoggingT simulate)
