{-# LANGUAGE TemplateHaskell #-}

{- |
Module      :  BMM simulator
Description :  Simulate populations using the boundary mutation model
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Enables simulation of sequence data for a predefined phylogeny using the
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

-- TODO: Read in specific mutation model.

data BMMParams = BMMParams
  { progName :: String
  , args     :: [String]
  , bmmArgs  :: Args.BMMArgs
  , gen      :: StdGen }

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
  return $ BMMParams p a bmmA g

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

getBMMInfoStr :: Int
              -> Double
              -> BMM.MutModel
              -> RM.StationaryDist
              -> Double
              -> Maybe Double
              -> Maybe [Double]
              -> String
getBMMInfoStr n h m f k ma mrs = unlines $
  [ "Population size: " ++ show n
  , "Heterozygosity: " ++ show h
  , "Mutation model matrix: "
  , show m
  , "This corresponds to state frequencies (A, C, G, T): " ++ show f
  , "And a kappa value of: " ++ show k
  , "Gamma rate heterogeneity: " ++ show (isJust ma) ]
  ++ gammaShape ++ gammaMeans
  where
    gammaShape = maybe [] (\a -> ["Shape parameter: " ++ show a]) ma
    gammaMeans = maybe [] (\rs -> ["This corresponds to uniformly distributed rates: " ++ show rs]) mrs

getTreeStr :: Tree.Scenario
           -> Tree.RTree String Tree.BranchLn
           -> Tree.RTree String Tree.BranchLn
           -> String
getTreeStr s trSubs trBMM = show s ++ unlines
  [ "Species tree in average number of substitutions: " ++ Tree.toNewick trSubs
  , "Species tree in  mutations and frequency shifts: " ++ Tree.toNewick trBMM ]

printTreeToFile :: Tree.RTree String Tree.BranchLn -> FilePath -> IO ()
printTreeToFile tree fn = do
  fh <- openFile fn WriteMode
  hPutStrLn fh $ Tree.toNewick tree
  hClose fh

simulate :: Simulation ()
simulate = do
  params <- lift get
  -- Mutation model.
  let bmmA       = bmmArgs params
      stateFreqs = Args.stateFreqs bmmA
      kappa      = Args.kappa bmmA
      hkyModel   = DNA.rateMatrix stateFreqs (DNA.HKY kappa)
  -- Gamma rate heterogeneity, handled with the Maybe monad.
  let gammaNCat  = Args.gammaNCat bmmA
      gammaShape = Args.gammaShape bmmA
      gammaMeans = liftM2 G.getMeans gammaNCat gammaShape
  -- Boundary mutation model.
  let popSize          = Args.popSize bmmA
      heterozygosity   = Args.heterozygosity bmmA
  -- TODO: This may be hidden in BndModel.hs, so that it is cleaner here.
      mutationModel    = BMM.normalizeToTheta
        (DNA.dnaRateMatrix hkyModel) stateFreqs popSize heterozygosity
      bmRateMatrix     = BMM.normalizedRateMatrix mutationModel stateFreqs popSize
      bmStationaryDist = BMM.stationaryDist mutationModel stateFreqs popSize
      bmStationaryGen  = Trans.stationaryDistToGenerator bmStationaryDist
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
      treeBMM     = BMM.scaleTreeToBMM popSize treeSubs
      popNames   = Tree.getLeaves treeBMM
      treePrb    = Trans.branchLengthsToTransitionProbs bmRateMatrix treeBMM
      treeGen    = Trans.treeProbMatrixToTreeGenerator treePrb
  -- Other options.
  let nSites       = Args.nSites bmmA
      fileName     = Args.outFileName bmmA
      treeFileName = fileName ++ ".tree"

  -- $(logError) $ pack "An error ocurred."

  -- Output.
  liftIO . putStr $ getCommandLineStr (progName params) (args params)

  liftIO . putStr $ getHeadlineStr "General options."
  liftIO . putStr $ getGeneratorStr (Args.seed bmmA) (gen params)
  liftIO . putStr $ getFileNamesStr fileName treeFileName
  liftIO . putStrLn $ "Number of simulated sites: " ++ show nSites

  liftIO . putStr $ getHeadlineStr "Boundary mutation model options."
  liftIO . putStr $ getBMMInfoStr popSize heterozygosity mutationModel stateFreqs kappa gammaShape gammaMeans

  liftIO . putStr $ getHeadlineStr "Tree options."
  liftIO . putStr $ getTreeStr scenario treeSubs treeBMM

  -- Also output tree to special file.
  liftIO $ printTreeToFile treeBMM treeFileName

  liftIO . putStr $ getHeadlineStr "Performing simulation."
  -- Prepare output file.
  fileHandle <- liftIO $ CF.open fileName
  liftIO $ CF.writeHeader fileHandle nSites popNames
  -- Simulation helpers.
  let toStates = map (BMS.rmStateToBMState popSize . snd) :: [(a, RM.State)] -> [BMS.State]
      writer   = CF.writeLine fileHandle "SIM"    :: CF.Pos -> CF.DataOneSite -> IO ()
  -- Simulation; loop over positions.
  let simAndPrintOneSite pos =
        if isNothing gammaShape then
          evalRandIO (Trans.simulateSite bmStationaryGen treeGen)
          >>= writer pos . toStates
        else
          -- Gamma rate heterogeneity is activated.
          evalRandIO (Trans.simulateSiteGen uniformGen bmStationaryGens treesGen)
          >>= writer pos . toStates
        where nCat              = fromMaybe (error "The number of gamma rate categories was not given.") gammaNCat
              uniformGen        = D.fromDistribution $ D.uniform [0 .. nCat - 1]
              means             = fromMaybe (error "No gamma shape parameter given.") gammaMeans
              mutationModels    = [ scale s mutationModel | s <- means ]
              bmRateMatrices    = [ BMM.normalizedRateMatrix m stateFreqs popSize | m <- mutationModels ]
              bmStationaryDists = [ BMM.stationaryDist m stateFreqs popSize | m <- mutationModels ]
              bmStationaryGens  = map Trans.stationaryDistToGenerator bmStationaryDists
              treesPrb          = [ Trans.branchLengthsToTransitionProbs b treeBMM | b <- bmRateMatrices ]
              treesGen          = map Trans.treeProbMatrixToTreeGenerator treesPrb
  liftIO $ forM_ [1..nSites] simAndPrintOneSite

  -- Done.
  liftIO $ CF.close fileHandle
  liftIO $ putStrLn "Done."

main :: IO BMMParams
main = getParams >>= execStateT (runStderrLoggingT simulate)
