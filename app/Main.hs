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
import qualified BndModel                         as BM
import qualified BndState                         as BS
import qualified CFWriter                         as CF
-- TODO: Use proper logger instead of a myriad of 'liftIO'.
-- import qualified Control.Monad.Logger        as Log
import           Control.Monad.Random.Strict
import           Control.Monad.Trans.State.Strict
import qualified Data.Distribution                as D
import           Data.Maybe
import           Data.Random
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

type Simulation = StateT BMMParams IO

getParams :: IO BMMParams
getParams = do
  p    <- Sys.getProgName
  a    <- Sys.getArgs
  bmmA <- Args.parseBMMArgs
  let seedArg = Args.seed bmmA
  unless (seedArg == "random") $ setStdGen (read $ seedArg ++ " 1")
  g    <- getStdGen
  return $ BMMParams p a bmmA g

simulate :: Simulation ()
simulate = do
  params <- get
  -- Mutation model.
  let bmmA       = bmmArgs params
      stateFreqs = Args.stateFreqs bmmA
      kappa      = Args.kappa bmmA
      hkyModel   = DNA.rateMatrixHKY stateFreqs kappa
  -- Gamma rate heterogeneity, handled with the Maybe monad.
  let gammaNCat  = Args.gammaNCat bmmA
      gammaShape = Args.gammaShape bmmA
      gammaMeans = liftM2 G.getMeans gammaNCat gammaShape
  -- Boundary mutation model.
  let popSize          = Args.popSize bmmA
      heterozygosity   = Args.heterozygosity bmmA
      mutationModel    = BM.normalizeToTheta
        hkyModel stateFreqs popSize heterozygosity
      bmRateMatrix     = BM.normalizedRateMatrix mutationModel stateFreqs popSize
      bmStationaryDist = BM.stationaryDist mutationModel stateFreqs popSize
      bmStationaryGen  = Trans.stationaryDistToGenerator bmStationaryDist
  -- Tree.
  let treeHeight             = Args.treeHeight bmmA
      treeType               = Args.treeType bmmA
      maybeTreeYuleRecipRate = Args.treeYuleRecipRate bmmA
      treeSubs = case treeType of
                   "ILS"  -> Tree.ils treeHeight
                   "Yule" -> fst $ sampleState (Tree.yule treeHeight recipRate) (gen params)
                     where recipRate = fromMaybe (error "No Yule reciprocal speciation rate specified.")
                                       maybeTreeYuleRecipRate
                   _      -> error $ "Tree type not recognized: " ++ treeType
      treeBM     = BM.scaleTreeToBMM popSize treeSubs
      popNames   = Tree.getLeaves treeBM
      treePrb    = Trans.branchLengthsToTransitionProbs bmRateMatrix treeBM
      treeGen    = Trans.treeProbMatrixToTreeGenerator treePrb
  -- Other options.
  let nSites       = Args.nSites bmmA
      fileName     = Args.outFileName bmmA
      treeFileName = fileName ++ ".tree"
  -- Output.
  liftIO $ putStrLn "Boundary mutation model simulator version 0.1.0.0."
  liftIO $ putStr "Command line: "
  liftIO $ putStrLn $ progName params ++ " " ++ unwords (args params)

  liftIO $ putStrLn ""
  liftIO $ putStrLn "--"
  liftIO $ putStrLn "General options."
  liftIO $ putStr $ "Seed: " ++ Args.seed bmmA
  liftIO $ putStrLn $ " (Generator: " ++ show (gen params) ++ ")"
  liftIO $ putStr "Number of simulated sites: "
  liftIO $ print nSites
  liftIO $ putStr "Data will be written to: "
  liftIO $ print fileName
  liftIO $ putStr "Species tree in mutations and frequency shifts will be written to:"
  liftIO $ print treeFileName

  liftIO $ putStrLn ""
  liftIO $ putStrLn "--"
  liftIO $ putStrLn "Boundary mutation model options."
  liftIO $ putStr "Population size: "
  liftIO $ print popSize
  liftIO $ putStr "Heterozygosity: "
  liftIO $ print heterozygosity
  liftIO $ putStrLn "Mutation model matrix:"
  liftIO $ print mutationModel
  liftIO $ putStr "This corresponds to state frequencies (A, C, G, T): "
  liftIO $ print stateFreqs
  liftIO $ putStr "And a kappa value of: "
  liftIO $ print kappa
  liftIO $ putStr "Gamma rate heterogeneity: "
  liftIO $ print $ isJust gammaShape
  liftIO $ when (isJust gammaShape) $ do
    putStr "Shape parameter: "
    print $ fromJust gammaShape
  liftIO $ when (isJust gammaMeans) $ do
    putStr "This corresponds to uniformly distributed rates: "
    print $ fromJust gammaMeans

  liftIO $ putStrLn ""
  liftIO $ putStrLn "--"
  liftIO $ putStrLn "Tree options."
  liftIO $ putStrLn $ "Tree type: " ++ treeType
  liftIO $ putStr "Tree height in average number of substitutions: "
  liftIO $ print treeHeight
  liftIO $ when (treeType == "Yule") $ do
    putStr "Reciprocal Yule speciation rate: "
    print $ fromMaybe (error "No reciprocal Yule speciation rate specified.") maybeTreeYuleRecipRate
  liftIO $ putStr "Species tree in average number of substitutions: "
  liftIO $ putStrLn $ Tree.toNewick treeSubs
  liftIO $ putStr "Species tree in  mutations and frequency shifts: "
  liftIO $ putStrLn $ Tree.toNewick treeBM

  -- Also output tree to special file.
  treeH <- liftIO $ openFile treeFileName WriteMode
  liftIO $ hPutStrLn treeH $ Tree.toNewick treeBM
  liftIO $ hClose treeH

  liftIO $ putStrLn ""
  liftIO $ putStrLn "--"
  liftIO $ putStrLn "Performing simulation."
  -- Prepare output file.
  fileHandle <- liftIO $ CF.open fileName
  liftIO $ CF.writeHeader fileHandle nSites popNames
  -- Simulation helpers.
  let toStates = map (BS.rmStateToBMState popSize . snd) :: [(a, RM.State)] -> [BS.State]
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
              bmRateMatrices    = [ BM.normalizedRateMatrix m stateFreqs popSize | m <- mutationModels ]
              bmStationaryDists = [ BM.stationaryDist m stateFreqs popSize | m <- mutationModels ]
              bmStationaryGens  = map Trans.stationaryDistToGenerator bmStationaryDists
              treesPrb          = [ Trans.branchLengthsToTransitionProbs b treeBM | b <- bmRateMatrices ]
              treesGen          = map Trans.treeProbMatrixToTreeGenerator treesPrb
  liftIO $ forM_ [1..nSites] simAndPrintOneSite
  -- Done.
  liftIO $ CF.close fileHandle
  liftIO $ putStrLn "Done."

main :: IO BMMParams
main = getParams >>= execStateT simulate
