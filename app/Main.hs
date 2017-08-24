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

import qualified ArgParse                    as Args
import qualified BndModel                    as BM
import qualified BndState                    as BS
import qualified CFWriter                    as CF
import           Control.Monad.Random.Strict
import qualified DNAModel                    as DNA
import qualified RateMatrix                  as RM
import qualified RTree                       as Tree
import qualified System.Environment          as Sys
import qualified Transition                  as Trans

-- import qualified Numeric.LinearAlgebra as L

-- Automatic version information does not work with flycheck ... ahhhh.
-- However, intero does not provide show warnings and so on.
-- import           Paths_bmm_simulate        (version)
-- import           Data.Version              (showVersion)

-- TODO: Read in tree type, specific mutation model.

main :: IO ()
main = do
  progName <- Sys.getProgName
  args <- Sys.getArgs
  bmSimArgs <- Args.parseBMSimArgs
  let seedArg = Args.seed bmSimArgs
  unless (seedArg == "random") $ setStdGen (read $ seedArg ++ " 1")
  generator <- getStdGen

  -- Mutation model.
  let stateFreqs = Args.stateFreqs bmSimArgs
      kappa      = Args.kappa bmSimArgs
      hkyModel   = DNA.rateMatrixHKY stateFreqs kappa
  -- Boundary mutation model.
  let popSize          = Args.popSize bmSimArgs
      heterozygosity   = Args.heterozygosity bmSimArgs
      mutationModel    = BM.normalizeToTheta
        hkyModel stateFreqs popSize heterozygosity
      bmRateMatrix     = BM.normalizedRateMatrix mutationModel stateFreqs popSize
      bmStationaryDist = BM.stationaryDist mutationModel stateFreqs popSize
      bmStationaryGen  = Trans.stationaryDistToGenerator bmStationaryDist
  -- Tree.
  let treeHeight    = Args.treeHeight bmSimArgs
      treeSubs      = Tree.ilsTree treeHeight
      treeBM        = BM.scaleTreeToBMM popSize treeSubs
      popNames      = Tree.getLeaves treeBM
      treeTransProb = Trans.branchLengthsToTransitionProbs bmRateMatrix treeBM
      treeGenerator = Trans.treeProbMatrixToTreeGenerator treeTransProb
  -- Other options.
  let nSites   = Args.nSites bmSimArgs
      fileName = Args.outFileName bmSimArgs

  -- Output.
  putStrLn "Boundary mutation model simulator version 0.1.0.0."
  putStr "Command line: "
  putStrLn $ progName ++ " " ++ unwords args

  putStrLn ""
  putStrLn "--"
  putStrLn "General options."
  putStr $ "Seed: " ++ seedArg
  putStrLn $ " (Generator: " ++ show generator ++ ")"
  putStr "Number of simulated sites: "
  print nSites
  putStr "Output will be written to: "
  print fileName

  putStrLn ""
  putStrLn "--"
  putStrLn "Boundary mutation model options."
  putStr "Population size: "
  print popSize
  putStr "Heterozygosity: "
  print heterozygosity
  putStrLn "Mutation model matrix:"
  print mutationModel
  putStr "This corresponds to state frequencies (A, C, G, T): "
  print stateFreqs
  putStr "And a kappa value of: "
  print kappa

  putStrLn ""
  putStrLn "--"
  putStrLn "Tree options."
  putStr "Species tree in average number of substitutions: "
  putStrLn $ Tree.toNewick treeSubs
  putStr "Species tree in  mutations and frequency shifts: "
  putStrLn $ Tree.toNewick treeBM

  putStrLn ""
  putStrLn "--"
  putStrLn "Performing simulation."
  -- Prepare output file.
  fileHandle <- CF.open fileName
  CF.writeHeader fileHandle nSites popNames
  -- Simulation helpers.
  let toStates = map (BS.idToState popSize . snd) :: [(a, RM.State)] -> [BS.State]
      writer   = CF.writeLine fileHandle "SIM"    :: CF.Pos -> CF.DataOneSite -> IO ()
  -- Simulation; loop over positions.
  forM_ [1..nSites] $ \pos -> evalRandIO (Trans.simulateSite bmStationaryGen treeGenerator)
                              >>= writer pos . toStates
  -- Done.
  CF.close fileHandle
  putStrLn "Done."
