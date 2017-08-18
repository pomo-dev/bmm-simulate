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
import qualified Control.Monad.Random.Strict as Rand
import qualified DNAModel                    as DNA
import qualified RTree                       as Tree
import qualified System.Environment          as Sys
import qualified Transition                  as Trans

import qualified Numeric.LinearAlgebra as L

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
  -- This is super complicated. Is there an easier way?
  generator <-
    if seedArg == "random"
    then
      Rand.getStdGen
    else
    return $ Rand.mkStdGen (read seedArg :: Int)
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
  -- Tree.
  let treeHeight    = Args.treeHeight bmSimArgs
      treeSubs      = Tree.ilsTree treeHeight
      treeBM        = BM.scaleTreeToBMM popSize treeSubs
      treeTransProb = Trans.branchLengthsToTransitionProbs bmRateMatrix treeBM
  -- Other options.
  let nSites = Args.nSites bmSimArgs
  -- Simulation.
  let leafs        = Rand.evalRand transition generator
        where transition = Trans.simulateNSites nSites bmStationaryDist treeTransProb
      popNames     = map fst $ head leafs
      dataAllSites = map (map (BS.idToState popSize . snd)) leafs
      fileName     = Args.outFileName bmSimArgs
  -- TODO: Gathering all data in a single vector needs a LOT OF MEMORY. It is
  -- better, to write the file line by line (during the simulation). I.e., get
  -- rid of Trans.simulateNSites and provide functions that prepare the
  -- generators as well as functions that simulate and write one site at a time.
  CF.write fileName nSites popNames dataAllSites

  -- Output.
  putStrLn "Boundary mutation model simulator version 0.1.0.0."
  putStr "Command line: "
  putStrLn $ progName ++ " " ++ unwords args

  putStrLn ""
  putStrLn "--"
  putStrLn "General options."
  putStrLn $ "Seed: " ++ seedArg
  putStr "Generator: "
  print generator
  putStr "Number of simulated sites: "
  print nSites
  putStr "Output written to: "
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
