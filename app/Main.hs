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
  generator <- case seedArg of
      "random" -> Rand.getStdGen
               -- This is a little tricky because it is hard to parse two
               -- integers on the command line (space messes things up and
               -- workarounds are not user friendly). Like this, the second seed
               -- value is always set to one (it turns out that getStdGen also
               -- always sets the second seed to one).
      _        -> return (read $ seedArg ++ " 1" :: Rand.StdGen)
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
  let nSites = Args.nSites bmSimArgs
  -- Output file.
  let fileName = Args.outFileName bmSimArgs
  fileHandle <- CF.open fileName
  CF.writeHeader fileHandle nSites popNames
  -- Simulation.

  -- Set chromosome name to "SIM".
  let leafs  = Rand.evalRand (Trans.simulateSite bmStationaryGen treeGenerator) generator
      states = map (BS.idToState popSize . snd) leafs
  CF.writeLine fileHandle "SIM" 1 states

  -- let leafs = Rand.evalRand transition generator
  --       where transition = Trans.simulateNSites nSites bmStationaryDist treeTransProb
  --     dataAllSites = map (map (BS.idToState popSize . snd)) leafs

  CF.close fileHandle

  -- TODO: Eventually move output before simulation.
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
