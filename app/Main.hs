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
  let stateFreqs      = Args.stateFreqs bmSimArgs
      popSize         = Args.popSize bmSimArgs
      kappa           = Args.kappa bmSimArgs
      heterozygosity  = Args.heterozygosity bmSimArgs
      treeHeight      = Args.treeHeight bmSimArgs
      mutationModel   = BM.normalizeToTheta
        hkyModel stateFreqs popSize heterozygosity
        where hkyModel = DNA.rateMatrixHKY stateFreqs kappa
      rateMatrix      = BM.rateMatrix mutationModel popSize
      treeSubs        = Tree.ilsTree treeHeight
      treeBM          = BM.scaleTreeToBMM popSize treeSubs
      treeTransProb   = Trans.branchLengthsToTransitionProbs rateMatrix treeBM
      stationaryDist  = BM.stationaryDist mutationModel stateFreqs popSize
      nSites          = Args.nSites bmSimArgs
      leafs           = Rand.evalRand transition generator
        where transition = Trans.simulateNSites nSites stationaryDist treeTransProb
      popNames     = map fst $ head leafs
      dataAllSites = map (map (BS.idToState popSize . snd)) leafs
      fileName     = Args.outFileName bmSimArgs
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
