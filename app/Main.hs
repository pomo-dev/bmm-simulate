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

import qualified ArgParse                  as Args
import qualified BndModel                  as BM
import qualified BndState                  as BS
import qualified CFWriter                  as CF
import qualified Control.Monad.Random.Lazy as Rand
import qualified DNAModel                  as DNA
import qualified ExampleRTrees             as Trees
import qualified Numeric.LinearAlgebra     as LinAlg
import qualified RTree                     as R
import qualified Transition                as Trans

-- Automatic version information does not work with flycheck ... ahhhh.
-- However, intero does not provide show warnings and so on.
-- import           Paths_bmm_simulate        (version)
-- import           Data.Version              (showVersion)

-- TODO: Read in the simulation parameters such as tree type, mutation model
-- with parameters, stationary distribution, heterozygosity, population size.

-- TODO: Allow manual seed.

main :: IO ()
main = do
  bmSimArgs <- Args.parseBMSimArgs
  seed <- Rand.getStdGen
  let -- Specification of the boundary mutation model.
      stateFreqs      = Args.stateFreqs bmSimArgs
      popSize         = 9
      -- A kappa value of 6.25 corresponds to a transition to transversion ratio of 3.0
      kappa           = 6.25
      -- The heterozygosity value.
      heterozygosity  = 0.0025
      -- The tree height.
      treeHeight      = 0.005
      mutationModel   = BM.normalizeToTheta hkyModel stateFreqs popSize heterozygosity
        where hkyModel= DNA.rateMatrixHKY stateFreqs kappa
      rateMatrix      = BM.rateMatrix mutationModel popSize
      treeSubs        = Trees.ilsTree treeHeight
      treeBM          = BM.scaleTreeToBMM popSize treeSubs
      treeTransProb   = Trans.branchLengthsToTransitionProbs rateMatrix treeBM
      stationaryDist  = BM.stationaryDist mutationModel stateFreqs popSize
      nSites          = 10000
  let leafs = Rand.evalRand (Trans.simulateNSites nSites stationaryDist treeTransProb) seed
  let popNames     = map fst $ head leafs
      dataAllSites = map (map (BS.idToBState popSize . snd)) leafs
      fileName     = Args.outFileName bmSimArgs
  CF.write fileName nSites popNames dataAllSites

  -- Output.
  putStrLn "Boundary mutation model simulator."
  putStrLn "--"
  putStr "Seed: "
  print seed
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
  putStr "Species tree in average number of substitutions: "
  putStrLn $ R.toNewick treeSubs
  putStr "Species tree in  mutations and frequency shifts: "
  putStrLn $ R.toNewick treeBM
  putStr "Number of simulated sites: "
  print nSites
  putStr "Output written to: "
  print fileName
