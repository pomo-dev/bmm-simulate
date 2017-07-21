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

import qualified BndModel                  as BM
import qualified BndState                  as BS
import qualified CFWriter                  as CF
import           Control.Monad.Random.Lazy
import qualified DNAModel                  as DNA
import qualified ExampleRTrees             as Trees
import           Numeric.LinearAlgebra
import qualified Transition                as Trans

-- TODO: Read in the simulation parameters such as tree type, mutation model
-- with parameters, stationary distribution, heterozygosity, population size.

-- TODO: Properly output simulation parameters and command line and the seed
-- (everything needed to repeat simulations and to access all relevant
-- parameters).

-- TODO: Testing (Quicktest?).

main :: IO ()
main = do
  print "Boundary mutation model simulator."
  let stateFreqs      = vector [0.3, 0.2, 0.2, 0.3]
      popSize         = 9
      -- The heterozygosity value.
      heterozygosity  = 0.0025
      -- A kappa value of 6.25 corresponds to a transition to transversion ratio of 3.0
      mutationModel   = BM.normalizeToTheta hkyModel stateFreqs popSize heterozygosity
        where hkyModel= DNA.rateMatrixHKY stateFreqs 6.25
      rateMatrix      = BM.rateMatrix mutationModel popSize
      treeBrLn        = BM.scaleTreeToBMM popSize Trees.ilsTree1Ne
      treeTransProb   = Trans.branchLengthsToTransitionProbs rateMatrix treeBrLn
      stationaryDist  = BM.stationaryDist mutationModel stateFreqs popSize
      nSites          = 100000
  print "The mutation model matrix is:"
  print mutationModel
  leafs <- evalRandIO (Trans.simulateNSites nSites stationaryDist treeTransProb)
  let popNames     = map fst $ head leafs
      dataAllSites = map (map (BS.idToBState popSize . snd)) leafs
      fileName     = "Test.cf"
  CF.write fileName nSites popNames dataAllSites
