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
import           Control.Monad.Random.Lazy
import qualified DNAModel                  as DNA
import qualified ExampleRTrees             as Trees
import           Numeric.LinearAlgebra
import qualified Transition                as Trans

-- TODO: Read in the simulation parameters such as tree type, mutation model
-- with parameters, stationary distribution, heterozygosity, population size.
main :: IO ()
main = do
  let stateFreqs      = vector [0.3, 0.2, 0.2, 0.3]
      popSize         = 10
      -- The heterozygosity value.
      heterozygosity  = 0.0025
      -- A kappa value of 6.25 corresponds to a transition to transversion ratio of 3.0
      mutationModel   = BM.normalizeToTheta hkyModel stateFreqs popSize heterozygosity
        where hkyModel= DNA.rateMatrixHKY stateFreqs 6.25
      rateMatrix      = BM.rateMatrix mutationModel popSize
      treeBrLn        = BM.scaleTreeToBMM popSize Trees.ilsTree1Ne
      treeTransProb   = Trans.branchLengthsToTransitionProbs rateMatrix treeBrLn
      stationaryDist  = BM.stationaryDist mutationModel stateFreqs popSize
  leafs <- evalRandIO (Trans.simulateNSites 1000 stationaryDist treeTransProb)
  print leafs
  return ()
