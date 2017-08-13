{- |
Module      :  Defaults
Description :  Default values and constants
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Default values and constants.

* Changelog

-}

module Defaults where

import qualified DNAModel              as DNA
import qualified Numeric.LinearAlgebra as LinAlg
import qualified System.Random         as Rand

-- Stationary distribution of the mutation model (or stationary frequencies).
stateFreqs :: DNA.StateFreqVec
stateFreqs = LinAlg.vector [0.3, 0.2, 0.2, 0.3]

-- Virtual population size.
popSize :: Int
popSize = 9

-- A kappa value of 6.25 corresponds to a transition to transversion ratio of
-- 3.0
kappa :: Double
kappa = 6.25

-- Heterozygosity value.
heterozygosity :: Double
heterozygosity = 0.0025

-- Tree height.
treeHeight :: Double
treeHeight = 0.005

-- Number of sites to simulate.
nSites :: Int
nSites = 1000000

-- By default, the seed is random.
seed :: IO Rand.StdGen
seed = Rand.getStdGen