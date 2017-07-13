{- |
Module      :  DNA Markov models
Description :  An implementation of DNA Markov models (substitution models or mutation models)
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

This module provides rate matrices for the most common DNA substitution models
(they will be added according to need). They are also used as mutation models by
the boundary mutation model.

* Changelog

-}

module DNAModel
  where

import Numeric.LinearAlgebra
import RateMatrix

-- The nucleotides.
data Nuc = A | C | G | T deriving (Eq, Show, Read, Ord, Bounded, Enum)

type Kappa      = Double

-- The matrix of exchangeabilities in the HKY model.
exchangeabilityMatrixHKY :: Kappa -> Matrix R
exchangeabilityMatrixHKY k = (4><4)
  [ 0.0, 1.0,   k, 1.0
  , 1.0, 0.0, 1.0, k
  ,   k, 1.0, 0.0, 1.0
  , 1.0,   k, 1.0, 0.0 ]

-- HKY model mutation matrix normalized so that one mutation happens per unit time.
rateMatrixHKY :: StateFreqVec -> Kappa -> RateMatrix
rateMatrixHKY f k = rateMatrixNormalize f $ rateMatrixSetDiagonal $ exch <> diag f
  where exch = exchangeabilityMatrixHKY k
