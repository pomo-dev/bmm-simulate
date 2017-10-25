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
  ( DNAModelSpec(..)
  , DNAModel(..)
  , StateFreqVec
  , Nuc
  , getDNAModelInfoStr
  , rateMatrix
  , dnaModelSpecGetStateFreqVec )
  where

import           Numeric.LinearAlgebra
import           RateMatrix

-- The nucleotides.
data Nuc = A | C | G | T deriving (Eq, Show, Read, Ord, Bounded, Enum)

nNucs :: Int
nNucs = 1 + fromEnum (maxBound :: Nuc)

-- The DNA model specifications. E.g., for the HKY model, we only have kappa.
data DNAModelSpec = JC |
                    -- HKY model with transition to transversion ratio kappa.
                    HKY Double StateFreqVec

instance Show DNAModelSpec where
  show (HKY k f) = "HKY[" ++ show k ++ "]" ++ show f
  show JC        = "JC"

dnaModelSpecGetStateFreqVec :: DNAModelSpec -> StateFreqVec
dnaModelSpecGetStateFreqVec (HKY _ f) = f
dnaModelSpecGetStateFreqVec JC = vector $ replicate nNucs 0.25

-- A rate matrix of a DNA models is called DNAModel.
data DNAModel = DNAModel { dnaRateMatrix :: RateMatrix
                         , dnaModelSpec  :: DNAModelSpec }

-- A stationary distribution of a DNA model is also called stationary frequency
-- vector.
type StateFreqVec  = StationaryDist

-- The matrix of exchangeabilities.
exchangeabilityMatrix :: DNAModelSpec -> Matrix R
exchangeabilityMatrix (HKY k _) = (4><4)
  [ 0.0, 1.0,   k, 1.0
  , 1.0, 0.0, 1.0, k
  ,   k, 1.0, 0.0, 1.0
  , 1.0,   k, 1.0, 0.0 ]
exchangeabilityMatrix JC = (4><4)
  [ 0.0, 1.0, 1.0, 1.0
  , 1.0, 0.0, 1.0, 1.0
  , 1.0, 1.0, 0.0, 1.0
  , 1.0, 1.0, 1.0, 0.0 ]
-- exchangeabilityMatrix _     = error "Model not yet supported."

-- | DNA model rate matrix normalized so that one mutation happens per unit time.
rateMatrix :: DNAModelSpec -> DNAModel
rateMatrix s =
  case s of
    (HKY _ _) -> DNAModel rm s
    JC        -> DNAModel rm s
    -- _ -> error "Model not yet supported."
  where
    f    = dnaModelSpecGetStateFreqVec s
    exch = exchangeabilityMatrix s
    rm   = normalizeRates f $ setDiagonal $ exch <> diag f
-- rateMatrix _ _       = error "Model not yet supported."

getDNAModelInfoStr :: DNAModel -> String
getDNAModelInfoStr (DNAModel m s) =
  case s of
    (HKY k f) -> unlines
      [ "HKY model with rate matrix"
      , show m
      , reportStateFreqVec f
      , "And a kappa value of: " ++ show k ]
    JC        -> unlines
      [ "JC model with rate matrix"
      , show m ]
      -- _ -> error "Model not yet supported."
  where reportStateFreqVec f = "This corresponds to state frequencies (A, C, G, T): " ++ show f
