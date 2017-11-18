{- |
Description :  DNA Markov models (substitution models or mutation models)
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
  , Nuc
  , getDNAModelInfoStr
  , rateMatrix
  , dnaModelSpecGetStateFreqVec )
  where

import           Data.List             (intercalate)
import           Numeric.LinearAlgebra
import           RateMatrix

-- | The nucleotides.
data Nuc = A | C | G | T deriving (Eq, Show, Read, Ord, Bounded, Enum)

nNucs :: Int
nNucs = 1 + fromEnum (maxBound :: Nuc)

-- | The DNA model specifications. E.g., for the HKY model, we only have kappa.
data DNAModelSpec = JC |
                    -- | HKY model with transition to transversion ratio kappa.
                    HKY Double StationaryDist |
                    -- | GTR model with five params and state frequency vector.
                    GTR Double Double Double Double Double StationaryDist

instance Show DNAModelSpec where
  show JC        = "JC"
  show (HKY k f) = "HKY[" ++ show k ++ "]" ++ show f
  show (GTR a b c d e f) = "GTR[" ++ paramsStr ++ "]" ++ show f
    where paramsStr = intercalate "," (map show [a,b,c,d,e])

-- | Extract stationary frequency vector of a DNA model specification 'DNAModelSpec'.
dnaModelSpecGetStateFreqVec :: DNAModelSpec -> StationaryDist
dnaModelSpecGetStateFreqVec JC                = vector $ replicate nNucs 0.25
dnaModelSpecGetStateFreqVec (HKY _ f)         = f
dnaModelSpecGetStateFreqVec (GTR _ _ _ _ _ f) = f

-- | A rate matrix of a DNA models is called DNAModel.
data DNAModel = DNAModel { dnaRateMatrix :: RateMatrix
                         , dnaModelSpec  :: DNAModelSpec }

-- | The matrix of exchangeabilities; diagonal entries are not set.
exchangeabilityMatrix :: DNAModelSpec -> Matrix R
exchangeabilityMatrix JC = (4><4)
  [ 0.0, 1.0, 1.0, 1.0
  , 1.0, 0.0, 1.0, 1.0
  , 1.0, 1.0, 0.0, 1.0
  , 1.0, 1.0, 1.0, 0.0 ]
exchangeabilityMatrix (HKY k _) = (4><4)
  [ 0.0, 1.0,   k, 1.0
  , 1.0, 0.0, 1.0,   k
  ,   k, 1.0, 0.0, 1.0
  , 1.0,   k, 1.0, 0.0 ]
exchangeabilityMatrix (GTR a b c d e _) = (4><4)
  [ 0.0,   a,   b,   c
  ,   a, 0.0,   d,   e
  ,   b,   d, 0.0, 1.0
  ,   c,   e, 1.0, 0.0 ]
-- exchangeabilityMatrix _     = error "Model not yet supported."

-- | DNA model rate matrix normalized so that one mutation happens per unit time.
rateMatrix :: DNAModelSpec -> DNAModel
rateMatrix s = DNAModel rm s
  where
    f    = dnaModelSpecGetStateFreqVec s
    exch = exchangeabilityMatrix s
    rm   = normalizeRates f $ setDiagonal $ exch <> diag f
-- rateMatrix _ _       = error "Model not yet supported."

-- | Report a specific model.
getDNAModelInfoStr :: DNAModel -> String
getDNAModelInfoStr (DNAModel m s) =
  case s of
    JC        -> unlines
      [ "JC model with rate matrix"
      , show m ]
    (HKY k f) -> unlines
      [ "HKY model"
      , reportRateMatrix
      , reportExchangeabilityMatrix
      , reportStateFreqVec f
      , "And a kappa value of: " ++ show k ]
    (GTR a b c d e f) -> unlines
      [ "GTR model"
      , reportRateMatrix
      , reportExchangeabilityMatrix
      , reportStateFreqVec f
      , "And parameters: " ++ intercalate ", " (map show [a,b,c,d,e]) ]
      -- _ -> error "Model not yet supported."
  where reportRateMatrix   = "Rate matrix " ++ show m
        reportStateFreqVec f = "This corresponds to state frequencies (A, C, G, T): " ++ show f
        reportExchangeabilityMatrix = "Exchangeability matrix (diagonal set to 0.0) " ++
                                      show (exchangeabilityMatrix s)
