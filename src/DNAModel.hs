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

import           Numeric.LinearAlgebra
import           RateMatrix

-- The nucleotides.
data Nuc = A | C | G | T deriving (Eq, Show, Read, Ord, Bounded, Enum)

-- The DNA model specifications. E.g., for the HKY model, we only have kappa.
newtype DNAModelSpec = HKY Double -- HKY model with transition to transversion
                                  -- ratio kappa.

-- TODO: Is too much work at the moment. Because log file parser also has to be changed.
-- instance Show DNAModelSpec where
--   show (HKY k) = "HKY substitution model"

-- A rate matrix of a DNA models is called DNAModel.
data DNAModel = DNAModel { dnaRateMatrix  :: RateMatrix
                         , dnaModelParams :: DNAModelSpec }

-- A stationary distribution of a DNA model is also called stationary frequency
-- vector.
type StateFreqVec  = StationaryDist

-- The matrix of exchangeabilities.
exchangeabilityMatrix :: DNAModelSpec -> Matrix R
exchangeabilityMatrix (HKY k) = (4><4)
  [ 0.0, 1.0,   k, 1.0
  , 1.0, 0.0, 1.0, k
  ,   k, 1.0, 0.0, 1.0
  , 1.0,   k, 1.0, 0.0 ]
-- exchangeabilityMatrix _     = error "Model not yet supported."

-- HKY model mutation matrix normalized so that one mutation happens per unit time.
rateMatrix :: StateFreqVec -> DNAModelSpec -> DNAModel
rateMatrix f s@(HKY _) = DNAModel rm s
  where exch = exchangeabilityMatrix s
        rm   = normalizeRates f $ setDiagonal $ exch <> diag f
-- rateMatrix _ _       = error "Model not yet supported."

getDNAModelInfoStr :: DNAModel -> StateFreqVec -> String
getDNAModelInfoStr (DNAModel m s) f = case s of
       (HKY k) -> unlines $
         ["HKY model."] ++
         reportMatrix ++
         reportStateFreqVec ++
         [ "And a kappa value of: " ++ show k ]
       -- _ -> error "Model not yet supported."
       where reportMatrix = [ "Rate matrix:", show m]
             reportStateFreqVec = [ "This corresponds to state frequencies (A, C, G, T): " ++ show f ]
