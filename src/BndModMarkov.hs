{- |
Module      :  Boundary model transition rate matrix
Description :  An implementation of the boundary model Markov process
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

This module defines the transition rate matrix of the boundary mutation model.

* Changelog

-}

module BndModMarkov
  where

import BndModState ( Allele
                   , PopSize
                   , BState(..)
                   , Bnd(..)
                   , Ply(..)
                   , stateSpaceSize
                   , indToBState
                   , connected)
import Numeric.LinearAlgebra

type StateFreqVec = Vector R
type RateMatrix   = Matrix R
type MutModel     = Matrix R

-- First we need to define a mutation model.
mutRate :: RateMatrix -> Allele -> Allele -> Double
mutRate m a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

matrixSetDiagToZero :: Matrix R -> Matrix R
matrixSetDiagToZero m = m - diag (takeDiag m)

-- Normalizes a Markov process generator such that one event happens per unit time.
rateMatrixNormalize :: StateFreqVec -> RateMatrix -> RateMatrix
rateMatrixNormalize f m = scale (1.0 / totalRate) m
  where totalRate = norm_1 $ f <# matrixSetDiagToZero m

-- Set the diagonal entries of a matrix such that the rows sum to 0.
rateMatrixSetDiagonal :: RateMatrix -> RateMatrix
rateMatrixSetDiagonal m = diagZeroes - diag (fromList rowSums)
  where diagZeroes = matrixSetDiagToZero m
        rowSums    = map norm_1 $ toRows diagZeroes

type Kappa = Double

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

-- The transition rate from one state to another.
bmRate :: MutModel -> BState -> BState -> Double
bmRate m s t
  | not $ connected s t = 0.0
  | otherwise           = rate s t
  where rate (BStatePly (Ply n i _ _)) _ = fromIntegral i * (fromIntegral n - fromIntegral i) / fromIntegral n
        rate (BStateBnd (Bnd _ a)) (BStatePly (Ply _ _ b c))
          | a == b    = mutRate m a c
          | a == c    = mutRate m a b
          | otherwise = error "Cannot compute rate between states."
        rate _ _ = error "Cannot compute rate between states."

-- The transition rate from one index to another.
bmRateByIndex :: MutModel -> PopSize -> Int -> Int -> Double
bmRateByIndex m n i j = bmRate m s t
  where s = indToBState n i
        t = indToBState n j

-- The build function (see below) has a weird way of assigning entries to
-- indices. The indices have to be the same data type as the entries. This is
-- just a helper function that changes the indices from Double to Int.
bmRateByDouble :: MutModel -> PopSize -> Double -> Double -> Double
bmRateByDouble m n x y = bmRateByIndex m n (round x) (round y)

rateMatrixBM :: MutModel -> PopSize -> RateMatrix
rateMatrixBM m n = rateMatrixSetDiagonal $ build (s,s) (bmRateByDouble m n)
  where s = stateSpaceSize n
