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

module BndModel
  where

import Numeric.LinearAlgebra
import BndState ( Allele
                , PopSize
                , BState(..)
                , stateSpaceSize
                , indToBState
                , connected)
import RateMatrix

-- The boundary mutation models uses an underlying mutation model.
type MutModel   = RateMatrix

-- First we need to define a mutation model.
mutRate :: MutModel -> Allele -> Allele -> Double
mutRate m a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

moranRate :: PopSize -> Int -> Double
moranRate n i = fromIntegral i * (fromIntegral n - fromIntegral i) / fromIntegral n

-- The transition rate from one state to another.
bmRate :: MutModel -> BState -> BState -> Double
bmRate m s t
  | not $ connected s t = 0.0
  | otherwise           = rate s t
  where rate (Ply n i _ _) _ = moranRate n i
        rate (Bnd _ a) (Ply _ _ b c)
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

type BranchLength = Double
type ProbMatrix   = Matrix R

-- The important matrix that gives the probabilities to move from one state to
-- another in a specific time (branch length).
probMatrix :: RateMatrix -> BranchLength -> ProbMatrix
probMatrix m t = expm $ scale t m
