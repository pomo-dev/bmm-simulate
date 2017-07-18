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

import           BndState              (Allele, BState (..), PopSize, connected,
                                        idToBState, stateSpaceSize)
import           DNAModel
import           Numeric.LinearAlgebra
import           RateMatrix
import           Tools
import           RTree

-- The boundary mutation models uses an underlying mutation model.
type MutModel         = RateMatrix
-- Let's just give rate matrices of the boundary mutation model a special name.
type BMModel          = RateMatrix

-- First we need to define a mutation model.
mutRate :: MutModel -> Allele -> Allele -> Double
mutRate m a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

stateFreq :: StateFreqVec -> Allele -> Double
stateFreq f a = f ! fromEnum a

moranCoef :: PopSize -> Int -> Double
moranCoef n i = iD * (nD - iD) / nD
  where iD = fromIntegral i
        nD = fromIntegral n

-- The transition rate from one state to another.
rate :: MutModel -> BState -> BState -> Double
rate m s t
  | not $ connected s t = 0.0
  | otherwise           = rate' s t
  where rate' (Ply n i _ _) _ = moranCoef n i
        rate' (Bnd _ a) (Ply _ _ b c)
          | a == b    = mutRate m a c
          | a == c    = mutRate m a b
          | otherwise = error "Cannot compute rate between states."
        rate' _ _ = error "Cannot compute rate between states."

-- The transition rate from one index to another.
rateById :: MutModel -> PopSize -> Int -> Int -> Double
rateById m n i j = rate m s t
  where s = idToBState n i
        t = idToBState n j

-- The build function (see below) has a weird way of assigning entries to
-- indices. The indices have to be the same data type as the entries. This is
-- just a helper function that changes the indices from Double to Int.
rateByDouble :: MutModel -> PopSize -> Double -> Double -> Double
rateByDouble m n x y = rateById m n (round x) (round y)

rateMatrix :: MutModel -> PopSize -> BMModel
rateMatrix m n = setDiagonal $ build (s,s) (rateByDouble m n)
  where s = stateSpaceSize n

-- Define a heterozygosity to make function definitions clearer.
type Heterozygosity = Double

-- Calculate the heterozygosity at stationarity.
theta :: MutModel -> StateFreqVec -> Heterozygosity
theta m f = f <.> rDiagZero #> f
  where e         = toExchMatrix m f
        (r, _)    = matrixSeparateSymSkew e
        -- The summation excludes the diagonal (a /= b).
        rDiagZero = matrixSetDiagToZero r

-- The normalization constant of the stationary distribution.
norm :: MutModel -> StateFreqVec -> PopSize -> Double
norm m f n = 1.0 + harmonic (n-1) * theta m f

-- Get entries of the stationary measure (not normalized) for a boundary model state.
stationaryMeasEntry :: MutModel -> StateFreqVec -> BState -> Double
stationaryMeasEntry _ f (Bnd _ a)     = piA
  where piA = f ! fromEnum a
stationaryMeasEntry m f (Ply n i a b) = piA * mAB / (fromIntegral n - fromIntegral i)
                                          + piB * mBA / fromIntegral i
  where piA = stateFreq f a
        mAB = mutRate m a b
        piB = stateFreq f b
        mBA = mutRate m b a

-- Only use this function when accessing single elements of the stationary
-- distribution. Otherwise, computation of the norm is unnecessarily repeated.
stationaryDistEntry :: MutModel -> StateFreqVec -> PopSize -> BState -> Double
stationaryDistEntry m f n s = stationaryMeasEntry m f s / norm m f n

stationaryMeasEntryById :: MutModel -> StateFreqVec -> PopSize -> Int -> Double
stationaryMeasEntryById m f n i = stationaryMeasEntry m f s
  where s = idToBState n i

stationaryMeasEntryByDouble :: MutModel -> StateFreqVec -> PopSize -> Double -> Double
stationaryMeasEntryByDouble m f n x = stationaryMeasEntryById m f n (round x)

stationaryDist :: MutModel -> StateFreqVec -> PopSize -> StationaryDist
stationaryDist m f n = scale (1.0/norm m f n) $
  build s (stationaryMeasEntryByDouble m f n)
  where s = stateSpaceSize n

-- Normalize the mutation coefficients such that the heterozygosity matches a
-- given level (see Eq. 12.14 in my thesis).
normalizeToTheta :: MutModel -> StateFreqVec -> PopSize -> Heterozygosity -> MutModel
normalizeToTheta m f n h = scale (h / (t * (1.0 - c * h))) m
  where
    -- The heterozygosity of the boundary mutation model.
    t = theta m f
    -- The branch length multiplicative factor introduced by the coalescent.
    c = harmonic (n-1)

-- The branch lengths of threes in the boundary mutation model are not measured
-- in average number of substitutions per site but in average number of
-- mutations or frequency shifts per site. The conversion factor is just the
-- square of the population size. This function converts the branch lengths of a
-- tree.
scaleTreeToBMM :: PopSize -> RTree a Double -> RTree a Double
scaleTreeToBMM n = fmap (* nSq)
  where nSq = fromIntegral n ** 2
