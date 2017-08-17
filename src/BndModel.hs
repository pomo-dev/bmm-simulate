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

import qualified BndState              as BS
import qualified DNAModel              as DNA
import           Numeric.LinearAlgebra
import qualified RateMatrix            as RM
import qualified RTree                 as Tree
import qualified Tools                 as T

-- The boundary mutation models uses an underlying mutation model.
type MutModel         = RM.RateMatrix
-- Let's just give rate matrices of the boundary mutation model a special name.
type BndMutModel      = RM.RateMatrix

-- First we need to define a mutation model.
mutRate :: MutModel -> BS.Allele -> BS.Allele -> Double
mutRate m a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

stateFreq :: DNA.StateFreqVec -> BS.Allele -> Double
stateFreq f a = f ! fromEnum a

moranCoef :: BS.PopSize -> Int -> Double
moranCoef n i = iD * (nD - iD) / nD
  where iD = fromIntegral i
        nD = fromIntegral n

-- The transition rate from one boundary state to another.
rate :: MutModel -> BS.State -> BS.State -> Double
rate m s t
  | not $ BS.connected s t = 0.0
  | otherwise           = rate' s t
  where rate' (BS.Ply n i _ _) _ = moranCoef n i
        rate' (BS.Bnd _ a) (BS.Ply _ _ b c)
          | a == b    = mutRate m a c
          | a == c    = mutRate m a b
          | otherwise = error "Cannot compute rate between states."
        rate' _ _ = error "Cannot compute rate between states."

-- The transition rate from one state (index) to another.
rateById :: MutModel -> BS.PopSize -> RM.State -> RM.State -> Double
rateById m n i j = rate m s t
  where s = BS.idToState n i
        t = BS.idToState n j

-- The build function (see below) has a weird way of assigning entries to
-- indices. The indices have to be the same data type as the entries. This is
-- just a helper function that changes the indices from Double to Int.
rateByDouble :: MutModel -> BS.PopSize -> Double -> Double -> Double
rateByDouble m n x y = rateById m n (round x) (round y)

rateMatrix :: MutModel -> BS.PopSize -> BndMutModel
rateMatrix m n = RM.setDiagonal $ build (s,s) (rateByDouble m n)
  where s = BS.stateSpaceSize n

-- Define a heterozygosity to make function definitions clearer.
type Heterozygosity = Double

-- Calculate the heterozygosity at stationarity.
theta :: MutModel -> DNA.StateFreqVec -> Heterozygosity
theta m f = f <.> rDiagZero #> f
  where e         = RM.toExchMatrix m f
        (r, _)    = T.matrixSeparateSymSkew e
        -- The summation excludes the diagonal (a /= b).
        rDiagZero = T.matrixSetDiagToZero r

-- The normalization constant of the stationary distribution.
norm :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> Double
norm m f n = 1.0 + T.harmonic (n-1) * theta m f

-- Get entries of the stationary measure (not normalized) for a boundary model state.
stationaryMeasEntry :: MutModel -> DNA.StateFreqVec -> BS.State -> Double
stationaryMeasEntry _ f (BS.Bnd _ a)     = piA
  where piA = f ! fromEnum a
stationaryMeasEntry m f (BS.Ply n i a b) = piA * mAB / (fromIntegral n - fromIntegral i)
                                           + piB * mBA / fromIntegral i
  where piA = stateFreq f a
        mAB = mutRate m a b
        piB = stateFreq f b
        mBA = mutRate m b a

-- Only use this function when accessing single elements of the stationary
-- distribution. Otherwise, computation of the norm is unnecessarily repeated.
stationaryDistEntry :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> BS.State -> Double
stationaryDistEntry m f n s = stationaryMeasEntry m f s / norm m f n

stationaryMeasEntryById :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> Int -> Double
stationaryMeasEntryById m f n i = stationaryMeasEntry m f s
  where s = BS.idToState n i

stationaryMeasEntryByDouble :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> Double -> Double
stationaryMeasEntryByDouble m f n x = stationaryMeasEntryById m f n (round x)

stationaryDist :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> RM.StationaryDist
stationaryDist m f n = scale (1.0/norm m f n) $
  build s (stationaryMeasEntryByDouble m f n)
  where s = BS.stateSpaceSize n

-- Normalize the mutation coefficients such that the heterozygosity matches a
-- given level (see Eq. 12.14 in my thesis).
normalizeToTheta :: MutModel -> DNA.StateFreqVec -> BS.PopSize -> Heterozygosity -> MutModel
normalizeToTheta m f n h = scale (h / (t * (1.0 - c * h))) m
  where
    -- The heterozygosity of the boundary mutation model.
    t = theta m f
    -- The branch length multiplicative factor introduced by the coalescent.
    c = T.harmonic (n-1)

-- The branch lengths of threes in the boundary mutation model are not measured
-- in average number of substitutions per site but in average number of
-- mutations or frequency shifts per site. The conversion factor is just the
-- square of the population size. This function converts the branch lengths of a
-- tree.
scaleTreeToBMM :: BS.PopSize -> Tree.RTree a Double -> Tree.RTree a Double
scaleTreeToBMM n = fmap (* nSq)
  where nSq = fromIntegral n ** 2
