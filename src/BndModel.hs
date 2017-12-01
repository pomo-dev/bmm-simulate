{- |
Description :  Boundary mutation model Markov process
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

This module defines the transition rate matrix of the boundary mutation model.

-}

module BndModel
  ( BMModel(..)
  , MutModel
  , Heterozygosity
  , createBMM
  , createBMMNormalized
  , scaleTreeToBMM
  , getBMMInfoStr
  ) where

import           BndState
import           Data.Maybe
import           DNAModel              hiding (rateMatrix)
import           Numeric.LinearAlgebra
import           RateMatrix            hiding (State)
import           RTree
import           Tools

-- | The boundary mutation models uses an underlying mutation model.
type MutModel = DNAModel

-- | Define a heterozygosity to make function definitions clearer.
type Heterozygosity = Double

-- | A boundary mutation model is defined by its rate matrix. However, it is
-- convenient to keep the underlying mutation model, the population size the
-- heterozygosity and the stationary distribution at hand.
data BMModel = BMModel { bmmRateMatrix     :: RateMatrix
                       , bmmMutModel       :: MutModel
                       , bmmPopSize        :: PopSize
                       , bmmHeterozygosity :: Heterozygosity
                       , bmmStationaryDist :: StationaryDist }

-- Get the mutation rate from one allele to another.
mutRate :: MutModel -> Allele -> Allele -> Double
mutRate (DNAModel m _) a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

-- Get state frequency for specific allele.
stationaryFreq :: StationaryDist -> Allele -> Double
stationaryFreq f a = f ! fromEnum a

-- Calculate the rate of frequency shifts.
moranCoef :: PopSize -> Int -> Double
moranCoef n i = iD * (nD - iD) / nD
  where iD = fromIntegral i
        nD = fromIntegral n

-- The transition rate from one boundary state to another.
rate :: MutModel -> State -> State -> Double
rate m s t
  | not $ connected s t = 0.0
  | otherwise           = rate' s t
  where rate' (Ply n i _ _) _ = moranCoef n i
        rate' (Bnd _ a) (Ply _ _ b c)
          | a == b    = mutRate m a c
          | a == c    = mutRate m a b
          | otherwise = error "Cannot compute rate between states."
        rate' _ _ = error "Cannot compute rate between states."

-- The transition rate from one state (index) to another.
rateById :: MutModel -> PopSize -> Int -> Int -> Double
rateById m n i j = rate m s t
  where s = idToState n i
        t = idToState n j

-- The build function (see below) has a weird way of assigning entries to
-- indices. The indices have to be the same data type as the entries. This is
-- just a helper function that changes the indices from Double to Int.
rateByDouble :: MutModel -> PopSize -> Double -> Double -> Double
rateByDouble m n x y = rateById m n (round x) (round y)

-- The rate matrix of the boundary mutation model. The dimension with four
-- alleles will be 4 + 6*(N-1).
rateMatrix :: MutModel -> PopSize -> RateMatrix
rateMatrix m n = setDiagonal $ build (s,s) (rateByDouble m n)
  where s = stateSpaceSize n

-- Normalize the rate matrix such that on average one event happens per unit
-- time.
normalizedRateMatrix :: MutModel -> PopSize -> RateMatrix
normalizedRateMatrix m n = normalizeRates f' m'
  where m' = rateMatrix m n
        -- f  = dnaModelSpecGetStateFreqVec (dnaModelSpec m)
        f' = stationaryDist m n

-- Calculate the heterozygosity at stationarity.
theta :: MutModel -> Heterozygosity
theta (DNAModel m s) = f <.> rDiagZero #> f
  where f         = dnaModelSpecGetStateFreqVec s
        e         = toExchMatrix m f
        (r, _)    = matrixSeparateSymSkew e
        -- The summation excludes the diagonal (a /= b).
        rDiagZero = matrixSetDiagToZero r

-- The normalization constant of the stationary distribution.
norm :: MutModel -> PopSize -> Double
norm m n = 1.0 + harmonic (n-1) * theta m

-- Get entries of the stationary measure (not normalized) for a boundary model state.
stationaryMeasEntry :: MutModel -> State -> Double
stationaryMeasEntry m s =
  let f   = dnaModelSpecGetStateFreqVec (dnaModelSpec m) in
  case s of
    (Bnd _ a)     -> f ! fromEnum a
    (Ply n i a b) -> piA * mAB / (fromIntegral n - fromIntegral i)
                                      + piB * mBA / fromIntegral i
      where piA = stationaryFreq f a
            mAB = mutRate m a b
            piB = stationaryFreq f b
            mBA = mutRate m b a

stationaryMeasEntryById :: MutModel -> PopSize -> Int -> Double
stationaryMeasEntryById m n i = stationaryMeasEntry m s
  where s = idToState n i

stationaryMeasEntryByDouble :: MutModel -> PopSize -> Double -> Double
stationaryMeasEntryByDouble m n x = stationaryMeasEntryById m n (round x)

-- Get the stationary distribution of a boundary mutation model.
stationaryDist :: MutModel -> PopSize -> StationaryDist
stationaryDist m n = scale (1.0/norm m n) $
  build s (stationaryMeasEntryByDouble m n)
  where s = stateSpaceSize n

-- Normalize the mutation coefficients such that the heterozygosity matches a
-- given level (see Eq. 12.14 in my thesis).
normalizeToTheta :: MutModel -> PopSize -> Heterozygosity -> MutModel
normalizeToTheta mo@(DNAModel m _) n h =
  mo { dnaRateMatrix = scale (h / (t * (1.0 - c * h))) m }
  where
    -- The heterozygosity of the boundary mutation model.
    t = theta mo
    -- The branch length multiplicative factor introduced by the coalescent.
    c = harmonic (n-1)

-- | Create a boundary mutation model using the minimal number of necessary
-- ingredients.
createBMM :: DNAModel -> PopSize -> Heterozygosity -> BMModel
createBMM m n h = BMModel rm m' n h f
  where m' = normalizeToTheta m n h
        rm = normalizedRateMatrix m' n
        f  = stationaryDist m n

-- | Create a boundary mutation model without providing a heterozygosity. This
-- means that the mutation model has to be normalized already.
createBMMNormalized :: DNAModel -> PopSize -> BMModel
createBMMNormalized m n = BMModel rm m n h f
  where rm = normalizedRateMatrix m n
        h  = theta m
        f  = stationaryDist m n

-- | The branch lengths of threes in the boundary mutation model are not
-- measured in average number of substitutions per site but in average number of
-- mutations or frequency shifts per site. The conversion factor is just the
-- square of the population size. This function converts the branch lengths of a
-- tree.
scaleTreeToBMM :: BMModel -> RTree a Double -> RTree a Double
scaleTreeToBMM bmm = fmap (* nSq)
  where nSq = fromIntegral (bmmPopSize bmm) ** 2

-- | Report the boundary mutation model specifications.
getBMMInfoStr :: BMModel
              -> Maybe Double
              -> Maybe [Double]
              -> String
getBMMInfoStr bmm ma mrs = unlines $
  [ "Population size: " ++ show (bmmPopSize bmm)
  , "Heterozygosity: " ++ show (bmmHeterozygosity bmm)
  , "Mutation model:"
  , getDNAModelInfoStr (bmmMutModel bmm)
  , "Gamma rate heterogeneity: " ++ show (isJust ma) ] ++
  gammaShape ++ gammaMeans
  where
    gammaShape = maybe [] (\a -> ["Shape parameter: " ++ show a]) ma
    gammaMeans = maybe [] (\rs -> ["This corresponds to uniformly distributed rates: " ++ show rs]) mrs
