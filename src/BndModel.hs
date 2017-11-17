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
  ( MutModel
  , Heterozygosity
  , normalizeToTheta
  , normalizedRateMatrix
  , stationaryDist
  , scaleTreeToBMM
  , getBMMInfoStr )
  where

import           BndState
import           Data.Maybe
import           DNAModel              hiding (rateMatrix)
import           Numeric.LinearAlgebra
import           RateMatrix            hiding (State)
import           RTree
import           Tools

-- TODO: Let's just give rate matrices of the boundary mutation model a special
-- name. data BMModel = BMModel { bmmRateMatrix :: RateMatrix , bmmMutMatrix ::
-- MutModel , bmmPopSize :: PopSize }

-- | The boundary mutation models uses an underlying mutation model.
type MutModel = DNAModel

-- | First we need to define a mutation model.
mutRate :: MutModel -> Allele -> Allele -> Double
mutRate (DNAModel m _) a b = m ! i ! j
  where i = fromEnum a
        j = fromEnum b

-- | Get state frequency for specific allele.
stateFreq :: StateFreqVec -> Allele -> Double
stateFreq f a = f ! fromEnum a

-- | Calculate the rate of frequency shifts.
moranCoef :: PopSize -> Int -> Double
moranCoef n i = iD * (nD - iD) / nD
  where iD = fromIntegral i
        nD = fromIntegral n

-- | The transition rate from one boundary state to another.
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

-- | The transition rate from one state (index) to another.
rateById :: MutModel -> PopSize -> Int -> Int -> Double
rateById m n i j = rate m s t
  where s = idToState n i
        t = idToState n j

-- | The build function (see below) has a weird way of assigning entries to
-- indices. The indices have to be the same data type as the entries. This is
-- just a helper function that changes the indices from Double to Int.
rateByDouble :: MutModel -> PopSize -> Double -> Double -> Double
rateByDouble m n x y = rateById m n (round x) (round y)

-- | The rate matrix of the boundary mutation model. The dimension with four
-- alleles will be 4 + 6*(N-1).
rateMatrix :: MutModel -> PopSize -> RateMatrix
rateMatrix m n = setDiagonal $ build (s,s) (rateByDouble m n)
  where s = stateSpaceSize n

-- | Normalize the rate matrix such that on average one event happens per unit
-- time.
normalizedRateMatrix :: MutModel -> PopSize -> RateMatrix
normalizedRateMatrix m n = normalizeRates f' m'
  where m' = rateMatrix m n
        -- f  = dnaModelSpecGetStateFreqVec (dnaModelSpec m)
        f' = stationaryDist m n

-- | Define a heterozygosity to make function definitions clearer.
type Heterozygosity = Double

-- | Calculate the heterozygosity at stationarity.
theta :: MutModel -> Heterozygosity
theta (DNAModel m s) = f <.> rDiagZero #> f
  where f         = dnaModelSpecGetStateFreqVec s
        e         = toExchMatrix m f
        (r, _)    = matrixSeparateSymSkew e
        -- | The summation excludes the diagonal (a /= b).
        rDiagZero = matrixSetDiagToZero r

-- | The normalization constant of the stationary distribution.
norm :: MutModel -> PopSize -> Double
norm m n = 1.0 + harmonic (n-1) * theta m

-- | Get entries of the stationary measure (not normalized) for a boundary model state.
stationaryMeasEntry :: MutModel -> State -> Double
stationaryMeasEntry m s =
  let f   = dnaModelSpecGetStateFreqVec (dnaModelSpec m) in
  case s of
    (Bnd _ a)     -> f ! fromEnum a
    (Ply n i a b) -> piA * mAB / (fromIntegral n - fromIntegral i)
                                      + piB * mBA / fromIntegral i
      where piA = stateFreq f a
            mAB = mutRate m a b
            piB = stateFreq f b
            mBA = mutRate m b a

-- -- Only use this function when accessing single elements of the stationary
-- -- distribution. Otherwise, computation of the norm is unnecessarily repeated.
-- stationaryDistEntry :: MutModel -> PopSize -> State -> Double
-- stationaryDistEntry m n s = stationaryMeasEntry m s / norm m n

stationaryMeasEntryById :: MutModel -> PopSize -> Int -> Double
stationaryMeasEntryById m n i = stationaryMeasEntry m s
  where s = idToState n i

stationaryMeasEntryByDouble :: MutModel -> PopSize -> Double -> Double
stationaryMeasEntryByDouble m n x = stationaryMeasEntryById m n (round x)

-- | Get the stationary distribution fir a boundary mutation model.
stationaryDist :: MutModel -> PopSize -> StationaryDist
stationaryDist m n = scale (1.0/norm m n) $
  build s (stationaryMeasEntryByDouble m n)
  where s = stateSpaceSize n

-- | Normalize the mutation coefficients such that the heterozygosity matches a
-- given level (see Eq. 12.14 in my thesis).
normalizeToTheta :: MutModel -> PopSize -> Heterozygosity -> MutModel
normalizeToTheta mo@(DNAModel m _) n h =
  mo { dnaRateMatrix = scale (h / (t * (1.0 - c * h))) m }
  where
    -- | The heterozygosity of the boundary mutation model.
    t = theta mo
    -- | The branch length multiplicative factor introduced by the coalescent.
    c = harmonic (n-1)

-- | The branch lengths of threes in the boundary mutation model are not
-- measured in average number of substitutions per site but in average number of
-- mutations or frequency shifts per site. The conversion factor is just the
-- square of the population size. This function converts the branch lengths of a
-- tree.
scaleTreeToBMM :: PopSize -> RTree a Double -> RTree a Double
scaleTreeToBMM n = fmap (* nSq)
  where nSq = fromIntegral n ** 2

-- | Report the boundary mutation model specifications.
getBMMInfoStr :: Int
              -> Double
              -> MutModel
              -> Maybe Double
              -> Maybe [Double]
              -> String
getBMMInfoStr n h m ma mrs = unlines $
  [ "Population size: " ++ show n
  , "Heterozygosity: " ++ show h
  , "Mutation model:"
  , getDNAModelInfoStr m
  , "Gamma rate heterogeneity: " ++ show (isJust ma) ] ++
  gammaShape ++ gammaMeans
  where
    gammaShape = maybe [] (\a -> ["Shape parameter: " ++ show a]) ma
    gammaMeans = maybe [] (\rs -> ["This corresponds to uniformly distributed rates: " ++ show rs]) mrs
