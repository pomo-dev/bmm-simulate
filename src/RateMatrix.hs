{- |
Module      :  Rate Matrix
Description :  Rate matrix helper functions
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable 
Portability :  non-portable (not tested)

Some helper functions that come handy when working with rate matrices of
continuous-time discrete-state Markov processes.

* Changelog

-}

module RateMatrix where

import Numeric.LinearAlgebra

-- A rate matrix is just a real matrix.
type RateMatrix     = Matrix R
-- A matrix of exchangeabilities, we have q = e * pi, where q is a rate matrix,
-- e is the exchangeability matrix and pi is the diagonal matrix containing the
-- stationary frequency distribution.
type ExchMatrix     = Matrix R
-- Stationary distribution of a rate matrix.
type StationaryDist = Vector R

matrixSetDiagToZero :: Matrix R -> Matrix R
matrixSetDiagToZero m = m - diag (takeDiag m)

-- Normalizes a Markov process generator such that one event happens per unit time.
rateMatrixNormalize :: StationaryDist -> RateMatrix -> RateMatrix
rateMatrixNormalize f m = scale (1.0 / totalRate) m
  where totalRate = norm_1 $ f <# matrixSetDiagToZero m

-- Set the diagonal entries of a matrix such that the rows sum to 0.
rateMatrixSetDiagonal :: RateMatrix -> RateMatrix
rateMatrixSetDiagonal m = diagZeroes - diag (fromList rowSums)
  where diagZeroes = matrixSetDiagToZero m
        rowSums    = map norm_1 $ toRows diagZeroes

-- Extract the exchangeability matrix from a rate matrix.
rateMatrixToExchMatrix :: RateMatrix -> StationaryDist -> ExchMatrix
rateMatrixToExchMatrix m f = m <> diag oneOverF
  where oneOverF = cmap (1.0/) f

-- This may be moved to a different module ProbMatrix or alike.
type BranchLength = Double
type ProbMatrix   = Matrix R

-- The important matrix that gives the probabilities to move from one state to
-- another in a specific time (branch length).
probMatrix :: RateMatrix -> BranchLength -> ProbMatrix
probMatrix m t = expm $ scale t m

