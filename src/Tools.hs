{- |
Description :  Some handy (mathematical, string manipulation, and so on) tools
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Please see function definitions and documentations.

* Changelog

-}

module Tools where

import Numeric.LinearAlgebra

-- | Get all values of a bounded enumerated type.
allValues :: (Bounded a, Enum a) => [a]
allValues = [minBound..]

-- | Calculate the nth harmonic number.
harmonic :: Int -> Double
harmonic 1 = 1.0
harmonic n = 1.0 / fromIntegral n + harmonic (n-1)

-- | Separate a square matrix into a symmetric and a skew-symmetric matrix.
matrixSeparateSymSkew :: Matrix R -> (Matrix R, Matrix R)
matrixSeparateSymSkew m = (mSym, mSkew)
  where trM = tr m
        mSym  = scale 0.5 $ m + trM
        mSkew = scale 0.5 $ m - trM

-- | Set the diagonal entries of a matrix to zero.
matrixSetDiagToZero :: Matrix R -> Matrix R
matrixSetDiagToZero m = m - diag (takeDiag m)

-- | Test for equality with tolerance (needed because of machine precision).
nearlyEq :: Double -> Double -> Double -> Bool
nearlyEq tol a b = tol > abs (a-b)
