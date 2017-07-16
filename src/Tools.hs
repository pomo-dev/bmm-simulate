{- |
Module      :  Tools
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

harmonic :: Int -> Double
harmonic 1 = 1.0
harmonic n = 1.0 / (fromIntegral n) + harmonic (n-1)

-- Separate a matrix into a symmetric and a skew-symmetric matrix.
separateMatrixSymSkew :: Matrix R -> (Matrix R, Matrix R)
separateMatrixSymSkew m = (sym, skew)
  where trM = tr m
        sym  = scale 0.5 $ m + trM
        skew = scale 0.5 $ m - trM
