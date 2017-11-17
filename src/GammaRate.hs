{- |
Description :  Convenience functions to simulate Gamma rate heterogeneity
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Sites evolve at different speeds. One way to model this phenomenon is to use
Gamma rate heterogeneity. This module provides the necessary definitions.

* Changelog

-}

module GammaRate where

import           Numeric.Integration.TanhSinh
import           Statistics.Distribution
import           Statistics.Distribution.Gamma

-- | For a given number of rate categories 'n' and a shape parameter 'alpha'
-- (the rate or scale is set such that the mean is 1.0), return a list of rates
-- that represent the respective categories. Use the mean rate for each
-- category.
getMeans :: Int -> Double -> [Double]
getMeans n alpha = means ++ lastMean
  where gamma = gammaDistr alpha (1.0/alpha)
        quantiles = [ quantile gamma (fromIntegral i / fromIntegral n) | i <- [0..n] ]
        -- Calculate the mean rate. Multiplication with the number of rate
        -- categories 'n' is necessary because in each n-quantile the
        -- probability mass is 1/n.
        meanFunc x = fromIntegral n * x * density gamma x
        -- Only calculate the first (n-1) categories with normal integration.
        means = [ intAToB meanFunc (quantiles !! i) (quantiles !! (i+1)) | i <- [0..n-2] ]
        -- The last category has to be calculated with an improper integration.
        lastMean = [intAToInf meanFunc (quantiles !! (n-1))]

-- | The error of integration.
eps :: Double
eps = 1e-6

-- | The integration method to use
method :: (Double -> Double ) -> Double -> Double -> [Result]
method = parSimpson

-- | Helper function for a normal integral from 'a' to 'b'.
intAToB :: (Double -> Double) -> Double -> Double -> Double
intAToB f a b = result . absolute eps $ method f a b

-- | Helper function for an improper integral from 'a' to infinity.
intAToInf :: (Double -> Double) -> Double -> Double
intAToInf f a = (result . absolute eps $ nonNegative method f) - intAToB f 0 a
