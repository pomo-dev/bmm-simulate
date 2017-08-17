{- |
Module      :  GammaRate
Description :  Convenience functions to simulate Gamma rate heterogeneity.
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

-- import qualified Distribution                  as D
import           Statistics.Distribution
import           Statistics.Distribution.Gamma

-- | For a given number of rate categories 'n' and a shape parameter 'alpha'
-- (the rate or scale is set such that the mean is 1.0), return a list of rates
-- that represent the respective categories. Use the mean rate for each
-- category.
getGammaRatesMean :: Int -> Double -> [Double]
getGammaRatesMean n alpha = quantiles
  where gamma = gammaDistr alpha (1.0/alpha)
        quantiles = [ quantile gamma (fromIntegral i / fromIntegral n) | i <- [0..n] ]
