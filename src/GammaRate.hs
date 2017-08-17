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

import Rand
import Statistics

-- divideDist :: ContDist -> Distribution
-- getGammaRates :: Int -> 
