{- |
Module      :  Distribution
Description :  Random values, distributions, etc.
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Impure.

Some functions to work with random computations such as a definition of a
distribution, and convenience functions to draw from distributions.

* Changelog


-}

module Distribution where

-- I had problems with using =scanl= for vectors, a function provided by
-- Data.Vector. However, hmatrix reexports some definitions of vector which
-- results in errors. The qualified imports solve this problems.
import           Control.Monad.Random.Strict
import qualified Data.Vector.Generic         as V
import qualified Numeric.LinearAlgebra       as L

-- A discrete distribution is a real vector.
type Distribution = L.Vector L.R

-- Randomly sample an index from a discrete probability distribution.
drawFromDist :: (RandomGen g) => Distribution -> Rand g Int
drawFromDist dist = do
  p <- getRandomR (0.0 :: Double, 1.0 :: Double)
  return $ fromDist p dist

-- A probability is just a double, an index is just an Int.
type Probability = Double

-- For a given value, return the index of the sample according to a given
-- distribution. In Control.Monad.Trans.Random, this is implemented for lists,
-- but I happen to use vectors in this program, so I had to rewrite this
-- function.
fromDist :: Probability -> Distribution -> Int
fromDist p dist = V.length $ V.takeWhile (<p) cums
  -- The cumulative distribution.
  where cums = V.scanl1' (+) dist

-- Get a uniform distribution with n categories.
uniformDist :: Int -> Distribution
uniformDist n = L.vector $ replicate n $ 1.0 / fromIntegral n
