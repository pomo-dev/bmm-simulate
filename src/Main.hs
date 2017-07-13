{- |
Module      :  BMM simulator
Description :  Simulate populations using the boundary mutation model
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

Enables simulation of sequence data for a predefined phylogeny using the
boundary mutation model. The output is a counts file.

* Changelog

-}


module Main where

import Numeric.LinearAlgebra
import BndModel
import DNAModel
import RateMatrix

main :: IO ()
main = do
  let f = vector [0.3, 0.2, 0.2, 0.3]
      m = rateMatrixHKY f 6.0
      q = rateMatrixBM m 9
      p = probMatrix q 0.3
  print q
  print p
  return ()
