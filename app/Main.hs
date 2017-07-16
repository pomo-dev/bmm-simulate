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

TODO DEBUG DEBUG

-}

import Numeric.LinearAlgebra
import BndModel
import DNAModel

main :: IO ()
main = do
  let f = vector (replicate 4 0.25)
      q = rateMatrixHKY f 6.0
      n = 9
      het = 0.0025
      qNorm = bmNormalizeToTheta q f n het
      m = rateMatrixBM qNorm n
  print m
  return ()
