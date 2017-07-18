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

import Numeric.LinearAlgebra
import qualified BndModel as BM
import qualified DNAModel as DNA

main :: IO ()
main = do
  let f = vector [0.3, 0.2, 0.2, 0.3]
      n = 9
      het = 0.0031
      q = BM.normalizeToTheta (DNA.rateMatrixHKY f 6.2157) f n het
      m = BM.rateMatrix q n
  print m
  return ()
