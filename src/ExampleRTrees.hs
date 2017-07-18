{- |
Module      :  Module Name
Description :  A collection of example trees
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

A collection of example trees to be used for simulating sequence data under the
boundary mutation model.

* Changelog

-}

module ExampleRTrees where

import RTree

-- The ILS tree.
-- Newick representation for 1Ne: (((s4:0.5,s3:0.5):0.1,s2:0.6):0.4,s1:1.0);.
ilsTree :: Double -> RTree String Double
ilsTree th = Node "root"
             (th/2-0.1) (Node "intern1"
                           0.1 (Node "intern2"
                                (th/2) (Leaf "s4")
                                (th/2) (Leaf "s3"))
                           (th/2+0.1) (Leaf "s2"))
             th (Leaf "s1")

ilsTree1Ne :: RTree String Double
ilsTree1Ne = ilsTree 1.0
