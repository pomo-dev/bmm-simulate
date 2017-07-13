{- |
Module      :  RTree
Description :  A binary, rooted tree structure with branch lengths
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

A binary, rooted tree with branch lengths is an important structure in
phylogenetic analysis.

* Changelog

-}

module RTree
  ( RTree(..)
  ) where

-- The tree data type with states of type a. The branch length of type b is the
-- length from the current to the left and right child.
data RTree a b = Node { state :: a
                      , lBrLn :: b
                      , lChld :: RTree a b
                      , rBrLn :: b
                      , rChld :: RTree a b }
               | Leaf { state :: a }
               deriving (Eq, Show, Read)

totalBrLn :: Num b => RTree a b -> b
totalBrLn (Leaf _) = 0
totalBrLn t = lBrLn t  + rBrLn t
                      + totalBrLn (lChld t)
                      + totalBrLn (rChld t)

-- In order to simulate data, I need the probabilities of finding specific
-- states at the tips of a given tree when starting at a given distribution of
-- states at the root. This distribution will be the stationary distribution in
-- my case.

-- preOrderTraversal :: Num b => (a -> a -> b -> Double) -> RTree a b -> Double
-- -- The likelihood of a leaf is 1.0.
-- preOrderTraversal _ (Leaf _) = 1.0
-- preOrderTraversal transitionProb (Node s lb lc rb rc)
--   = (transitionProb s (getState lc) lb) * (transitionProb s (getState rc) rb)
--     * preOrderTraversal transitionProb lc
--     * preOrderTraversal transitionProb rc

getLeaves :: RTree a b -> [RTree a b]
getLeaves (Node _ _ lc _ rc) = getLeaves lc ++ getLeaves rc
getLeaves leaf = [leaf]

-- Some examples and tests.
myLeftLeaf :: RTree Char b
myLeftLeaf = Leaf 'l'

myRightLeaf :: RTree Char b
myRightLeaf = Leaf 'r'


myTree :: RTree Char Integer
myTree = Node 'r' 1 myLeftLeaf 2 myRightLeaf

