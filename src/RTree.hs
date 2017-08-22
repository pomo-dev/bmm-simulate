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
  ( BranchLn
  , RTree(..)
  , totalBrLn
  , getLeaves
  , toNewick
  , ilsTree
  ) where

-- Branch lengths on trees are measured in Double.
type BranchLn = Double

-- The strict tree data type with node names or states of type a. The branch
-- length of type b is the length from the current to the left and right child.
data RTree a b = Node { state :: !a
                      , lBrLn :: !b
                      , lChld :: RTree a b
                      , rBrLn :: !b
                      , rChld :: RTree a b }
               | Leaf { state :: !a }
               deriving (Eq, Show, Read)

-- Make (RTree a) a functor. Like this, we can scale branch lengths or convert
-- them to transition probability matrices.
instance Functor (RTree a) where
  fmap _ (Leaf a) = Leaf a
  fmap f (Node a lb lc rb rc) = Node a (f lb) (fmap f lc) (f rb) (fmap f rc)

-- The total branch length; only works when the branch lengths are numbers.
totalBrLn :: RTree a BranchLn -> BranchLn
totalBrLn (Leaf _) = 0
totalBrLn t = lBrLn t + rBrLn t
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

-- -- Some examples and tests.
-- myLeftLeaf :: RTree Char b
-- myLeftLeaf = Leaf 'l'

-- myRightLeaf :: RTree Char b
-- myRightLeaf = Leaf 'r'


-- myTree :: RTree Char Integer
-- myTree = Node 'r' 1 myLeftLeaf 2 myRightLeaf

toNewick :: RTree String BranchLn -> String
toNewick t = toNewick' t ++ ";"
  where
    toNewick' (Leaf a) = a
    toNewick' (Node _ lb lc rb rc) = "(" ++
                                     toNewick' lc ++ ":" ++ show lb ++ "," ++
                                     toNewick' rc ++ ":" ++ show rb ++ ")"

-- The ILS tree with tree height `th`.
-- Newick representation for height 1.0: (((s4:0.5,s3:0.5):0.1,s2:0.6):0.4,s1:1.0);.
ilsTree :: Double -> RTree String Double
ilsTree th = Node "root"
             (th/2.0 - th/10.0) (Node "intern1"
                                 (th/10.0) (Node "intern2"
                                            (th/2) (Leaf "s4")
                                            (th/2) (Leaf "s3"))
                                 (th/2.0 + th/10.0) (Leaf "s2"))
             th (Leaf "s1")
