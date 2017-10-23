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
  , Scenario
  , totalBrLn
  , getLeaves
  , toNewick
  , ils
  , yule
  ) where

import Data.Random
import Data.Random.Distribution.Exponential

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

-- | The simulation scenario. It is a tree type with parameters.
--
-- At the moment, ILS (incomplete lineage sorting) and Yule trees are supported.
data Scenario = ILS  { ilsTreeHeight  :: Double }
              | Yule { yuleTreeHeight :: Double
                     , yuleRate       :: Double }

instance Show Scenario where
  show s = case s of
    ILS h -> unlines [ showType "ILS"
                     , showHeight h ]
    Yule h r -> unlines [ showType "Yule"
                        , showHeight h
                        , showRate r ]
    -- _ -> error "Tree type not supported."
    where showType tt = "Tree type: " ++ tt
          showHeight th = "Tree height in average number of substitutions: " ++ show th
          showRate rt = "Yule speciation rate: " ++ show rt

-- The total branch length; only works when the branch lengths are numbers.
totalBrLn :: RTree a BranchLn -> BranchLn
totalBrLn (Leaf _) = 0
totalBrLn t = lBrLn t + rBrLn t
                      + totalBrLn (lChld t)
                      + totalBrLn (rChld t)

getLeaves :: RTree a b -> [a]
getLeaves (Node _ _ lc _ rc) = getLeaves lc ++ getLeaves rc
getLeaves leaf = [state leaf]

toNewick :: RTree String BranchLn -> String
toNewick t = toNewick' t ++ ";"
  where
    toNewick' (Leaf a) = a
    toNewick' (Node _ lb lc rb rc) = "(" ++
                                     toNewick' lc ++ ":" ++ show lb ++ "," ++
                                     toNewick' rc ++ ":" ++ show rb ++ ")"

-- | The ILS tree with tree height 'th'.
-- Newick representation for height 1.0: (((s4:0.5,s3:0.5):0.1,s2:0.6):0.4,s1:1.0);.
ils :: Double -> (RTree String BranchLn, Scenario)
ils th = (tree, ILS th)
  where tree = Node "root"
               (th/2.0 - th/10.0) (Node "intern1"
                                   (th/10.0) (Node "intern2"
                                              (th/2) (Leaf "s4")
                                              (th/2) (Leaf "s3"))
                                   (th/2.0 + th/10.0) (Leaf "s2"))
               th (Leaf "s1")

-- | Yule tree with height 'th' and speciation rate 'rate'.
yule :: Double -> Double -> RVar (RTree String BranchLn, Scenario)
yule th rate = do
  let scenario = Yule th rate
  tree <- yule' th rate
  return (tree, scenario)
  where yule' h r = do
          lBrLnSample <- exponential (1.0/r)
          rBrLnSample <- exponential (1.0/r)
          let lBrLen = if lBrLnSample >= h then h else lBrLnSample
              rBrLen = if rBrLnSample >= h then h else rBrLnSample
          lChild <- if lBrLnSample >= h then pure (Leaf "") else
                      yule' (h - lBrLen) r
          rChild <- if rBrLnSample >= h then pure (Leaf "") else
                      yule' (h - rBrLen) r
          return $ labelLeaves (Node "" lBrLen lChild rBrLen rChild)

-- | Set all node states to the empty string "" and label the leaves from 0 to
-- the number of leaves.
labelLeaves :: RTree String BranchLn -> RTree String BranchLn
labelLeaves t = fst $ labelLeaves' t 0

labelLeaves' :: RTree a b -> Int -> (RTree String b, Int)
labelLeaves' (Leaf _) n = (Leaf (show n), n+1)
labelLeaves' (Node _ lB lC rB rC) n = (Node "" lB lC' rB rC', n'')
  where (lC', n' ) = labelLeaves' lC n
        (rC', n'') = labelLeaves' rC n'
