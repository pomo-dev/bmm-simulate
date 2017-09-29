{-# LANGUAGE BangPatterns #-}
{- |
   Module      :  Transition
   Description :  Functions to work with transition probability matrices on rooted trees
   Copyright   :  (c) Dominik Schrempf 2017
   License     :  GPLv3

   Maintainer  :  dominik.schrempf@gmail.com
   Stability   :  unstable
   Portability :  non-portable (not tested)

Calculate transition probability matrices, map rate matrices on trees, populate
a tree with states according to a stationary distribution, etc.

The implementation of the Markov process is more than basic and can be improved in a lot of ways.

* Changelog

-}

module Transition where

import           Control.Monad.Random.Strict
import qualified Data.Distribution           as D
import qualified Numeric.LinearAlgebra       as L
import qualified RateMatrix                  as RM
import qualified RTree                       as Tree

-- This may be moved to a different module ProbMatrix or alike.
type ProbMatrix   = L.Matrix L.R

-- The important matrix that gives the probabilities to move from one state to
-- another in a specific time (branch length).
probMatrix :: RM.RateMatrix -> Tree.BranchLn -> ProbMatrix
probMatrix m t = exponential
  where !exponential = L.expm $ L.scale t m

-- Convert a tree with branch lengths into a tree that has transition
-- probability matrices assigned to each of its branches.
branchLengthsToTransitionProbs :: RM.RateMatrix -> Tree.RTree a Double -> Tree.RTree a ProbMatrix
branchLengthsToTransitionProbs m = fmap (probMatrix m)

-- Move from a given state to a new one according to a transition probability
-- matrix (for performance reasons this probability matrix needs to be given as
-- a list of generators, see
-- https://hackage.haskell.org/package/distribution-1.1.0.0/docs/Data-Distribution-Sample.html).
jump :: (RandomGen g) => RM.State -> [D.Generator RM.State] -> Rand g RM.State
jump (RM.State s) p = target
  where !target = D.getSample $ p !! s

-- Perform N jumps from a given state and according to a transition probability
-- matrix transformed to a list of generators. This implementation uses `foldM`
-- and I am not sure how to access or store the actual chain. This could be done
-- by an equivalent of `scanl` for general monads, which I was unable to find.
-- This function is neat, but will most likely not be needed. However, it is
-- instructive and is left in place.
jumpN :: (RandomGen g) => RM.State -> [D.Generator RM.State] -> Int -> Rand g RM.State
jumpN s p n = foldM jump s (replicate n p)

-- This is the heart of the simulation. Take a tree (most probably with node
-- labels of type a, but this could be anything) and a root state. If there is a
-- split (i.e., if we have a node and not a leaf), jump down the left branch and
-- populate and flatten the tree and jump down the right branch and populate and
-- flatten the tree and append the two results. This is what `liftM2` is doing,
-- it lifts the append function of lists (++) to the `Rand g` monad. If we
-- encounter a leaf, just return the node label (or whatever the type a is) and
-- the state that we ended up at.
populateAndFlattenTree :: (RandomGen g) => Tree.RTree a [D.Generator RM.State] -> RM.State -> Rand g [(a, RM.State)]
populateAndFlattenTree (Tree.Leaf a) s = return [(a, s)]
populateAndFlattenTree (Tree.Node _ lp lc rp rc) s = liftM2 (++) (jumpDownBranch lp lc) (jumpDownBranch rp rc)
  where jumpDownBranch p t = jump s p >>= populateAndFlattenTree t

stationaryDistToGenerator :: RM.StationaryDist -> D.Generator RM.State
stationaryDistToGenerator f = fG
  -- TODO: This is a little complicated. I need to convert the vector to a list
  -- to be able to create a distribution.
  where !fL = L.toList f
        !fD = D.fromList $ zip (map RM.State [0..]) fL
        !fG = D.fromDistribution fD

treeProbMatrixToTreeGenerator :: Tree.RTree a ProbMatrix -> Tree.RTree a [D.Generator RM.State]
treeProbMatrixToTreeGenerator t = tG
  where
    -- Create a tree with the probability matrices as list of row vectors.
    !tL = fmap L.toLists t
    -- A complicated double map. We need to create generators for each branch on
    -- the tree (fmap) and for each target state on each branch (map).
    !tG = (fmap . map) (D.fromDistribution . D.fromList . zip (map RM.State [0..])) tL

-- Simulate data (states at the leaves) for a tree with transition probabilities
-- on its branches and with the stationary distribution of states at the root.
-- This function has to be impure because the state at the root is randomly
-- chosen from the stationary distribution and the states at the nodes and
-- leaves are randomly chosen according to the transition probabilities.
simulateSite :: (RandomGen g) =>
                D.Generator RM.State
             -> Tree.RTree a [D.Generator RM.State]
             -> Rand g [(a, RM.State)]
simulateSite f t = do
  !rootState <- D.getSample f
  populateAndFlattenTree t rootState

-- Randomly draw an index according to a given generator. Use the stationary
-- distribution and rooted tree at the drawn index to simulate a site. This is
-- useful for simulation, e.g., Gamma rate heterogeneity models.
simulateSiteGen :: (RandomGen g) =>
                   D.Generator Int
                -> [D.Generator RM.State]
                -> [Tree.RTree a [D.Generator RM.State]]
                -> Rand g [(a, RM.State)]
simulateSiteGen gen fs trs = do
  !i <- D.getSample gen
  let !f = fs  !! i
      !t = trs !! i
  simulateSite f t
