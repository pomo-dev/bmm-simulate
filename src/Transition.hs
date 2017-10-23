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

import           Control.Monad.Random.Strict hiding (fromList, toList)
import           Data.Distribution           hiding (toList)
import           Numeric.LinearAlgebra       hiding (fromList)
import           RateMatrix
import           RTree

-- This may be moved to a different module ProbMatrix or alike.
type ProbMatrix   = Matrix R

-- The important matrix that gives the probabilities to move from one state to
-- another in a specific time (branch length).
probMatrix :: RateMatrix -> BranchLn -> ProbMatrix
probMatrix m t = exponential
  where !exponential = expm $ scale t m

-- Convert a tree with branch lengths into a tree that has transition
-- probability matrices assigned to each of its branches.
branchLengthsToTransitionProbs :: RateMatrix -> RTree a Double -> RTree a ProbMatrix
branchLengthsToTransitionProbs m = fmap (probMatrix m)

-- Move from a given state to a new one according to a transition probability
-- matrix (for performance reasons this probability matrix needs to be given as
-- a list of generators, see
-- https://hackage.haskell.org/package/distribution-1.1.0.0/docs/Data-Distribution-Sample.html).
-- This function is the bottleneck of the simulator and takes up most of the
-- computation time. However, I was not able to find a faster implementation
-- that the one from Data.Distribution.
jump :: (MonadRandom m) => State -> [Generator State] -> m State
jump (State s) p = target
  where !target = getSample $ p !! s

-- Perform N jumps from a given state and according to a transition probability
-- matrix transformed to a list of generators. This implementation uses `foldM`
-- and I am not sure how to access or store the actual chain. This could be done
-- by an equivalent of `scanl` for general monads, which I was unable to find.
-- This function is neat, but will most likely not be needed. However, it is
-- instructive and is left in place.
jumpN :: (MonadRandom m) => State -> [Generator State] -> Int -> m State
jumpN s p n = foldM jump s (replicate n p)

-- This is the heart of the simulation. Take a tree (most probably with node
-- labels of type a, but this could be anything) and a root state. If there is a
-- split (i.e., if we have a node and not a leaf), jump down the left branch and
-- populate and flatten the tree and jump down the right branch and populate and
-- flatten the tree and append the two results. This is what `liftM2` is doing,
-- it lifts the append function of lists (++) to the `Rand g` monad. If we
-- encounter a leaf, just return the node label (or whatever the type a is) and
-- the state that we ended up at.
populateAndFlattenTree :: (MonadRandom m) => RTree a [Generator State] -> State -> m [(a, State)]
populateAndFlattenTree (Leaf a) s = return [(a, s)]
populateAndFlattenTree (Node _ lp lc rp rc) s = liftM2 (++) (jumpDownBranch lp lc) (jumpDownBranch rp rc)
  where jumpDownBranch p t = jump s p >>= populateAndFlattenTree t

stationaryDistToGenerator :: StationaryDist -> Generator State
stationaryDistToGenerator f = fG
  -- TODO: This is a little complicated. I need to convert the vector to a list
  -- to be able to create a distribution.
  where !fL = toList f
        !fD = fromList $ zip (map State [0..]) fL
        !fG = fromDistribution fD

treeProbMatrixToTreeGenerator :: RTree a ProbMatrix -> RTree a [Generator State]
treeProbMatrixToTreeGenerator t = tG
  where
    -- Create a tree with the probability matrices as list of row vectors.
    !tL = fmap toLists t
    -- A complicated double map. We need to create generators for each branch on
    -- the tree (fmap) and for each target state on each branch (map).
    !tG = (fmap . map) (fromDistribution . fromList . zip (map State [0..])) tL

-- Simulate data (states at the leaves) for a tree with transition probabilities
-- on its branches and with the stationary distribution of states at the root.
-- This function has to be impure because the state at the root is randomly
-- chosen from the stationary distribution and the states at the nodes and
-- leaves are randomly chosen according to the transition probabilities.
simulateSite :: (MonadRandom m) =>
                Generator State
             -> RTree a [Generator State]
             -> m [(a, State)]
simulateSite f t = do
  !rootState <- getSample f
  populateAndFlattenTree t rootState

-- Randomly draw an index according to a given generator. Use the stationary
-- distribution and rooted tree at the drawn index to simulate a site. This is
-- useful for simulation, e.g., Gamma rate heterogeneity models.
simulateSiteGen :: (MonadRandom m) =>
                   Generator Int
                -> [Generator State]
                -> [RTree a [Generator State]]
                -> m [(a, State)]
simulateSiteGen gen fs trs = do
  !i <- getSample gen
  let !f = fs  !! i
      !t = trs !! i
  simulateSite f t
