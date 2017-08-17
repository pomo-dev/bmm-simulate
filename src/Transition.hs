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

Modules to look at:
- Data.Distribution
- Control.Monad.Trans.State.Lazy

* Changelog

-}

module Transition where

import           Control.Monad.Random.Strict
import qualified Numeric.LinearAlgebra       as L
import           Rand
import           RateMatrix
import           RTree

-- This may be moved to a different module ProbMatrix or alike.
type ProbMatrix   = L.Matrix L.R

-- The important matrix that gives the probabilities to move from one state to
-- another in a specific time (branch length).
probMatrix :: RateMatrix -> BranchLn -> ProbMatrix
probMatrix m t = L.expm $ L.scale t m

-- Convert a tree with branch lengths into a tree that has transition
-- probability matrices assigned to each of its branches.
branchLengthsToTransitionProbs :: RateMatrix -> RTree a Double -> RTree a ProbMatrix
branchLengthsToTransitionProbs m = fmap (probMatrix m)

-- Move from a given state to a new one according to a transition probability matrix.
jump :: (RandomGen g) => State -> ProbMatrix -> Rand g State
jump s p = drawFromDist (L.flatten $ p L.? [s])

-- Perform N jumps from a given state and according to a transition probability
-- matrix. This implementation uses `foldM` and I am not sure how to access or
-- store the actual chain. This could be done by an equivalent of `scanl` for
-- general monads, which I was unable to find. This function is neat, but will
-- most likely not be needed. However, it is instructive and is left in place.
jumpN :: (RandomGen g) => State -> ProbMatrix -> Int -> Rand g State
jumpN s p n = foldM jump s (replicate n p)

-- This is the heart of the simulation. Take a tree (most probably with node
-- labels of type a, but this could be anything) and a root state. If there is a
-- split (i.e., if we have a node and not a leaf), jump down the left branch and
-- populate and flatten the tree and jump down the right branch and populate and
-- flatten the tree and append the two results. This is what `liftM2` is doing,
-- it lifts the append function of lists (++) to the `Rand g` monad. If we
-- encounter a leaf, just return the node label (or whatever the type a is) and
-- the state that we ended up at.
populateAndFlattenTree :: (RandomGen g) => RTree a ProbMatrix -> State -> Rand g [(a, State)]
populateAndFlattenTree (Leaf a) s = return [(a, s)]
populateAndFlattenTree (Node _ lp lc rp rc) s = liftM2 (++) (jumpDownBranch lp lc) (jumpDownBranch rp rc)
  where jumpDownBranch p t = jump s p >>= populateAndFlattenTree t

-- Simulate data (states at the leaves) for a tree with transition probabilities
-- on its branches and with the stationary distribution of states at the root.
-- This function has to be impure because the state at the root is randomly
-- chosen from the stationary distribution and the states at the nodes and
-- leaves are randomly chosen according to the transition probabilities.
simulateSite :: (RandomGen g) => StationaryDist -> RTree a ProbMatrix -> Rand g [(a, State)]
simulateSite f t = do
  rootState <- drawFromDist f
  populateAndFlattenTree t rootState
-- Short, but less verbose:
-- simulateTree f t = drawFromDist f >>= populateAndFlattenTree t

-- Simulate n sites.
simulateNSites :: (RandomGen g) => Int -> StationaryDist -> RTree a ProbMatrix -> Rand g [[(a, State)]]
simulateNSites n f t = replicateM n $ simulateSite f t

-- Randomly draw an index according to a given distribution. Use the stationary
-- distribution and rooted tree at the drawn index to simulate a site. This is
-- useful for simulation, e.g., Gamma rate heterogeneity models.
simulateSiteDistr :: (RandomGen g) => Distribution -> [StationaryDist] -> [RTree a ProbMatrix] -> Rand g [(a, State)]
simulateSiteDistr dist fs trs = do
  i <- drawFromDist dist
  let f = fs  !! i
      t = trs !! i
  simulateSite f t
