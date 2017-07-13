{- |
Module      :  Boundary mutation model states
Description :  An implementation of the state space of the boundary mutation model.
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

The boundary mutation model is a discrete-state, continuous-time Markov process
that allows mutations only when the population is monomorphic.

* Changelog

-}

module BndModState
  ( Nuc(..)
  , Allele
  , PopSize
  , AlleleCount
  , Bnd(..)
  , Ply(..)
  , BState(..)
  , stateSpace
  , stateSpaceSize
  , bStateInd
  , indToBState
  , connected
  ) where

import Data.List

-- First, we need to define the state space.
-- The nucleotides.
data Nuc = A | C | G | T deriving (Eq, Show, Read, Ord, Bounded, Enum)
-- Alleles are just nucleotides at the moment.
type Allele = Nuc
-- The population size; has to be larger than one, otherwise there be dragons.
type PopSize = Int
-- The absolute frequency of an allele.
type AlleleCount = Int

-- The number of alleles.
nAlleles :: Int
nAlleles = 1 + fromEnum (maxBound :: Allele)

-- A boundary mutation model state is either fixed (Bnd for boundary) or
-- polymorphic (Ply).
data Bnd = Bnd { bndN :: PopSize
               , bndA :: Allele }
           deriving (Eq, Show, Read)

validBnd :: Bnd -> Bool
validBnd (Bnd n _)
  | n <= 1    = False
  | otherwise = True

filterValidBnd :: [Bnd] -> [Bnd]
filterValidBnd = filter validBnd

-- A total order on the boundary states. The sorting order is N > A.
instance Ord Bnd where
  b1 <= b2
    | b1 == b2           = True
    | bndN b1 /= bndN b2 = bndN b1 <= bndN b2
    | otherwise          = bndA b1 <= bndA b2

data Ply = Ply { plyN :: PopSize
               , plyI :: AlleleCount
               , plyA :: Allele
               , plyB :: Allele }
           deriving (Eq, Show, Read)

validPly :: Ply -> Bool
validPly (Ply n i a b)
  | n <= 1    = False
  | a >= b    = False
  | i <= 0    = False
  | i >= n    = False
  | otherwise = True

filterValidPly :: [Ply] -> [Ply]
filterValidPly = filter validPly

-- A total order on the polymorphic states. The sorting order is N > A > B > I.
instance Ord Ply where
  p1 <= p2
    | p1 == p2            = True
    | plyN p1 /= plyN p2  = plyN p1 <= plyN p2
    | plyA p1 < plyA p2   = True
    | plyA p1 > plyA p2   = False
    -- We can be sure that plyA p1 == plyA p2 now.
    | plyB p1 < plyB p2   = True
    | plyB p1 > plyB p2   = False
    -- Now we can be sure that both nucleotides are the same.
    | otherwise           = plyI p1 < plyI p2

-- A boundary mutation model state is either a boundary state or a polymorphic
-- state. This also automatically defines a total order.
data BState = BStateBnd Bnd | BStatePly Ply
             deriving (Eq, Show, Read, Ord)

-- validBState :: BState -> Bool
-- validBState (BStateBnd b) = validBnd b
-- validBState (BStatePly p) = validPly p

getPopSize :: BState -> PopSize
getPopSize (BStateBnd b)  = bndN b
getPopSize (BStatePly p)  = plyN p

-- | Sorted list of all possible PoMo states for a specific population size.
stateSpace :: PopSize -> [BState]
stateSpace n
  | n <= 1    = error "The population size has to be larger than one."
  | otherwise = sort ( bndStates ++ plyStates )
  where allBndStates = [ Bnd n a |
                         a <- [minBound .. maxBound] :: [Allele] ]
        bndStates    = map BStateBnd $ filterValidBnd allBndStates
        allPlyStates = [ Ply n i a b |
                         i <- [0..n],
                         a <- [minBound .. maxBound] :: [Allele],
                         b <- [minBound .. maxBound] :: [Allele] ]
        plyStates    = map BStatePly $ filterValidPly allPlyStates

stateSpaceSize :: PopSize -> Int
stateSpaceSize n = k + k*(k-1) `div` 2 * (n-1)
  where k = nAlleles

bStateInd :: BState -> Maybe Int
bStateInd s = elemIndex s (stateSpace $ getPopSize s)

indToBState :: PopSize -> Int -> BState
indToBState n i = stateSpace n !! i

-- Check if two states are connected. By definition, states are NOT connected
-- with themselves.
connected :: BState -> BState -> Bool
connected s t = s `elem` getNeighbors t

getNeighbors :: BState -> [BState]
getNeighbors (BStateBnd b) = map BStatePly $ getBndNeighbors b
getNeighbors (BStatePly p) = getPlyNeighbors p

-- Neighbors of boundary states are always polymorphic states.
getBndNeighbors :: Bnd -> [Ply]
getBndNeighbors (Bnd n a) = filterValidPly allNeighbors
  where allNeighbors = [ Ply n (n-1) a b |
                         b <- [minBound .. maxBound] :: [Allele] ]
                       ++
                       [ Ply n 1 b a |
                         b <- [minBound .. maxBound] :: [Allele] ]

-- Neighbors of polymorphic states can be polymorphic or boundary states.
getPlyNeighbors :: Ply -> [BState]
getPlyNeighbors (Ply n i a b)
  -- Careful when the population size is two, because then each polymorphic
  -- states has two boundary states as neighbors.
  | i == 1 && n == 2  = BStateBnd (Bnd n a) : [BStateBnd (Bnd n b)]
  | i == 1            = BStateBnd (Bnd n b) : [BStatePly (Ply n 2 a b)]
  | i == (n-1)        = BStateBnd (Bnd n a) : [BStatePly (Ply n (n-2) a b)]
  | otherwise         = BStatePly (Ply n (i+1) a b) : [BStatePly (Ply n (i-1) a b)]

-- -- Examples and tests.
-- bndOne = Bnd 9  A
-- bndTwo = Bnd 9  C
-- bndThr = Bnd 10 T

-- plyOne = Ply 9  2 A C
-- plyTwo = Ply 9  3 A C
-- plyThr = Ply 9  7 C T
-- plyFou = Ply 10 8 A C

-- stateOne = BStatePly plyOne
-- stateTwo = BStatePly plyTwo
-- stateThr = BStateBnd bndOne
