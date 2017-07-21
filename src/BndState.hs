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

module BndState
  ( Allele
  , PopSize
  , AlleleCount
  , BState(..)
  , stateSpace
  , stateSpaceSize
  , bStateId
  , idToBState
  , connected
  ) where

import           Data.List
import           DNAModel
import           RateMatrix
import           Tools

-- First, we need to define the state space.

-- Alleles are just nucleotides at the moment. However, I want to keep the code
-- such that it can be extended easily to codons or amino acids.
type Allele = Nuc
-- The population size; has to be larger than one, otherwise there be dragons.
type PopSize = Int
-- The absolute frequency of an allele.
type AlleleCount = Int

-- The number of alleles.
nAlleles :: Int
nAlleles = 1 + fromEnum (maxBound :: Allele)

-- A boundary mutation model state is either a boundary state or a polymorphic
-- state. This also automatically defines a total order.
data BState = Bnd { bndN :: PopSize
                  , bndA :: Allele }
            | Ply { plyN :: PopSize
                  , plyI :: AlleleCount
                  , plyA :: Allele
                  , plyB :: Allele }
             deriving (Eq, Read)

instance Show BState where
  show (Bnd n a) =  foldl1 (++) $ intersperse "," $ map toCounts allValues
    where toCounts b
            | a == b    = show n
            | otherwise = "0"
  show (Ply n i a b) = foldl1 (++) $ intersperse "," $ map toCounts allValues
    where toCounts c
            | c == a    = show i
            | c == b    = show (n-i)
            | otherwise = "0"

-- A total order on the boundary mutation model states. In general, Bnd < Ply.
-- Then, sorting happens according to the order population size, first allele,
-- second allele, allele count. It may be beneficial to reverse the allele count
-- order (i.e., make a polymorphic state with higher allele count show up before
-- a polymorphic state with lower allele count, this would move some polymorphic
-- states closer to their respective boundaries),
instance Ord BState where
  Bnd {} <= Ply {}            = True
  Ply {} <= Bnd {}            = False
  s@(Bnd n a) <= t@(Bnd m b)
    | s == t                  = True
    | n /= m                  = n <= m
    | otherwise               = a <= b
  s@(Ply n i a b) <= t@(Ply m j c d)
    | s == t                  = True
    | n /= m                  = n <= m
    | a < c                   = True
    | a > c                   = False
    -- We can be sure that a  = c now.
    | b < d                   = True
    | b > d                   = False
    -- Now we can be sure that both nucleotides are the same.
    | otherwise               = i <= j

validBState :: BState -> Bool
validBState (Bnd n _)
  | n <= 1    = False
  | otherwise = True
validBState (Ply n i a b)
  | n <= 1    = False
  | a >= b    = False
  | i <= 0    = False
  | i >= n    = False
  | otherwise = True

filterValidBState :: [BState] -> [BState]
filterValidBState = filter validBState

getPopSize :: BState -> PopSize
getPopSize (Bnd n _)     = n
getPopSize (Ply n _ _ _) = n

-- | Sorted list of all possible PoMo states for a specific population size.
stateSpace :: PopSize -> [BState]
stateSpace n
  | n <= 1    = error "The population size has to be larger than one."
  | otherwise = sort $ filterValidBState ( allBndStates ++ allPlyStates )
  where allBndStates = [ Bnd n a |
                         a <- [minBound .. maxBound] :: [Allele] ]
        allPlyStates = [ Ply n i a b |
                         i <- [0..n],
                         a <- [minBound .. maxBound] :: [Allele],
                         b <- [minBound .. maxBound] :: [Allele] ]

stateSpaceSize :: PopSize -> Int
stateSpaceSize n = k + k*(k-1) `div` 2 * (n-1)
  where k = nAlleles

bStateId :: BState -> Maybe Int
bStateId s = elemIndex s (stateSpace $ getPopSize s)

idToBState :: PopSize -> State -> BState
idToBState n i = stateSpace n !! i

-- Check if two states are connected. By definition, states are NOT connected
-- with themselves.
connected :: BState -> BState -> Bool
connected s t = s `elem` getNeighbors t

getNeighbors :: BState -> [BState]
getNeighbors (Bnd n a) = filterValidBState allNeighbors
  where allNeighbors = [ Ply n (n-1) a b |
                         b <- [minBound .. maxBound] :: [Allele] ]
                       ++
                       [ Ply n 1 b a |
                         b <- [minBound .. maxBound] :: [Allele] ]
getNeighbors (Ply n i a b)
  -- Careful when the population size is two, because then each polymorphic
  -- states has two boundary states as neighbors.
  | i == 1 && n == 2  = Bnd n a : [Bnd n b]
  | i == 1            = Bnd n b : [Ply n 2 a b]
  | i == (n-1)        = Bnd n a : [Ply n (n-2) a b]
  | otherwise         = Ply n (i+1) a b : [Ply n (i-1) a b]

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
