{- |
Description :  Write a counts file
Copyright   :  (c) Dominik Schrempf 2017
License     :  GPLv3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  non-portable (not tested)

After simulating data according to a tree and a mutation model (and a population
size, etc.), write the data to file using counts file format.

* Changelog

-}

module CFWriter
  ( Pos
  , DataOneSite
  , NSites
  , PopulationNames
  , Chrom
  , writeHeader
  , writeLine
  , open
  ) where

import           BndState
import           System.IO
import           Tools     (left, right)

-- | The number of sites that will be printed.
type NSites = Int

-- | The names of the populations.
type PopulationNames = [String]

-- The column width of the counts file.
colW :: Int
colW = 11

-- | Compose the header using the number of sites and the population names.
header :: NSites -> PopulationNames -> String
header nSites popNames = lineOne ++ "\n" ++ lineTwo ++ "\n"
  where nPop = length popNames
        lineOne = "COUNTSFILE NPOP " ++ show nPop ++ " NSITES " ++ show nSites
        lineTwo = left colW "CHROM" ++ right colW "POS" ++ unwords (map (right colW) popNames)

-- | The chromosome name.
type Chrom = String

-- | The position on the chromosome.
type Pos   = Int

-- | The set of boundary states for one site.
type DataOneSite = [State]

-- | Get a data line in the counts file.
getDataLine :: Chrom -> Pos -> DataOneSite -> String
getDataLine chrom pos bstates = left colW chrom ++ right colW (show pos) ++ bStatesString ++ "\n"
  where bStatesString = unwords $ map (right colW . show) bstates

-- | I am not sure why this is needed. Maybe to make imports redundant? Can be
-- removed.
open :: FilePath -> IO Handle
open file = openFile file WriteMode

-- | Write a counts file. For now, the chromosome name is set to SIM for
-- simulated and the position is just a counter starting at 1.
writeHeader :: Handle -> NSites -> PopulationNames -> IO ()
writeHeader handle nSites popNames =
  hPutStr handle $ header nSites popNames

-- | Write a counts file line including chromosome (which will be SIM for the
-- simulation), position (consecutively numbered) and the data.
writeLine :: Handle -> Chrom -> Pos -> DataOneSite -> IO ()
writeLine handle chr pos bStates =
  hPutStr handle $ getDataLine chr pos bStates
