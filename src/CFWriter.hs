{- |
Module      :  CFWriter
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

module CFWriter where

import           BndState  as BS
import           System.IO

-- The number of populations (leafs) on the tree.
type NPop = Int

-- The number of sites that will be printed.
type NSites = Int

-- The names of the populations.
type PopulationNames = [String]

-- | Compose the header using the number of sites and the population names.
header :: NSites -> PopulationNames -> String
header nSites popNames = lineOne ++ "\n" ++ lineTwo ++ "\n"
  where nPop = length popNames
        lineOne = "COUNTSFILE NPOP " ++ show nPop ++ " NSITES " ++ show nSites
        lineTwo = "CHROM POS " ++ unwords popNames

-- The chromosome name.
type Chrom = String

-- The position on the chromosome.
type Pos   = Int

-- The set of boundary states for one site.
type DataOneSite = [BS.State]

-- Get a data line in the counts file.
getDataLine :: Chrom -> Pos -> DataOneSite -> String
getDataLine chrom pos bstates = chrom ++ " " ++ show pos ++ " " ++ bStatesString ++ "\n"
  where bStatesString = unwords $ map show bstates

-- At the moment, all data is passed around in one go, if I have problems, I
-- have to separate writing on a line per line basis.
type DataAllSites = [DataOneSite]

-- Write a counts file. For now, the chromosome name is set to SIM for simulated
-- and the position is just a counter starting at 1.
write :: FilePath -> NSites -> PopulationNames -> DataAllSites -> IO ()
write fileName nSites popNames bStatesArray = do
  handle <- openFile fileName WriteMode
  hPutStr handle $ header nSites popNames
  let dataLines = zipWith (getDataLine "SIM") [1..] bStatesArray
  mapM_ (hPutStr handle) dataLines
  hClose handle
