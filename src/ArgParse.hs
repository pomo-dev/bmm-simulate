{- |
   Module      :  ArgParse
   Description :  Parse arguments needed for the bmm simulator.
   Copyright   :  (c) Dominik Schrempf 2017
   License     :  GPLv3

   Maintainer  :  dominik.schrempf@gmail.com
   Stability   :  unstable
   Portability :  non-portable (not tested)

Provides a data type with all options and a parser.

* Changelog

-}

module ArgParse where

import qualified Data.Attoparsec.Text  as A
import           Data.Semigroup        ((<>))
import qualified Data.Text             as T
import qualified Defaults              as Def
import qualified DNAModel              as DNA
import qualified Numeric.LinearAlgebra as L
import           Options.Applicative

-- Convenience function to read in more complicated command line options with
-- attoparsec and optparse
-- (https://github.com/pcapriotti/optparse-applicative#option-readers).
attoReadM :: A.Parser a -> ReadM a
attoReadM p = eitherReader (A.parseOnly p . T.pack)

data BMSimArgs = BMSimArgs
  { outFileName    :: String
  , stateFreqs     :: L.Vector L.R
  , kappa          :: Double
  , popSize        :: Int
  , heterozygosity :: Double
  , treeHeight     :: Double
  , nSites         :: Int
  , seed           :: String }

-- The impure IO action that reads the arguments and prints out help if needed.
-- Maybe put this into Main.hs?
parseBMSimArgs :: IO BMSimArgs
parseBMSimArgs = execParser $ info (helper <*> bmSimOptions)
  ( fullDesc
    <> progDesc "Simulate count files using the boundary mutation model."
    <> header "Boundary mutation model simulator")


-- General things and options.
outFileNameOpt :: Parser String
outFileNameOpt = strOption
  ( long "output"
    <> short 'o'
    <> metavar "FILENAME"
    <> value Def.outFileName
    <> showDefault
    <> help "Write output to FILENAME (counts file format)")

-- TODO: Define default values at the top (probably move all this to Main.hs).
-- Option to input the stationary frequencies of the mutation model.
stateFreqsOpt :: Parser (L.Vector L.R)
stateFreqsOpt = option (attoReadM parseStateFreq)
  ( long "freq"
    <> short 'f'
    <> metavar "pi_A,pi_C,pi_G,pi_T"
    <> value Def.stateFreqs
    <> showDefault
    <> help "Set the stationary frequencies of the nucleotides")

-- Read a stationary frequency of the form `pi_A,pi_C,pi_G,...`.
parseStateFreq :: A.Parser DNA.StateFreqVec
parseStateFreq = do
  f <- L.vector . take nAlleles <$> A.sepBy A.double (A.char ',')
  if L.norm_1 f == 1.0 then return f
    else error $ "Stationary frequencies sum to " ++ show (L.norm_1 f) ++ " but should sum to 1.0."
    -- TODO: Nucleotide count hard coded. See `BndState.nAlleles`. However, I cannot
    -- take different numbers depending on the mutation model because optparse has
    -- an applicative interface and not a monadic one. Also, allowing vectors of
    -- arbitrary size is undesirable because of meaningless error messages.
    where nAlleles = 4

kappaOpt :: Parser Double
kappaOpt = option auto
  ( long "kappa"
    <> short 'k'
    <> metavar "VALUE"
    <> value Def.kappa
    <> showDefault
    <> help "Set the kappa value for the HKY model" )

popSizeOpt :: Parser Int
popSizeOpt = option auto
  ( short 'N'
    <> metavar "VALUE"
    <> value Def.popSize
    <> showDefault
    <> help "Set the virtual population size" )

heterozygosityOpt :: Parser Double
heterozygosityOpt = option auto
  ( long "heterozygosity"
    <> short 't'
    <> metavar "VALUE"
    <> value Def.heterozygosity
    <> showDefault
    <> help "Set the heterozygosity" )

treeHeightOpt :: Parser Double
treeHeightOpt = option auto
  ( long "treeheight"
    <> short 'h'
    <> metavar "VALUE"
    <> value Def.treeHeight
    <> showDefault
    <> help "Set the tree height [average number of substitutions]" )

nSitesOpt :: Parser Int
nSitesOpt = option auto
  ( long "nsites"
    <> short 'n'
    <> metavar "VALUE"
    <> value Def.nSites
    <> showDefault
    <> help "Set the number of sites to simulate" )

seedOpt :: Parser String
seedOpt = strOption
  ( long "seed"
    <> short 's'
    <> metavar "VALUE"
    <> value "random"
    <> showDefault
    <> help "Set the seed for the random number generator" )

-- Composition of all options.
bmSimOptions :: Parser BMSimArgs
bmSimOptions = BMSimArgs
  <$> outFileNameOpt
  <*> stateFreqsOpt
  <*> kappaOpt
  <*> popSizeOpt
  <*> heterozygosityOpt
  <*> treeHeightOpt
  <*> nSitesOpt
  <*> seedOpt
