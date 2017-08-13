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
import qualified DNAModel              as DNA
import qualified Numeric.LinearAlgebra as LinAlg
import           Options.Applicative

-- Convenience function to read in more complicated command line options with
-- attoparsec and optparse
-- (https://github.com/pcapriotti/optparse-applicative#option-readers).
attoReadM :: A.Parser a -> ReadM a
attoReadM p = eitherReader (A.parseOnly p . T.pack)

data BMSimArgs = BMSimArgs
  { outFileName :: String
  , stateFreqs  :: LinAlg.Vector LinAlg.R}

-- General things and options.
parseBMSimArgs :: IO BMSimArgs
parseBMSimArgs = execParser $ info (helper <*> bmSimOptions)
  ( fullDesc
    <> progDesc "Simulate count files using the boundary mutation model."
    <> header "Boundary mutation model simulator")

outFileNameOpt :: Parser String
outFileNameOpt = strOption
  ( long "output"
    <> short 'o'
    <> metavar "FILE"
    <> help "Write output to FILE (counts file format)")

-- TODO: Define default values at the top (probably move all this to Main.hs).
-- Option to input the stationary frequencies of the mutation model.
stateFreqsOpt :: Parser (LinAlg.Vector LinAlg.R)
stateFreqsOpt = option (attoReadM parseStateFreq)
  ( long "freq"
    <> short 'f'
    <> metavar "pi_A,pi_C,pi_G,pi_T]"
    <> value (LinAlg.vector [0.3, 0.2, 0.2, 0.3])
    <> help "Set the stationary frequencies of the nucleotides")

-- Read a stationary frequency of the form `pi_A,pi_C,pi_G,...`.
parseStateFreq :: A.Parser DNA.StateFreqVec
-- TODO: Nucleotide count hard coded. See `BndState.nAlleles`. However, I cannot
-- take different numbers depending on the mutation model because optparse has
-- an applicative interface and not a monadic one. Also, allowing vectors of
-- arbitrary size is undesirable because of meaningless error messages.
parseStateFreq = do
  f <- LinAlg.vector <$> take 4 <$> (A.sepBy A.double (A.char ','))
  if LinAlg.norm_1 f == 1.0 then return f
    else error $ "Stationary frequencies sum to " ++ show (LinAlg.norm_1 f) ++ " but should sum to 1.0."

-- Composition of all options.
bmSimOptions :: Parser BMSimArgs
bmSimOptions = BMSimArgs
  <$> outFileNameOpt
  <*> stateFreqsOpt
