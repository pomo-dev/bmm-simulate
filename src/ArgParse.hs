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

import Options.Applicative
import Data.Semigroup ((<>))

newtype BMSimArgs = BMSimArgs
  { outFileName :: String }

parseBMSimArgs :: IO BMSimArgs
parseBMSimArgs = execParser $ info (bmSimOptions <**> helper)
  ( fullDesc
    <> progDesc "Simulate count files using the boundary mutation model."
    <> header "Boundary mutation model simulator")

outFileNameP :: Parser String
outFileNameP = strOption
  ( long "output"
    <> metavar "OUTPUT"
    <> help "Output file name")

bmSimOptions :: Parser BMSimArgs
bmSimOptions = BMSimArgs
  <$> outFileNameP
