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

module ArgParse
  ( BMMArgs
  , parseBMMArgs
  , seed
  , dnaModelSpec
  , gammaNCat
  , gammaShape
  , popSize
  , heterozygosity
  , treeHeight
  , treeType
  , treeYuleRate
  , nSites
  , outFileName )
where

import qualified Data.Attoparsec.Text         as A
import           Data.Semigroup               ((<>))
import           Data.Text                    (pack, unpack)
import qualified Defaults                     as Def
import           DNAModel                     (DNAModelSpec (..), StateFreqVec)
import           Numeric.LinearAlgebra        (norm_1, size, vector)
import           Options.Applicative
import qualified Text.PrettyPrint.ANSI.Leijen as Doc
import           Tools                        (nearlyEq)

-- Convenience function to read in more complicated command line options with
-- attoparsec and optparse
-- (https://github.com/pcapriotti/optparse-applicative#option-readers).
attoReadM :: A.Parser a -> ReadM a
attoReadM p = eitherReader (A.parseOnly p . pack)

data BMMArgs = BMMArgs
  { outFileName    :: String
  -- , stateFreqs     :: Maybe (Vector R)
  , dnaModelSpec   :: DNAModelSpec
  , gammaShape     :: Maybe Double
  , gammaNCat      :: Maybe Int
  , popSize        :: Int
  , heterozygosity :: Double
  , treeHeight     :: Double
  , treeType       :: String
  , treeYuleRate   :: Maybe Double
  , nSites         :: Int
  , seed           :: String }

-- Composition of all options.
bmSimOptions :: Parser BMMArgs
bmSimOptions = BMMArgs
  <$> outFileNameOpt
  -- <*> stateFreqsOpt
  <*> dnaModelOpt
  <*> gammaShapeOpt
  <*> gammaNCatOpt
  <*> popSizeOpt
  <*> heterozygosityOpt
  <*> treeHeightOpt
  <*> treeTypeOpt
  <*> treeYuleRateOpt
  <*> nSitesOpt
  <*> seedOpt

-- The impure IO action that reads the arguments and prints out help if needed.
-- Maybe put this into Main.hs?
parseBMMArgs :: IO BMMArgs
parseBMMArgs = execParser $
  info (helper <*> bmSimOptions)
  (fullDesc
    <> progDesc "Simulate count files using the boundary mutation model."
    <> header "Boundary mutation model simulator"
    <> footerDoc models )
  where
    models = Just $ foldl1 (Doc.<$>) (map Doc.text strs)
    strs   = [ "Available mutation models:"
             , "  - HKY model with transition to transversion ratio kappa."
             , "    Specified with \"-m HKY[DOUBLE] -f [DOUBLE,DOUBLE,DOUBLE,DOUBLE]\"." ]

-- General things and options.
outFileNameOpt :: Parser String
outFileNameOpt = strOption
  ( long "output"
    <> short 'o'
    <> metavar "FILEPATH"
    <> value Def.outFileName
    <> showDefault
    <> help "Write output to FILEPATH in counts file format" )

-- -- Option to input the stationary frequencies of the mutation model.
-- stateFreqsOpt :: Parser (Maybe (Vector R))
-- stateFreqsOpt = option (attoReadM parseStateFreq)
--   ( long "freq"
--     <> short 'f'
--     <> metavar "[DOUBLE,DOUBLE,DOUBLE,DOUBLE]"
--     <> value Def.stateFreqs
--     <> showDefault
--     <> help msg )
--   where
--     msg = "Set the stationary frequencies of the nucleotides in order A, C, G and T; " ++
--           "not applicable if mutation model does not support different nucleotide frequencies"

-- Read a stationary frequency of the form `pi_A,pi_C,pi_G,...`.
parseStateFreq :: A.Parser StateFreqVec
parseStateFreq = do
  _ <- A.char '['
  f <- vector <$> A.sepBy A.double (A.char ',')
  _ <- A.char ']'
  if size f /= nAlleles
    then error "Length of stationary frequency vector is faulty, only DNA models are supported."
  else if nearlyEq 1e-6 (norm_1 f) 1.0 then return f
    else error $ "Stationary frequencies sum to " ++ show (norm_1 f) ++ " but should sum to 1.0."
    -- Nucleotide count hard coded. See `BndState.nAlleles`. However, I cannot
    -- take different numbers depending on the mutation model because optparse has
    -- an applicative interface and not a monadic one. Also, allowing vectors of
    -- arbitrary size is undesirable because of meaningless error messages.
    where nAlleles = 4

dnaModelOpt :: Parser DNAModelSpec
dnaModelOpt = option (attoReadM parseDNAModelSpec)
  ( long "mutation-model"
    <> short 'm'
    <> metavar "MODEL"
    <> value Def.dnaModelSpec
    <> showDefault
    <> help "Set the mutation model; see available models below" )

parseParams :: A.Parser [Double]
parseParams = do
  _ <- A.char '['
  params <- A.sepBy1 A.double (A.char ',')
  _ <- A.char ']'
  return params

parseDNAModelSpec :: A.Parser DNAModelSpec
parseDNAModelSpec = do
  m  <- A.takeWhile (/= '[')
  case unpack m of
       "HKY" -> do
         ps <- parseParams
         f  <- parseStateFreq
         if length ps /= 1
           then error "HKY model only has one parameter, kappa."
           else return $ HKY (head ps) f
       "JC" -> return JC
       _ -> error "Model string could not be parsed."

gammaShapeOpt :: Parser (Maybe Double)
gammaShapeOpt = optional $ option auto
  ( long "gamma-shape"
    <> metavar "DOUBLE"
    <> help "Set gamma shape parameter; activate gamma rate heterogeneity (default: off)" )

gammaNCatOpt :: Parser (Maybe Int)
gammaNCatOpt = optional $ option auto
  ( long "gamma-ncat"
    <> metavar "INT"
    <> help "Set the number of gamma rate categories (no default value)" )

popSizeOpt :: Parser Int
popSizeOpt = option auto
  ( short 'N'
    <> metavar "DOUBLE"
    <> value Def.popSize
    <> showDefault
    <> help "Set the virtual population size" )

heterozygosityOpt :: Parser Double
heterozygosityOpt = option auto
  ( long "heterozygosity"
    <> short 't'
    <> metavar "DOUBLE"
    <> value Def.heterozygosity
    <> showDefault
    <> help "Set heterozygosity" )

treeHeightOpt :: Parser Double
treeHeightOpt = option auto
  ( long "tree-height"
    <> short 'h'
    <> metavar "DOUBLE"
    <> value Def.treeHeight
    <> showDefault
    <> help "Set tree height [average number of substitutions]" )

treeTypeOpt :: Parser String
treeTypeOpt = strOption
  ( long "tree-type"
    <> metavar "TYPE"
    <> value Def.treeType
    <> showDefault
    <> help "Set tree type; ILS or Yule")

treeYuleRateOpt :: Parser (Maybe Double)
treeYuleRateOpt = optional $ option auto
  ( long "tree-yule-rate"
  <> metavar "DOUBLE"
  <> help "Set the speciation rate of Yule tree (no default value)")

nSitesOpt :: Parser Int
nSitesOpt = option auto
  ( long "nsites"
    <> short 'n'
    <> metavar "INT"
    <> value Def.nSites
    <> showDefault
    <> help "Set number of sites to simulate" )

seedOpt :: Parser String
seedOpt = strOption
  ( long "seed"
    <> short 's'
    <> metavar "INT"
    <> value "random"
    <> showDefault
    <> help "Set seed for the random number generator" )
