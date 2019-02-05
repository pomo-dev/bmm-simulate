# A boundary mutation model simulator

For a given species tree (topology and branch lengths; Newick format), simulate
evolutionary sequences under the boundary mutation model which is used by
polymorphism-aware phylogenetic models (PoMo). The output is a [counts
file](http://www.iqtree.org/doc/Polymorphism-Aware-Models#counts-files).

The [PoMo reference
implementation](http://www.iqtree.org/doc/Polymorphism-Aware-Models) is part of
the [IQ-TREE software package](http://www.iqtree.org/).

For a detailed description of this software, please refer to

- Schrempf, D., Minh, B. Q., von Haeseler, A., & Kosiol, C., Polymorphism-aware
  species trees with advanced mutation models, bootstrap and rate heterogeneity,
  bioRxiv, (2019). http://dx.doi.org/10.1101/483479

For a description of the boundary mutation model and PoMo, please refer to the
following papers in reverse chronological order.

- Schrempf, D., & Hobolth, A., An alternative derivation of the stationary
  distribution of the multivariate neutral wright–fisher model for low mutation
  rates with a view to mutation rate estimation from site frequency data,
  Theoretical Population Biology, 114(), 88–94 (2017).
  http://dx.doi.org/10.1016/j.tpb.2016.12.001

- Schrempf, D., Minh, B. Q., De Maio, Nicola, von Haeseler, A., & Kosiol, C.,
  Reversible polymorphism-aware phylogenetic models and their application to
  tree inference, Journal of Theoretical Biology, 407(), 362–370 (2016).
  http://dx.doi.org/10.1016/j.jtbi.2016.07.042

- De Maio, N., Schrempf, D., & Kosiol, C., Pomo: an allele frequency-based
  approach for species tree estimation, Systematic Biology, 64(6), 1018–1031
  (2015). http://dx.doi.org/10.1093/sysbio/syv048

# Installation
`bmm-simulate` is written in [Haskell](https://www.haskell.org/) and can be
installed with [Stack](https://docs.haskellstack.org/en/stable/README/).

1. Install Stack with your package manager, or directly from the web page.
```
curl -sSL https://get.haskellstack.org/ | sh
```

2. Clone the `bmm-simulate` repository.
```
git clone https://github.com/pomo-dev/bmm-simulate.git
```
    
3. Navigate to the newly created `bmm-simulate` folder and build the simulator.
   This will take a while.
```
stack build
```

4. Run `bmm-simulate` from within the project directory.
```
stack exec bmm-simulate -- --help
```

5. If needed, install the simulator.
```
stack install
```

Now, the binary can be directly used.

# Help text
A copy of the output of `bmm-simulate --help` is pasted below.

```
Boundary mutation model simulator

Usage: bmm-simulate [-o|--output NAME] [-m|--mutation-model MODEL]
                    [--gamma-shape DOUBLE] [--gamma-ncat INT] [-N DOUBLE]
                    [-t|--heterozygosity DOUBLE] [-h|--tree-height DOUBLE]
                    [--tree-type TYPE] [--tree-yule-rate DOUBLE]
                    [-n|--nsites INT] [-s|--seed INT]
  Simulate count files using the boundary mutation model.

Available options:
  -h,--help                Show this help text
  -o,--output NAME         Write output files to
                           NAME.[cf.gz|log|tree] (default: "Test")
  -m,--mutation-model MODEL
                           Set the mutation model; available models are shown
                           below (default: HKY[6.25][0.3,0.2,0.2,0.3])
  --gamma-shape DOUBLE     Activate gamma rate heterogeneity and set gamma shape
                           parameter (default: off)
  --gamma-ncat INT         Set the number of gamma rate categories (no default
                           value)
  -N DOUBLE                Set the virtual population size (default: 9)
  -t,--heterozygosity DOUBLE
                           Set heterozygosity (default: 2.5e-3)
  -h,--tree-height DOUBLE  Set tree height [average number of
                           substitutions] (default: 5.0e-3)
  --tree-type TYPE         Set tree type; ILS or Yule (default: "ILS")
  --tree-yule-rate DOUBLE  Set the speciation rate of Yule tree (no default
                           value)
  -n,--nsites INT          Set number of sites to simulate (default: 1000000)
  -s,--seed INT            Set seed for the random number
                           generator (default: "random")

Available mutation models:
  - HKY model with transition to transversion ratio kappa and a state frequency vector.
    Specified with "-m HKY[DOUBLE][DOUBLE,DOUBLE,DOUBLE,DOUBLE]".
  - GTR model with five rate parameters and state frequency vector.
    Specified with "-m HKY[DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE][DOUBLE,DOUBLE,DOUBLE,DOUBLE]".

Note: The state frequency vector has to sum up to 1.0 and only has three free parameters.
```
