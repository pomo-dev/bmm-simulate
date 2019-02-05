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

1. Install stack with your package manager, or directly from the web page.

    curl -sSL https://get.haskellstack.org/ | sh

2. Clone the `bmm-simulate` repository.

    git clone https://github.com/pomo-dev/bmm-simulate.git
    
3. Navigate to the newly created `bmm-simulate` folder and build the simulator.
   This takes a while.

    stack build

4. Run `bmm-simualate` from within the project directory.

    stack exec bmm-simulate -- --help

5. If needed, install the simulator.

    stack install
