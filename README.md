# Code repo for "The evolution of thermal performance curves in fungi farmed by attine ant mutualists in aboveground or belowground microclimates"
This repository contains code for the publication "The evolution of thermal performance curves in fungi farmed by attine ant mutualists in aboveground or belowground microclimates" by Oscar Hess, Tiem van der Deure, Mille Bolander, Caio Leal Dutra, and Jonathan Z. Shik. 

In this paper, we analyze the performance of fungus symbionts cultivated by six fungus-farming ants across a temperature gradient.

## Reproducibility
### Phylogenetic trees
The shell script was run using MAFFT v7.490 (available as a Tar file [here](https://mafft.cbrc.jp/alignment/software/mafft-7.490-without-extensions-src.tgz)) and IQ-TREE v 1.6.12 (available from GitHub [here](https://github.com/Cibiv/IQ-TREE/releases/tag/v1.6.12))

After installing the packages, run the script `phylotree.sh` in the shell to generate the tree. The sequences and alignment are available as .fasta files in this repository.

### Data analysis
All data analysis was done in Julia v.1.10. We recommend using [juliaup](https://github.com/JuliaLang/juliaup) to install this version of Julia.

After installing Julia and cloning this repository, start Julia in the directory that contains this code. To replicate the environment that were used to run the analysis, type `]` to enter pkg mode, then `activate .` to activate the local environment, and then  `instantiate` to download the packages specified in `Manifest.toml`.

The julia script `01_sampling.jl` generates and saves the posterior distributions. The figures and tables in the paper are generated in the script `02_sampling.jl`.

Note that the sampling script makes heavy use of multi-threading to speed up computations. It is highly recommended to start julia with mutiple threads. See the official [manual](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to do this.
