# Code repo for "The evolution of thermal performance curves in fungi farmed by attine ant mutualists in aboveground or belowground microclimates"
This repository contains code for the publication "The evolution of thermal performance curves in fungi farmed by attine ant mutualists in aboveground or belowground microclimates" by Oscar Hess, Tiem van der Deure, Mille Bolander, Caio Leal Dutra, and Jonathan Z. Shik. 

In this paper, we analyze the performance of fungus symbionts cultivated by six fungus-farming ants across a temperature gradient.

## Corresponding author information
Name: Jonathan Shik
Institution: University of Copenhagen
Address: Universitetsparken 15, 2100 Copenhagen
Email: jonathan.shik@bio.ku.dk


## Usage notes: 
These datasets were generated in at the University of Copenhagen based on colonies attine ant colonies (Atta colombica, Acromyrmex echinatior, Paratrachymyrmex cornetzi, Cyphomyrmex costatus, Mycocepurus smithii, Apterostigma dentigerum) collected at the Smithsonian Tropical Research Institute in Panama (Soberenia Park).

## Description of data and code
This archive contains all data and code to reproduce all statistical analyses and produce raw versions of the figures, as described in the Methods and reported in the Results of the associated paper.

Datasets included:
1) data.csv
2) metadata.csv
3) alignment.fasta
4) sequences.fasta

### data.csv: 
This file contains results of the laboratory in vitro experiment with fungi grown at different temperatures
* Number of variables: 8
* Number of cases/rows: 358
* Variable List:
  * species: The attine ant species from which the fungus was isolated 
  * colony:  Colony IDs representing individual colonies (biological replicates) from which the fungal sample was isolated 
  * temperature: the temperature treatment (degrees C) at which the Petri dish was incubated (10, 15, 20, 25, 30)
  * area_d0: the area of the fungus (cm^2) measured on Day 
  * area_d7: the area of the fungus (cm^2) measured on Day 7
  * area_d14: the area of the fungus (cm^2) measured on Day 14
  * area_d21: the area of the fungus (cm^2) measured on Day 7
  * area_d28: the area of the fungus (cm^2) measured on Day 28

### metadata.csv: 
This file contains metadata for data.csv.
 * Number of variables: 4
 * Number of cases/rows: 19
 * Variable List:
  * species: The attine ant species from which the fungus was isolated 
  * colony:  Colony IDs representing individual colonies (biological replicates) from which the fungal sample was isolated 
  * Ecology: Whether the fungus is farmed aboveground or belowground (Categorical variable)

### alignment.fasta: 
Aligned ITS sequences for the barcoding analyses and the tree in Figure 1 of the publication
* Number of colonies seqenced: 20

### sequences.fasta: 
Raw ITS sequences for the barcoding analyses and the tree in Figure 1 of the publication
* Number of colonies seqenced: 20

## Reproducibility
### Phylogenetic trees
The shell script was run using MAFFT v7.490 (available as a Tar file [here](https://mafft.cbrc.jp/alignment/software/mafft-7.490-without-extensions-src.tgz)) and IQ-TREE v 1.6.12 (available from GitHub [here](https://github.com/Cibiv/IQ-TREE/releases/tag/v1.6.12))

After installing the packages, run the script `phylotree.sh` in the shell to generate the tree. The script uses `sequences.fasta` and `alignment.fasta`.

### Data analysis
All data analysis was done in Julia v.1.10. We recommend using [juliaup](https://github.com/JuliaLang/juliaup) to install this version of Julia.

After installing Julia and cloning this repository, start Julia in the directory that contains this code. To replicate the environment that were used to run the analysis, type `]` to enter pkg mode, then `activate .` to activate the local environment, and then  `instantiate` to download the packages specified in `Manifest.toml`.

The julia script `01_sampling.jl` generates and saves the posterior distributions. The figures and tables in the paper are generated in the script `02_sampling.jl`. These scripts load data from `data.csv` and `metadata.csv`, as well as helper functions and additional code in `data.jl`, `bayesian_model.jl`, and `response_functions`.

Note that the sampling script makes heavy use of multi-threading to speed up computations. It is highly recommended to start julia with mutiple threads. See the official [manual](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to do this.
