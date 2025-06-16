
# G4 Discovery Pipeline

Scripts for annotating/predicting G-Quadruplexes (G4s) in a genome sequence, combining `pqsfinder` and `G4Hunter`.  

## Table of Contents


- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Features](#features)
- [Usage](#usage)
- [Notes](#notes)
- [References](#references)
- [Citation](#citation)

  

## Overview

This repository provides a Python script for predicting G-quadruplex (G4) structures in any FASTA sequence. The tool processes the FASTA file format and outputs a BED file containing non-overlapping G4s for each strand.

### Workflow:

1.  The script first runs the `pqsfinder` tool on the user-provided sequence to identify all potential overlapping G4s. It filters the G4s based on a user-defined `pqsfinder` score threshold (e.g., 30), and generates an output file with the extension `.fa.pqs`.
2.  The script, takes the output of the previous step as input and then calculates the `G4Hunter` score for each identified G4, ensuring that the G4 motif is consistent with those defined in the `g4DiscoveryFuncs.py` script.
3.  G4s with fewer than a specified number of tetrads (e.g., 3) and scores below the set thresholds for both `pqsfinder` (e.g., 40) and `G4Hunter` (e.g., 1.5) are filtered out.
4. Finally, G4s are grouped by starting position, with the highest-scoring, shortest G4 selected from each region to ensure non-overlapping, stable G4s.

## Prerequisites

Before using this package, ensure the following prerequisites are met: 
1. **Docker Installed**: 
	- Install Docker if it has not already been installed on your system. 
	- Refer to the [Docker Installation Guide](https://docs.docker.com/get-docker/) for platform-specific instructions. 
2. **Required Docker Container**: 
	- The package requires the container `kxk302/pqsfinder:1.0.0`. 
	- To pull the container, after installation, run the following command: `docker pull kxk302/pqsfinder:1.0.0 ` 
	- To verify the installation of the required container, run: `docker images`

*For more information on the dockerized version of pqsfinder, please refer to the repository at: [kxk302/pqsfinder-docker](https://github.com/kxk302/PqsFinder_Docker)*

## Features

-   **Dockerized Execution**: Fully containerized to run independently without requiring R language/packages.
-   **Flexible Motif Detection**: Supports both standard `((G{3,}[ATCG]{1,12}){3,}G{3,})` and bulged `((G([ATC]{0,1})G([ATC]{0,1})G([ATCG]{1,3})){3,}G([ATC]{0,1})G([ATC]{0,1})G)` G4 motifs.
-   **Non-overlapping G4 Detection**: Identifies non-overlapping G4 motifs on a given strand and prioritizes the most stable G4s within a region. 

## Usage

### Command-line Usage
**Running G4 Discovery**:

Use case: `g4Discovery.py [-h] -fa FASTA_FILE -chr CHROMOSOME -o OUTPUT [-t TETRAD] [-ps PQSSCORE] [-hs G4HUNTER] [-psd DOCKER_MIN_PQSSCORE]`

```
options:
-h, --help  show this help message and exit
-fa FASTA_FILE, --fasta_file FASTA_FILE
			Path to the input FASTA file
-chr CHROMOSOME, --chromosome CHROMOSOME
			Chromosome identifier, either an integer or a single-letter
-o OUTPUT, --output OUTPUT
			Path to the output BED file
-t TETRAD, --tetrad TETRAD
			Minimum number of tetrads for a G4 to be considered
-ps PQSSCORE, --pqsscore PQSSCORE
			Minimum pqsfinder score for a G4 to be considered
-hs G4HUNTER, --g4hunter G4HUNTER
			Minimum absolute G4Hunter score for a G4 to be considered
-psd DOCKER_MIN_PQSSCORE, --docker_min_pqsscore DOCKER_MIN_PQSSCORE
			Minimum pqsfinder score for the docker to run
```

### Example Use Case
`python3 g4Discovery.py -fa ../test/test.fa -chr 1 -o ../output/out.bed`

## Notes 

  - The input FASTA file should contain only one sequence (e.g. sequence from one chromosome), with a single identifier that starts with the `>` symbol (e.g. `>chr1 human CHM13`).
  - The docker daemon must be active in the background for the python script to run successfully.

## References
1. Hon, J., Martínek, T., Zendulka, J., & Lexa, M. (2017). [pqsfinder: an exhaustive and imperfection-tolerant search tool for potential quadruplex-forming sequences in R](https://doi.org/10.1093/bioinformatics/btx413). _Bioinformatics_, _33_(21), 3373-3379. `doi: 10.1093/bioinformatics/btx413`
2. Bedrat, A., Lacroix, L., & Mergny, J. L. (2016). [Re-evaluation of G-quadruplex propensity with G4Hunter](https://doi.org/10.1093/nar/gkw006). _Nucleic acids research_, _44_(4), 1746-1759. `doi: 10.1093/nar/gkw006`

## Citation
If you use this tool in your research, please cite the following paper:

> Mohanty, S. K., Chiaromonte F., & Makova, K. D. (2025). [Evolutionary dynamics of predicted G-quadruplexes in human and other great apes](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03635-1). *Genome Biology*, Vol. 26(161), `DOI: 10.1186/s13059-025-03635-1`
