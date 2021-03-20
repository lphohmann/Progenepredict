# README file for application example 1 of Progenepredict

## Introduction
In this example Progenepredict was utilised to predict genes in the genome of Escherichia coli str. K-12 substr. MG1655.

## Files and usage
***Input file:***

**1.** Escherichia_coli_str._K-12_substr._MG1655.fasta

***Output files:***

**1.** Escherichia_coli_str._K-12_substr._MG1655_predictions

**2.** Escherichia_coli_str._K-12_substr._MG1655_predictions.fasta

***Usage:***
```bash
progenepredict.py -g Escherichia_coli_str._K-12_substr._MG1655.fasta -o Escherichia_coli_str._K-12_substr._MG1655_predictions -ml 300 -sc ATG GTG TTG -ec TAA TAG TGA -sd AGGAGG -me 1 -f
```

## Data

The input genome of this example was obtained from the National Center for Biotechnology Information (NCBI Reference Sequence: NC_000913.3) at https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3?report=fasta.
