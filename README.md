# Progenepredict readme file

**Progenepredict software basis:**
- Python programming language (v3.9.2)
- package for array computing NumPy (v1.20.1)

## Introduction and workflow
### **Introduction:**

Progenepredict is a tool for *de novo* gene prediction in prokaryotic genomes. Its prediction process is based on investigating all reading frames for the presence of open reading frames (ORFs) with Shine-Dalgarno (SD) sequences in the vicinity of their start codon.

### **Workflow:**

***Input file:*** A FASTA file containing the sequence corresponding to the forward strand (5' to 3') of a prokaryotic genome.

***Procedure:***

> **1.** The reverse complement strand is created based on the leading strand
>
>  **2.** Iteration in steps of three over the leading strand employing the following mechanism to scan for ORFs. For each reading-frame:
>
>>  **2.1** The current codon is created
>>
>>  **2.2** The codon is assessed to see if it is:
>>
>>>	* A start codon preceded by a SD sequence (assessed by employing the Smith-Waterman sequence alignment algorithm)
>>>
>>>	* A stop codon
>>>
>>>	* None of the above

Step two is continuously executed until a start codon is found that identifies the start of an ORF associated with a coding sequence. Subsequently:

>	**3.** The iteration over the strand and codon assessment continues until a stop codon is found (any additional start codons do not reinitialize the current ORF)
>
>  **4.** When a stop codon is found the length of the complete ORF is assessed and based on the user-specified criterion of minimum length the predicted gene is either:
>>  *criterion fulfilled:* Written to the output file with associated information
>>
>> *criterion not fulfilled:* Discarded
>
> **5.** The process continues from **(1)** until the end of the input sequence is reached
>
> **6.** The same process takes place for the reverse complement strand


***Output file(s):***
1. A **tab-separated output file** that contains the start and end position of each predicted gene, the start codon that initialized the ORF, as well as the reading frame and strand it was found on

2. A **FASTA file** containing the nucleotide sequences of the predicted genes *(optional)*

## Installation and setup
Recommended installation of required software with conda:

**1. Installing conda**
```bash
# download the installer corresponding to your os from https://docs.conda.io/en/latest/miniconda.html
# in this example miniconda3 is installed for linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# run script
Miniconda3-latest-Linux-x86_64.sh
# update and configure conda
conda update conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda -forge
```
**2. Installing NumPy**
```bash
# create new environment with python version 3.9.2
conda create -n progenepredict python=3.9.2
# activate it
conda activate progenepredict
# install NumPy version 1.20.1
conda install NumPy=1.20.1
# check that the correct software version were installed
conda list | grep -E "numpy|python"
```

Alternative ways to install the required software are found at:

Python: https://www.python.org

NumPy: https://numpy.org

## Usage
Now, Progenepredict can be run the conda environment progenepredict that was created in the previous step. The output files will be created in the directory the program is run from. 

**1. Usage instructions**
```bash
progenepredict.py [-h] -g GENOMEFILE [-o OUTPUTFILE] [-ml MINORFLENGTH] [-sc STARTCODONS] [-ec STOPCODONS] [-sd SHINEDALGARNOSEQUENCE] [-me MAXALIGNERRORS] [-f]
```

* **-h**	show help message with usage instructions
* **-g**	the file containing the genome sequence in FASTA format
* **-o**	the desired name of the output file containing the predicted genes
* **-ml**	the minimum ORF length of predicted genes *(default=300)*
* **-sc**	the possible start codons *(default=ATG)*
* **-ec**	the possible stop codons *(default=TAA TAG TGA)*
* **-sd**	the Shine-Dalgarno sequence *(default=AGGAGG)*
* **-me**	the max. number of errors that may be present in the pairwise alignment when checking for a Shine-Dalgarno sequence *(default=1)*
* **-f**	output a file with the nucleotide sequences of predicted genes in FASTA format

**Note about the default settings:**
The default values for arguments were set based on the features of the E. coli genome and are therefore recommended to be adapted to the species the respective input genome originates from.
