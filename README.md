# EchoSV
## Overview
EchoSV is a versatile tool for comparing and merging structural variant (SV) call sets that were generated using different reference genomes. It studies how SVs "echo" across these references through a hybrid workflow that combines lift-over and graph-based matching.
![EchoSV Workflow](echosv_workflow.pdf)
Given two or more SV call sets from the same sample—each aligned to a different reference—EchoSV can perform two primary operations:
- **Merge**: Consolidates multiple SV call sets into a single, unified output. For example, it can merge two DSA haplotype–based call sets into one consolidated file.
- **Compare**: Generates a detailed comparison identifying overlapping variants and those exclusive to a specific reference, such as when analyzing calls across GRCh38, CHM13, and DSA.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)  
- [Usage](#usage)  
- [License](#license)  
- [Contact](#contact)  

## Requirements

EchoSV depends on the following Python packages:
- **pysam**: Read and write BAM/CRAM files and VCF records for variant processing.  
- **intervaltree**: Efficiently store and query genomic intervals to detect overlapping SVs.  
- **Biopython (Bio)**: Parse and manipulate sequence data during liftover steps.  
- **scipy**: Perform statistical analysis and numerical computations on variant metrics.  
- **networkx**: Construct and traverse graphs that model structural variant matches.  

## Installation

Option 1: From GitHub (recommended)

```bash
# Clone the repo and install
git clone git@github.com:parklab/EchoSV.git
cd EchoSV
pip install .
```

Option 2: Via PyPI

```bash
pip install echosv
```

## Usage

EchoSV workflow consists of three main steps: **chain**, **genotype**, and **match**. Below are detailed instructions and examples using the test data under `test_data/input_data`.

### Step 0: Download and uncompress test data
```bash
wget http://genomebrowser-uploads.hms.harvard.edu/data/yuz006/test_data.tar.gz
tar -xzvf test_data.tar.gz
```

### Step 1: Generate chains

The `chain` command generates liftover chain files for mapping SV coordinates across reference assemblies.  
Before running it, please create a ref2-to-ref1 alignment using [minimap2](https://github.com/lh3/minimap2)'s asm-to-asm mapping:

```bash
minimap2 -a -x asm5 --cs <ref1.fa> <ref2.fa> | samtools view -hSb - | samtools sort -O BAM -o ref2_to_ref1.bam
```

```bash
# Generate chain file between two references
echosv chain -b test_data/input_data/chm13_to_grch38.bam -f test_data/input_data/chm13.fa -o test_data/chm13_to_grch38.chain.gz 
```

**Parameters:**
- `-b`: Path to the ref2-to-ref1 alignment (BAM format)
- `-f`: Path to the ref2 reference genome (FASTA format)
- `-o`: Output chain file for coordinate mapping (and bed file for alignment coverage)

### Step 2: Collect supporting reads

The `genotype` command collects supporting reads for each SV in the corresponding BAMs and prepares them for matching by graph-based matching.

```bash
# Genotype SVs from the first call set
echosv genotype --longread -i test_data/input_data/grch38_colo829_somatic_svs.vcf.gz -b BAM [BAMs...] -o test_data/grch38_colo829_genotyped.vcf.gz
```

**Parameters:**
- `--longread`: Genotyping SVs from long-read alignments
- `--shortread`: Genotyping SVs from short-read alignments
- `-i`: Input SV vcf file
- `-b`: Bam file(s)
- `-o`: Output VCF with genotyped info

### Step 3: Match SVs

The `match` command compares SV call sets and either merges them or identifies overlapping and reference-unique variants.

```bash
# Compare and merge two SV call sets
echosv match -i test_data/test_colo829_config.json --merge

# Or compare to find unique and shared variants
echosv match -i test_data/test_colo829_config.json 
```

**Parameters:**
- `-i`: Input config file, see example as ./test_data/test_colo829_config.json
- `--merge`: If given, merge concordant SVs across references and derive a single VCF.
- `--min_echo_score`: Minimum echo score to consider an SV for matching (default: 0.5).

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.

## Contact
Feel free to open an issue in Github or contact Yuwei Zhang ([yuwei_zhang@hms.harvard.edu](mailto:yuwei_zhang@hms.harvard.edu)) if you have any problem in using EchoSV.