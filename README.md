# COmapper

**COmapper: High-resolution mapping of meiotic crossovers by long-read sequencing in _Arabidopsis thaliana_**

### Contributions

- The project was initiated by Prof. Kyuha Choi, and the code was written by Dohwan Byun at the Plant Genomic Recombination (PGR) Laboratory of Pohang University of Science and Technology (POSTECH).
- This research was conducted in collaboration with the following individuals: Dohwan Byun¹, Namil Son¹, Heejin Kim¹, Jaeil Kim¹, Jihye Park¹, Sang-jun Park¹, Hyein Kim¹, Seula Lee², Youbong Hyun², Piotr A. Ziolkowski³, Ian R. Henderson⁴, and Kyuha Choi¹†.

  ¹Department of Biological Sciences, Pohang University of Science and Technology, Pohang, Gyeongbuk, Republic of Korea

  ²School of Biological Sciences, Seoul National University, Seoul, Republic of Korea

  ³Laboratory of Genome Biology, Institute of Molecular Biology and Biotechnology, Adam Mickiewicz University, Poznań, Poland

  ⁴Department of Plant Sciences, University of Cambridge, Cambridge, UK

- The manuscript is available at

## System requirements

linux machine with more than 2TB storage space, 6-core cpu, 32Gb ram

## Getting started

First, clone this git repository.

    $ git clone https://github.com/KyuhaChoi-Lab/COmapper.git
    
### Setting conda environment

**Prerequisites**
- minimap2 v2.28 (https://github.com/lh3/minimap2)
- samtools v1.20 (https://www.htslib.org/)
- sambamba v.1.0.1 (https://lomereiter.github.io/sambamba/)
- python v3.9.19
- pandas v1.1.3
- typing-extensions v.4.12.2

We recommend setting up conda environment to manage the dependencies for COmapper.
You can use the provided environment.yaml to create conda environment with all the required packagers.

```
    $ conda env create --name COmapper --file environment.yaml
```

### Preparing index file of TAIR10 for minimap2
1. Download TAIR10 reference genome sequence from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.3/
2. Create index file

```
    $ minimap2 -d TAIR10.mmi TAIR10.fasta
```

### Explanation of the files in resources directory

  `masksam.awk`: script for basefiltering
    
  `TAIR10.mmi`: minimap index file

  `coller_marker_v4.6.tsv`: SNP list file

    
### Input data
COmapper requires Nanopore sequencing reads in `.fq.gz` (or `.fq`) format.

A test input file is provided in the `raw` directory (`test.fq.gz`).

## Working procedure

1.  Navigate to the home directory of the cloned Git repository.

2.	Place the Nanopore sequencing results (in fq.gz format) into the `raw` directory.

3.	Open the `variables.sh` file and configure the following variables:
    - `dir_home`: The home directory of the cloned repository. Ensure at least 2TB of free space.
  	- `n_parallel`: The number of threads to use.
  	- `input`: The name of the `.fq.gz` input file (e.g., `test.fq.gz`).
    - `awkloc`: location of `masksam.awk` file.
    - `ref`: location of the TAIR10 minimap index file (`TAIR10.mmi`).

3. Activte the conda environment

```
    $ conda activate COmapper
```

4. Run the `COmapper_data_process_ver5.sh` script:

```
    $ source variables.sh
    $ bash COmapper_data_process_ver5.sh
```

5. The output of `COmapper_data_process_ver5.sh` is stored in the `output/tsv` directory. At this stage, high-quality Nanopore reads are retained, and low-quality bases are masked. The next step is to detect crossover molecules from the .tsv files using `COmapper_ver1.2.5_git_upload.py`.

6. The `COmapper_ver1.2.5_git_upload.py` script requires four options

   --inputfolder (-i) : input folder location (including .tsv files), default: current location
   
   --snpfile (-s) : snp file location, default: current location/collerF2.masked.tiger.txt

   --threads (-t) : thread number, default: 6

   --output (-o) : output file name, default: result.csv

    You can set these options via command line arguments or by editing the Python script directly. Then, run `COmapper_ver1.2.5_git_upload.py`:

```
       $ python COmapper_ver1.2.5_git_upload.py -i INPUTFOLDER -s SNPFILE -t THREADS -o OUTPUTFILE
```

8. COmapper will return genotyping and classification information via standard output and crossover molecule information in a .csv file.
