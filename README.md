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

4. Activte the conda environment

```
    $ conda activate COmapper
```

5. Run the `COmapper_data_process_ver5.sh` script:

```
    $ source variables.sh
    $ bash COmapper_data_process_ver5.sh
```

6. The output of `COmapper_data_process_ver5.sh` is stored in the `output/tsv` directory. At this stage, high-quality Nanopore reads are retained, and low-quality bases are masked. The next step is to detect crossover molecules from the .tsv files using `COmapper_ver1.2.5_git_upload.py`.

7. The `COmapper_ver1.2.5_git_upload.py` script requires four options

   --inputfolder (-i) : input folder location (including .tsv files), default: current location
   
   --snpfile (-s) : snp file location, default: current location/collerF2.masked.tiger.txt

   --threads (-t) : thread number, default: 6

   --output (-o) : output file name, default: result.csv

    You can set these options via command line arguments or by editing the Python script directly. Then, run `COmapper_ver1.2.5_git_upload.py`:

```
       $ python COmapper_ver1.2.5_git_upload.py -i INPUTFOLDER -s SNPFILE -t THREADS -o OUTPUTFILE
```

8. COmapper will return genotyping and classification information via standard output and crossover molecule information in a .csv file.

   .csv file provides information of crossover reads in each rows. Each column provides the following information:
   - read_ID: Unique number for each read, sorted by position
   - chr_num: The number of the chromosomes (e.g. 1,2,3,4,5)
   - pos: The start position of the read in the chromosome.
   - length: The length of the read
   - num_of_snp: The number of SNPs in the read
   - snp_pos: The matrix with the position of the SNPs in the read
   - bases: The matrix with the base of each SNP position
   - types: The string with SNP genotype excluding 'N’ types
   - SNP_genotype: The string with the genotyping result of all SNPs in the read
   - haplotype: The type of read as classified by COmapper.
       1: C-L simple crossover
       2: L-C simple crossover
       11~14: C-L gene conversion associated crossover
       15~18: L-C gene conversion associated crossover
   - CO_info: The information of the crossover site in the read with five elements. [chromosome number, nearest SNP position before crossover site, nearest SNP position after crossover site, mid-point of two nearest SNP positions, width of two nearest SNP positions
output file example:
      ,read_id,chr_num,pos,length,num_of_snp,snp_pos,bases,types,SNP_genotype,haplotype,CO_info
      0,178316,3,15109545,7676,30,"[15113299, 15114197, 15114345, 15114454, 15114543, 15114782, 15114786, 15114795, 15114932, 15114953, 15114962, 15115000, 15115018, 15115074, 15115081, 15115174, 15115207, 15115423, 15115472, 15115724, 15115780, 15116162, 15116318, 15116369, 15116394, 15116704, 15116983, 15117110, 15117127, 15117215]","['C', 'N', 'N', 'N', 'N', 'G', 'G', 'N', 'G', 'N', 'N', 'C', 'N', 'N', 'N', 'A', 'N', 'C', 'N', 'N', 'T', 'N', 'N', 'N', 'C', 'C', 'G', 'N', 'G', 'G']",CCLCCCCCLLLL,CNNNNCLNNNNCNNNCNCNNCNNNCLLNLL,1,"[3, 15116394, 15116704, 15116549.0, 310]"
      1,181082,3,15139535,7267,28,"[15140158, 15140378, 15140612, 15140768, 15141097, 15142546, 15142572, 15142588, 15142597, 15142801, 15142830, 15142868, 15142920, 15142959, 15143104, 15143162, 15143812, 15143837, 15143927, 15144535, 15144797, 15145079, 15145100, 15145113, 15146701, 15146706, 15146741, 15146789]","['G', 'N', 'N', 'N', 'G', 'N', 'G', 'N', 'A', 'N', 'A', 'N', 'G', 'C', 'N', 'N', 'N', 'G', 'T', 'C', 'C', 'N', 'N', 'T', 'G', 'G', 'N', 'C']",LLCLCLLLLCLCCCC,LNNNLNCNLNCNLLNNNLLCLNNCCCNC,11,"[3, 15144535, 15144797, 15144666.0, 262]"



