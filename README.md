# INITIAL DEVELOPMENT RELEASE

For questions, bugs, and suggestions, please contact bolukolu@utk.edu.

<p align="right">
<img src="https://github.com/bodeolukolu/ngsComposer/blob/master/misc/logo_ngsComposer.png">
</p>


# ngsComposer: empirically different

Base-call error-filtering and read preprocessing pipeline designed by biologists

## Features

- Full start-to-finish pipeline for many library types
- Few dependencies (Python3 and R)
- Easy to learn
- Supports variable length barcodes and dual-indexing
- Trims buffer sequences and quality filters on a read-by-read basis
- Accepts project directory of multiple libraries
- Designed by biologists (please don't run away!)

## Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Demultiplexing](#demultiplexing)
  - [Standalone Tools](#standalone)
- [Troubleshooting](#troubleshooting)
- [Versioning](#versioning)
- [License](#license)


## Installation

Currently, ngsComposer is only available for unix-based systems (i.e. macOS and linux).

Clone or download the Git repository to your desired folder

```bash
git clone https://github.com/bodeolukolu/ngsComposer.git
```

Dependencies:
- Python3 version 3.5 or above
- R, R-ggplot2
- pigz (not required, but recommended for parallel gzip and gunzip i.e. faster)

For help troubleshooting installation, see the troubleshooting section

## Usage

### Basic usage

Set up your project directory containing the following:
- a folder named <samples>, which contains fastq file(s). Multiple libraries or demultiplexed fastq files can be included.
- file named <barcodes_lib1.txt>, which contains barcodes and associated sample IDs. Additional library fastq files can be included, e.g. second library will have <barcodes_lib2.txt> file corresponding to lib2_R1/lib2_R2 variable (holds fastq file names) in the config.sh.
- 2 files containing adapter sequences of R1/P5/forward (<adapters.R1.txt>) and R2/P7/reverse (<adapters.R2.txt>) reads.
- conf.sh (see "Configuration" below for detailed instructions on creating this file)
                                     .

<a><img src="https://imgur.com/SYiOxfR.png" title="source: imgur.com" width=400 /></a>


From command line, run ngsComposer with the specified directory of your project
```bash
$ bash <path_to_ngsComposer_directory>/ngsComposer <path_to_project_directory>
```

If this is the first time running the pipeline, you may need to wait for R to install the appropriate packages and dependencies.

Several **example datasets** are included in the "examples" directory. Users are encouraged to examine and run these small projects to assist in understanding pipeline functionality.

***

### Overview
The order of steps in the ngsComposer pipeline are outlined in the following figure:

<a><img src="https://i.imgur.com/99BsbsJ.png" title="source: imgur.com" /></a>

The steps implemented are first specified in a configuration file.

***

### Configuration
Using a text editor, save a file containing any of the following variables as a bash script called 'conf.sh' (includes '.sh' as file extension) and include it in your project directory.

**General parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|threads|2|choose maximum number of subprocesses that can run simultaneously|integer|optional|
|walkaway|True|run from beginning to end without pausing at qc steps |True or False|optional|
|cluster|False|run on compute cluster node (default: slurm) or workstation|True or False|optional|
|samples_alt_dir|False|input files stored in different different from project directory|True or False|optional|
|rm_transit|True|remove each transitional file folder to save space|True or False|optional|



**Input files**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|lib1_R1|na|input fastq file name for R1/P5/forward reads|string|required|
|lib1_R2|na|input fastq file name for R2/P7/reverse reads|string|optional|
|lib1_bc|na|name of file containing barcodes|string|required|
|lib2_R1|na|additional input fastq file name for R1/P5/forward reads|string|required|
|lib2_R2|na|additional input fastq file name for R2/P7/reverse reads|string|optional|
|lib2_bc|na|additional name of file containing barcodes|string|required|


**Tool Parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|front_trim|0|number of bases in buffer sequence to trim|integer|optional|
|mismatch|1|number of mismatches (hamming distance) allowed in barcodes|integer|optional|
|R1_motif|na|motif filtering for R1/P5/forward reads|string or list of comma-separated strings|optional|
|R2_motif|na|motif filtering for R2/P7/reverse reads|string or list of comma-separated strings|optional|
|end_score|20|end-trim once entire window >= this Q score|integer|optional|
|window|10|size of window to test for >= end_trim|integer|optional|
|min_len|64|minimum read length to retain after end-trimming and adapter removal|integer|optional|
|adapter_match|12|number of base matches to identify adapters|integer|optional|
|q_min|20|Q score minimum (Phred value 0-40) applied to q_percent variable|integer|optional|
|q_percent|80|percentage of basses in read >= q_min Q scores|integer|optional|

**Visualizations**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|QC_all|na||summary and/or full|optional|
|QC_demultiplexed|na||summary and/or full|optional|
|QC_motif_validated|na||summary and/or full|optional|
|QC_end_trimmed|na||summary and/or full|optional|
|QC_adapter_removed|na||summary and/or full|optional|
|QC_final|summary||summary and/or full|optional|


**Note: na indicates that variable is user-defined or hard-coded/computed intuitively, as well as a function of ploidy.*



An example configuration file may look like this:

**config.sh**

```
#General_parameters
###################################################
threads=24
walkaway=True
cluster=True
samples_alt_dir=False
rm_transit=True

#Input_files
###################################################
lib1_R1=test1_R1.fastq.gz
lib1_R2=test1_R2.fastq.gz
lib1_bc=barcodes_lib1.txt
lib2_R1=test2_R1.fastq.gz
lib2_R2=test2_R2.fastq.gz
lib2_bc=barcodes_lib2.txt

#Tool_parameters
###################################################
front_trim=6
mismatch=1
R1_motif=TGCATA,TGCATC,TGCATT
R2_motif=CATG
end_score=20
window=10
min_len=64
adapter_match=12
q_min=20
q_percent=80

#Visualizations
###################################################
QC_all=summary,full
QC_demultiplexed=summary,full
QC_motif_validated=summary,full
QC_end_trimmed=summary,full
QC_adapter_removed=summary,full
QC_final=summary,full
```

*In the above example, the maximum number of subprocesses spawned will be 24 (**threads = 24**).  The pipeline will pause after relevant steps (**walkaway = False**) so users can view qc plots and have the option of modifying or bypassing the step. To save disk space, transitional directories will be removed (**rm_transit = True**) and only the final filtered data and any qc stats created in the pipeline will remain. Regardless of if walkaway if True or False, pipeline will ask if initial_qc should be generated.*

*A buffer sequence of length 6 (**front_trim = 6**) will be trimmed before demultiplexing, which will allow mismatch at a hamming distance of 1 (**mismatch=1**). For variable length barcodes, the same number of proximal bases (based on the minimum barcode length) are used for demultiplexing, while the additional distal bases in barcodes are trimmed off.*

*In this case, samples were double-digested with AluI and HaeIII and A-tailed before adapter ligation (**R1_motif=TCC,TCT** and **R2_motif=TCC,TCT**). Only reads containing these motifs will pass to subsequent steps. For A-tailed libraries, an A can be appended to the R1_motif and R2_motif strings.*

*Automatic end-trimming will be performed based on Q score. Here, groups of bases are considered within a moving window of 10 bases at a time (**window=10**) until that window consists only of the desired Q score at or above 20 (**end_score=20**). It is at this point that the read is trimmed. Reads that are less than 64 bp will be discarded (**min_len=64**)*

*Only reads that have a Q score of 20 (**q_min=20**) acrosss at least 95 percent of the read (**q_percent=80**) will pass to subsequent steps. If a R1 read or an R2 read passes while its partner fails, it will be placed into a single-end read subfolder and the failing read will be discarded.*


Alternatively, a configuration file may only need to include necessary components for a run:

**conf.py**

```
#Input_files
###################################################
lib1_R1=test1_R1.fastq.gz
lib1_bc=barcodes_lib1.txt
```
*Since most of the parameters are hard-coded in an intuitive manner, by specifying only the fastq file name (at least single-end data) and associated barcode (only required for demultiplexing), the pipelines determines the other parameters as some stated below:<br />
- threads: computes available number of cores (n) and uses n-2 threads
- defaults: walkaway=True, cluster=False, samples_alt_dir=False, rm_transit=True, front_trim=0, mismatch=1, no motif filtering, end_score=20, min_len=64, adapter_match=12, q_min=20, q_percent=80, only summary final QC, and initial QC will be determined based on a prompt before submitting job.

***

### Demultiplexing
#### Barcodes file(s)
Optionally, one or more barcode files may be included in the project directory for demultiplexing. The following files are required at minimum:
- barcodes_1.txt
- index.txt

<a><img src="https://imgur.com/VlLCeY4.png" title="source: imgur.com" width=400 /></a>

Naming conventions: "index.txt" is required, the barcodes file can be named as desired (see "Index file for directing multiple barcode files")

The barcodes file is a tab or space delimited file with no spaces in sample names (or, copy directly from your favorite spreadsheet program into a text file). Forward barcodes begin each row and reverse barcodes begin each column with the desired sample names indicated in the interior of the matrix. For example, the following would be required for a dual-indexed library:

**barcodes_1.txt**
```
	A	C	G	T
A	sample1	sample5	sample6	sample10
C	sample2	sample5	sample7	sample10
G	sample3	sample5	sample8	sample10
T	sample4	sample5	sample9	sample10
```
*Note that in the example above the reverse barcode "C" corresponds with multiple identical sample names (sample5). While not common practice, ngsComposer accomodates repeated sample names and concatenates accordingly.*

If reverse barcodes do not require demultiplexing, the barcode file can be set up as follows with "NA" or any other text used as a header in the first row:

**barcodes_1.txt**
```
	NA
A	sample1
C	sample2
G	sample3
T	sample4
```

#### Index file for directing multiple barcode files
The index file is a tab delimited file required to associate the barcodes file with a specific library in your project directory. It must include the filename of the forward read (R1) followed by the appropriate barcodes file. Reverse reads (R2), if present, will automatically be detected and are not indicated in this file.

**index.txt**
```
1_R1.fastq  barcodes_1.txt
```

Alternatively, multiple barcoding schemes may be included to accomodate multiple libraries. For example:

**index.txt**
```
1_R1.fastq  1_bcs.txt
2_R1.fastq  2_bcs.txt
```

<a><img src="https://imgur.com/mnrrviL.png" title="source: imgur.com" width=600 /></a>

*In this example, sample "1_R1.fastq" and "1_R2.fastq" correspond with "1_bcs.txt" and "2_R1.fastq" and "2_R2.fastq" correspond with "2_bcs.txt"*

***

### Adapters
#### Adapters file(s)
Optionally, 'adapters.R2.txt' and 'adapters.R1.txt' may be included in the project directory for recognition and removal of adapters. The 'adapters.R2.txt' file contains the adapters expected to appear in the R1 reads. Adapter sequences should be newline-separated and be in 5' to 3' orientation. If libraries are barcoded, users are encouraged to provide adapter sequences that contain the corresponding barcodes expected in the opposing end of the read's adapter.

<a><img src="https://i.imgur.com/3In8TX0.png" title="source: imgur.com" width=300 /></a>


**adapters.R2.txt**
```
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGCTCAGTTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTATCTGACCT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTATATGAGACG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTTATGGAAT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTAATCTCGTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCGCGATGTT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAGAGCACTAG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCCTTGATC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTACTCAGTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCGTCTGACT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGAACATACGG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCTATGACTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTAATGGCAAG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTGCCGCTTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGGCAATGGA
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCCGTAACCG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAACCATTCTC
```
*Each of the above sample adapters is presented in 5' to 3' orientation and shares a common 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' adapter sequence followed by expected barcodes. Adapter sequences may also include restriction motifs for greater detection, but these sequences will also be removed. Porifera.py creates all reverse-complements before detection.*

<a><img src="https://i.imgur.com/mYBQWIv.png" title="source: imgur.com" width=600 /></a>

*When paired end data is used, as above, 'adapters.R1.txt' and 'adapters.R2.txt' must be provided. Adapters are tested for the inclusion of barcodes and only those combinations of R1/R2 barcodes leading to a given sample will be used to search for adapters quickly and with a lower false positive rate.*

### Standalone
All tools available in the ngsComposer pipeline can be called individually from the command line. Please see the <a href="https://github.com/bodeolukolu/ngsComposer/blob/master/tools/README.md">ngsComposer Standalone Tools page</a> for usage.


## Troubleshooting

### Python installation
To view Python version, from the terminal type:

```bash
$ python3 --version
```

If python3 is not found, you can try one of the python3 releases from the Python Software Foundation <a href="https://www.python.org/downloads/">downloads page</a>.

Alternatively, a package manager is an easy way to install Python from the terminal. For Ubuntu, Python can be installed directly using apt (replace 'X' with an existing version in the apt repository):

```bash
$ sudo apt-get update
$ sudo apt-get install python3.X
```

...or with homebrew on macOS using:

```
brew install python3
```

After installation please check that the newest version is present in your current environment (i.e.; $PATH).

### R installation
To view R version, from the terminal type:

```bash
$ R --version
```

To install the newest version of R, see the releases available at the Comprehensive R Archive Network <a href="https://cran.r-project.org/">downloads page</a>.

For Ubuntu, R can be installed directly using apt:

```bash
$ sudo apt update
$ sudo apt install r-base
```

...or with homebrew on macOS using:

```
brew install r
```

**Notes on ggplot2 installation:**

ngsComposer requires the R package ggplot2 and its dependencies. ngsComposer will attempt to automatically download these packages to the local ngsComposer repo (/ngsComposer/tools/helpers/R_packages).

The installation of ggplot2 and dependencies may take some time during the first use.
If package installation fails, manual installation within R may be necessary.
It may be beneficial to install R packages as root and if using macOS, ensure Xcode toolkit is up to date.

## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License

<a href="https://github.com/bodeolukolu/ngsComposer/blob/master/misc/LICENSE">Apache License Version 2.0</a>
