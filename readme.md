### InFlo Version 1.0

##### For Academic/Non-Profit use
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>

##### For Technical Support:
Send us email at varadanlab@gmail.com

________________________________________________________________
### Reference

***Manuscript in press, NPG Oncogene***

#### InFlo: A Novel Systems Biology Framework Identifies cAMP-CREB1 Axis as a Key Modulator of Platinum Resistance in Ovarian Cancer
***Nevenka Dimitrova, Anil Belur Nagaraj, Abolfazl Razi, Salendra Singh, Sitharthan Kamalakaran, Nilanjana Banerjee, Peronne Joseph, Alexander Mankovich, Prateek Mittal, Analisa DiFeo, Vinay Varadan***

### Introduction
________________________________________________________________
<p>InFlo is a novel systems biology approach for characterizing a complex cellular processes using a unique multidimensional framework integrating transcriptomic, genomic and/or epigenomic profiles for any given biological sample. InFlo robustly characterizes tissue-specific differences in activities of signaling networks on a genome scale using unique probabilistic models of molecular interactions on a per-sample basis. InFlo is proved to be a robust framework for both discovery of evidence-based biomarkers and therapeutic targets, as well as for facilitating selection of tailored therapies in individual patients.</p>

### System Requirements
________________________________________________________________

### Platforms and System Requirements

Hard Disk Space : The module needs around 300 MB of Hard Disk space.

The performance of the tool depends on the processing capabilities of the machine. 

The module performance was tested on the following platform

  * 32-core  2.10 GHz Intel(R) Xeon(R) CPU E5-2450 with 64 GB linux machine

### Software Architecture

The script was developed on the R platform with version 3.1, and it is assumed that the user also uses the same for running the script. 

### Software Requirements

 * [**R**](http://www.r-project.org/) : Free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms.
<br> preferred Version : >= 3.1 

* [**GCC**](https://gcc.gnu.org/) : The GNU Compiler Collection includes front ends for C, C++, Objective-C, Fortran, Java, Ada, and Go, as well as libraries for these languages (libstdc++, libgcj,...). 
<br> preferred Version : > 4.8

* [**Boost**](http://www.boost.org/) : Boost provides free peer-reviewed portable C++ source libraries.
<br> preferred Version : 1.47 (Provided with the package)

* [**libDAI**](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/) : A free and open source C++ library for Discrete Approximate Inference in graphical models.
<br> preferred Version : 0.2.5 (Provided with the package)

________________________________________________________________
## Downloading InFlo
InFlo can be downloaded from github by recursively cloning the git repository
   
   $ git clone https://github.com/VaradanLab/InFlo.git
________________________________________________________________
## Preparing Data for the run

##### Input Files : All the files should be Independent matrices(Genes as rows and samples as columns)
 * Gene Expression File
 * Somatic Copy Number File
 * Methylation Data File
#### Users must note that for running the InFlo All the three Input files should have same rows and columns in the same order.

 * SAMPLE_INFORMATION : This File consisting of information if sample is "Tumor" or "Normal"
________________________________________________________________

#### Pathway Files
 InFlo requires the the pathways information to be provided in a tab delimited file. This file consists of the Nodes Information followed  by the Interactions. Kindly check the example Pathways provided the <Home>/pathways.zip File. 

* User will also need to provide a list of Pathways they want to analyse. This file consists of the following columns to be Provided by User. 
  * PID [Numeric]	: Its a unique ID for the pathway
  * FileName [Character] : Filename for the Pathway	
  * Pathway_Name :  Name of the Pathway
  
 Kindly check the example Pathway Information File provided the <Home>/examples_files/TEST_Pathway_Info_Short.txt File.
________________________________________________________________

## Running InFlo


You will need to run the "InFlo.R" script included in the InFlo package.

    $ EXPORT num_of_cores==8 #User have to provide the number of cores available for the process. 
    $ RScript InFlo.R [-R] [location to InFlo_PROJ_Config.txt]
for help :

    $ RScript InFlo.R [-H]

To run InFlo on HPC. Slurm Script should look the following. 
________________________________________________________________
    #SBATCH -J InFlo_Job
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=8
    #SBATCH --output=RD_1ST.out
    #SBATCH --mem=20gb

    module load gcc
    export Number_of_Cores==cpus-per-task
    RScript InFlo.R [-R] [location to InFlo_PROJ_Config.txt]
_________________________________________________________________


The R script does the following jobs

 * Read the initial files
 * Check for the duplicate genes, and eliminate redundancy using maximum absolute deviation, and select only genes that are present in all three data types. 
 * Check for sample names and eliminate the samples which are not present in all provided data types. 
 * Compares gene expression (either microarray or RNASeq), copy-number (SNP 6.0 arrays or whole-exome sequencing), and DNA Methylation profiles of tumors versus controls and generates probability of each gene being up/down regulated in the individual tumor sample as compared to the control samples
 
 <b>The required condition for above tests is that a Data type should have atleast 3 Normal Samples.</b>

 * Using the results generated by wilcox test and Pathway files, Create Pathways specific *_mRNA.tab files and *_Genome.tab files.   
 * Remove all the files without any information
 * Run the Inflo core on the remaining files.
 * Combine All the result files in a text files consisting of mean Probabalistic scores on a per-tumor sample Basis. 

<p> InFlo generates the probability of activation or inactivation of specific interactions in the pathway file on a per-tumor sample basis by modeling the pathway activity using the gene expression and copy-number data processed</p>

<p>The user can do required downstream analysis based on its requirement.</p>
________________________________________________________________

### For any queries kindly mail us at varadanlab@gmail.com
