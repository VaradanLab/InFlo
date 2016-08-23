##InFlo Version 1.0

##### For Academic/Non-Profit use
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>


##### For Non-Academic/For-Profit use:
Please contact us at varadanlab@gmail.com for licensing options. 

##### For Technical Support:
Send us email at varadanlab@gmail.com

==============================================================================================================================
###Reference

#### *Under Review

#### *InFlo: A Novel Systems Biology Framework Identifies cAMP-CREB1 Axis as a Key Modulator of Platinum Resistance in Ovarian Cancer
<p>Nevenka Dimitrova(1), Anil Belur Nagaraj(3), Abolfazl Razi(3), Salendra Singh(3), Sitharthan Kamalakaran(1), Nilanjana Banerjee(1), Peronne Joseph(3), Alexander Mankovich(1), Prateek Mittal(2),#, Analisa DiFeo(3) Vinay Varadan(3)<\p>
===============================================================================================================================
###Introduction
===============================================================================================================================
<p>InFlo is a novel systems biology approach for characterizing a complex cellular processes using a unique multidimensional framework integrating transcriptomic, genomic and/or epigenomic profiles for any given biological sample. InFlo robustly characterizes tissue-specific differences in activities of signaling networks on a genome scale using unique probabilistic models of molecular interactions on a per-sample basis. InFlo is proved to be a robust framework for both discovery of evidence-based biomarkers and therapeutic targets, as well as for facilitating selection of tailored therapies in individual patients.</p>
==============================================================================================================================
###System Requirements
==============================================================================================================================
###Platforms and System Requirements

Hard Disk Space : The module needs around 300 MB of Hard Disk space.

The performance of the tool depends on the processing capabilities of the machine. 

The module performance was tested on two platforms 
  * 2.8 GHz Intel Core i7 processor with 8 GB of RAM iMac machine
  * 32-core  2.10 GHz Intel(R) Xeon(R) CPU E5-2450 with 64 GB linux machine

###Software Architecture

The script was developed on the R platform with version 3.1, and it is assumed that the user also uses the same for running the script. 

###Software Requirements

 * [**R**](http://www.r-project.org/) : Free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms.
<br> preferred Version : >= 3.1 

* [**GCC**](https://gcc.gnu.org/) : The GNU Compiler Collection includes front ends for C, C++, Objective-C, Fortran, Java, Ada, and Go, as well as libraries for these languages (libstdc++, libgcj,...). 
<br> preferred Version : > 4.8

* [**Boost**](http://www.boost.org/) : Boost provides free peer-reviewed portable C++ source libraries.
<br> preferred Version : 1.47 (Provided with the package)

* [**libDAI**](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/) : A free and open source C++ library for Discrete Approximate Inference in graphical models.
<br> preferred Version : 0.2.5 (Provided with the package)

==============================================================================================================================
##Running InFlo
##### Input Files : All the files should be Independent matrices(Genes as rows and samples as columns)
 * Gene Expression File
 * Somatic Copy Number File
 * Methylation Data File

You will need to run the "InFlo.R" script included in the InFlo package.

    $ cd directory/to/InFlo
    $ RScript InFlo.R [-R] [location to InFlo_PROJ_Config.txt]
for help :

    $ RScript InFlo.R [-H]

The R script does the following jobs

 * Read the initial files
 * Check for the duplicate genes, and eliminate redundancy using maximum absolute deviation, and select only genes that are present in all three data types. 
 * Check for sample names and eliminate the samples which are not present in all provided data types. 
 * Run the gene wise Wilcoxon test or DeSeq Test(For RNASeqV2 Data) on the Tumor Gene Expression V/s the Normal Gene Expression and generates probability of each gene being up/down regulated in the individual tumor sample as compared to the control samples
 * Run the gene wise Wilcoxon test on the Tumor Somatic Copy Number V/s the Normal Somatic Copy Number and generates a probability of copy-number alteration (amplification or deletion) in the tumor sample as compared to the control samples on a per-gene basis.
 * Run the gene wise Wilcoxon test on the Tumor Methylation Data  V/s the Normal Methylation data and generates a probability of Methylation events (Methylated or Demethylated) in the tumor sample as compared to the control samples on a per-gene basis
 
 <b>The required condition for above tests is that a Data type should have atleast 3 Normal Samples.</b>

 * Using the results generated by wilcox test and Pathway files, Create Pathways specific *_mRNA.tab files and *_Genome.tab files.   
 * Remove all the files without any information
 * Run the Inflo core on the remaining files.
 * Combine All the result files in a text files consisting of mean Probabalistic scores on a per-tumor sample Basis. 

<p> InFlo generates the probability of activation or inactivation of specific interactions in the pathway file on a per-tumor sample basis by modeling the pathway activity using the gene expression and copy-number data processed</p>

<p>As soon InFlo finishes producing the data user have to use the "Post_Process.R" scripts, which takes up the InFlo results files and Finds out the significant interactions among various Pathways.</p>
================================================================================================================================

###For any queries kindly mail us at varadanlab@gmail.com
