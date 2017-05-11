# sv-phasing_linkedreads
This repository contains scripts to phase SVs using linked read WGS data from the 10X Genomics platform.

For all of the scripts, please see the import packages section for package dependencies. 

######################################################################################
FIRST SCRIPT TO RUN:
phase_svs.py: This script phases the structural variant output of longranger (10X Genomics)

If you have run Long Ranger on linked read sequencing data then you have all of the files required to run this script. 

Input:
- VCF file generated by Long Ranger for tumor sample
- VCF file generated by Long Ranger for normal sample (if no normal, then can use tumor VCF here as well)
- File of structural variants called by Long Ranger (bedpe)
- Bam file generated by Long Ranger for tumor sample
- Window size to consider surrounding each SV breakpoint (i.e. barcodes will be considered within this window) -- 100 kb is default and recommended

Output: 
- *.sv_bcs.txt
  - summarizes the number of barcodes in each breakpoint window, the number of barcodes shared between breakpoints of SV event 
  (i.e. SV-specific barcodes), and the number of SV-specific barcodes shared between SV events 
- *.sv_haps.txt
  - for each breakpoint, all assiated phase_ids are reported (i.e. within specified window) along with the barcode SNV overlap 
  for haplotype 1 and haplotype 2
  
######################################################################################
NEXT SCRIPT TO RUN:
count_bcs_in_windows.py: This script counts the number of unique barcodes in windows around the SV breakpoints

Input:
- Bam file generated by Long Ranger for tumor sample
- One of the output files generated from above script (phase_svs.py): *.sv_bcs.txt (from here, we want the SVs and their   
  SV-specific barcodes)
- Window size to consider around each breakpoint
- Window size to consider within above window area
  - i.e. default is 1000 bp windows within 500,000 bp total window

Output:
- For each SV in the *.sv_bcs.txt file, 2 files will be generated:
  1) *_1.bc_windows.txt
  2) *_2.bc_windows.txt
  Each file contains a row for each tandem 1 kb window across the 500 kb region (i.e. 500 rows), and a column for each 
  SV-specific barcode for the SV; the cells contain either a '1' indicating that a read with that barcode mapped in that 
  window or a '0' indicating that it did not

######################################################################################
TO PLOT THE OUTPUT OF THE ABOVE SCRIPTS:
plot_bcs_across_bkpts.R: This script plots the location of reads with HMW barcodes (using output of: count_bcs_in_windows.py)

NOTE: Before running this script, a small amount of intervention is required. Based on the results generated by phase_svs.py in the file *.sv_haps.txt, generate a file (*_hap_bcs.txt) for each SV with 2 tab-separated columns, that looks like, for ex:
bcs                     haps
GAAACTCCACCACGTG        1
CAACCTCCAAGCCGCT        1
GTCCCATCACAAACTC        1

The first column contains the names of all the barcodes that could be assigned to a haplotype and the second column contains the haplotype assignment (i.e. either 1 or 2). This information will be used to color the resulting plot.

Output:
- *_bkpt_region.pdf: a pdf of a plot where each SV-specific molecule is mapped across the breakpoints 

######################################################################################
IF YOU HAPPEN TO HAVE ALL THE REQUIRED INPUTS...:
filt_svs.py: This script annotates the structural variant output of longranger (10X Genomics)


  


