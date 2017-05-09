#!/usr/bin/env python


"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgreer2@stanford.edu
:Creation date: 24.11.2016
:Description: 

This script annotates the structural variant output of longranger (10X Genomics)

This script requires:
- all of the python packages listed (imported) below
	
Revisions:
None to date

CURRENT VERSION: 1.0

"""

cur_version = 1.0

### LOAD THE NECESSARY PACKAGES ###

import sys
import os
import __main__ as main
import argparse

import pandas as pd
import pybedtools
from pybedtools import BedTool
from collections import defaultdict

pd.options.mode.chained_assignment = None


#################################################################
################                                 ################
################        PARSE THE ARGUMENTS      ################
################                                 ################
#################################################################

### ARGPARSE PARSING ###

def usage():
	print "Usage examples:"
	print os.path.basename(main.__file__) + " --help"
	print os.path.basename(main.__file__) + " -v longranger_svs_tumor.bedpe -n longranger_svs_normal.bedpe -l lumpy_svs.vcf -b bicseq_cnvs.txt -g /path/to/file/of/genes.txt -p 5000 -q 1000000 -out output_prefix"
	sys.exit(0)

def parse_args():
	parser = argparse.ArgumentParser(description = "A Python script for annotating the SVs called by longranger")
	parser.add_argument("--usage", help="usage example", dest="usage", action='store_true')
	parser.add_argument("-v", help="BEDPE file of tumor SVs called by longranger (REQUIRED)", dest="lr_tum_in")
	parser.add_argument("-n", help="BEDPE file of normal SVs called by longranger (REQUIRED)", dest="lr_norm_in")
	parser.add_argument("-l", help="VCF file of SVs called by lumpy (REQUIRED)", dest="lmpy_in")
	parser.add_argument("-b", help="BED file of CNVs called by BICseq (REQUIRED)", dest="bic_in")
	parser.add_argument("-g", help="BED file of genes of interest (REQUIRED)", dest="gene_in")
	parser.add_argument("-p", help="base pairs of padding around regions for intersection", dest="padding", default=0)
	parser.add_argument("-q", help="base pairs of padding around regions for gene intersection", dest="g_padding", default=0)
	parser.add_argument("-out", help="prefix for output files (REQUIRED)", dest="out_in")
	parser.add_argument("--version", action='version', version='%(prog)s ' + str(cur_version))	
	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()
	if(args.usage):
		usage()
	if(not args.lr_tum_in or not args.lr_norm_in or not args.lmpy_in or not args.bic_in or not args.out_in):
		print os.path.basename(main.__file__) + " missing a required input file\n"
		usage()
		sys.exit(1)

### SET THE ARGUMENTS ###

tum_sv = args.lr_tum_in
norm_sv = args.lr_norm_in
lumpy_vcf = args.lmpy_in
bic_file = args.bic_in
gene_file = args.gene_in
pad_bp = args.padding
gpad_bp = args.g_padding
out_prefix = str(args.out_in)


#################################################################
################                                 ################
################         DEFINE FUNCTIONS        ################
################                                 ################
#################################################################

## Function to parse data from SV bedpe files -- create df for all breakpoint 1's and df for all breakpoint 2's, then stack df's and pad coordinates

def bedpe_parse(d, pad):
	df_1 = d[['#chrom1','start1','stop1','name']]
	df_1.columns = ['#chrom','start','stop','name']
	df_1['bkpt'] = '1'
	df_2 = d[['chrom2','start2','stop2','name']]
	df_2.columns = ['#chrom','start','stop','name']
	df_2['bkpt'] = '2'
	df_both = pd.concat([df_1, df_2], ignore_index=True)
	df_both['start_pad'] = df_both['start'].apply(lambda x: 0 if x-int(pad)<0 else x-int(pad))
	df_both['stop_pad'] = df_both['stop'] + int(pad)
	cols = df_both.columns.tolist()
	cols = cols[:1] + cols[-2:] + cols[1:-2]
	df_both = df_both[cols]
	return df_both

### Function to turn INFO field of vcf file into a dictionary
def makeDictInfo(r):
	d = defaultdict(dict)
	info_list=r.split(';')
	for i in range(len(info_list)):
		if '=' in info_list[i]:
			key, value = info_list[i].split("=")
			d[key]=value
	return d

### Function to turn FORMAT/SAMPLE fields of vcf file into a dictionary
def makeDictSample(r):
	format_list=r['format'].split(':')
	sample_list=r['sample'].split(':')
	d = dict(zip(format_list,sample_list))
	return d

### Function to generate end position 
def getEnd(r):
	cur_dict = r['info_dict']
	if "END" in cur_dict:
		return cur_dict['END']
	else:
		end_pos = 'na'
		return end_pos


#################################################################
################                                 ################
################       DETERMINE SOMATIC SVs     ################
################                                 ################
#################################################################

### PARSE NORMAL SV FILE (i.e. put each breakpoint on its own line + pad the coordinates)

df_norm = pd.read_table(norm_sv, sep="\t", skiprows=1)
df_norm_both = bedpe_parse(df_norm, pad_bp)
header = list(df_norm_both.columns)
df_norm_both.columns = [h+"_norm" for h in header]


### PARSE TUMOR SV FILE

df_tum = pd.read_table(tum_sv, sep="\t", skiprows=1)
df_tum_both = bedpe_parse(df_tum, pad_bp)


### PERFORM INTERSECTION OF NORMAL + TUMOR FILES

df_norm_both_str = df_norm_both.to_string(index=False)							## Convert df to string (seems to be only way to coerce pandas df into bedtool object)
df_norm_bed = BedTool(df_norm_both_str, from_string=True) 						## Convert df as str to bedtool object

df_tum_both_str = df_tum_both.to_string(index=False) 							## Convert df to string (seems to be only way to coerce pandas df into bedtool object)
df_tum_bed = BedTool(df_tum_both_str, from_string=True) 						## Convert df as str to bedtool object
som_header = list(df_tum_both.columns) + list(df_norm_both.columns) 			## Create header
som_header = [s.strip('#') for s in som_header] 								## Remove comment from header
som_isect = df_tum_bed.intersect(df_norm_bed, wa=True, wb=True) 				## Intersect files
som_isect = som_isect.to_dataframe(names=som_header) 							## Convert result to pandas data frame


### ADD SOMATIC INFO TO ORIGINAL TUMOR SV FILE

## Determine which tumor breakpoints are germline and add info to original file + add a column to indicate if somatic

som_isect = som_isect[['name','bkpt','name_norm','chrom_norm','start_norm','stop_norm']] 															## Extract + reorder columns
som_isect[['name_norm','chrom_norm','start_norm','stop_norm']] = som_isect[['name_norm','chrom_norm','start_norm','stop_norm']].astype(str)					## Change cols to strings in prep for join
som_isect['info_lr_norm(name, chrom, start, end)'] = som_isect[['name_norm','chrom_norm','start_norm','stop_norm']].apply(lambda x: ','.join(x), axis=1)	## Merge columns of normal info
som_isect = som_isect[['name','bkpt','info_lr_norm(name, chrom, start, end)']]																		## Extract columns to retain in final table
som_isect[['bkpt']] = som_isect[['bkpt']].astype(str) 																								## Change col to string in prep for join
som_isect_summ = som_isect.groupby('name').agg({'bkpt': lambda x: ', '.join(x), 'info_lr_norm(name, chrom, start, end)': lambda x: '; '.join(x)}).reset_index() ## Group info for final table 																	## Rename columns

som_merge = pd.merge(df_tum, som_isect_summ, how='left', on='name') 																						## Merge somatic summary table with original longranger calls for tumor
som_merge['som'] = som_merge['bkpt'].apply(lambda x: 'X' if str(x)=='nan' else '') 																	## Add column to indicate if SV event is somatic
som_merge=som_merge.rename(columns = {'bkpt':'brkpt_in_norm'})

### FINAL OUTPUT OF THIS STAGE
som_merge.to_csv(out_prefix + ".som.bedpe", sep="\t", index=False)


#################################################################
################                                 ################
################  DETERMINE SVs CALLED BY LUMPY  ################
################                                 ################
#################################################################

### OPEN LUMPY VCF FILE

df_lmpy = pd.read_table(lumpy_vcf, sep="\t", comment='#', header=None, names=['chr','pos','id','ref','alt','qual','filter','info','format','sample'])

### APPLY DICT FUNCTIONS

df_lmpy['info_dict']=df_lmpy['info'].apply(lambda x: makeDictInfo(x)) 	## Turn info field into dictionary
df_lmpy['sample_dict']=df_lmpy.apply(lambda row: makeDictSample(row), axis=1)	## Turn format/sample fields into dictionary

## EXTRACT VALUES OF INTEREST FROM DICTS (+ CREATE NEW COLUMNS)

df_lmpy['sv_type']=df_lmpy['info_dict'].apply(lambda x: x['SVTYPE'])	## Make a column of SVtype
df_lmpy['end'] = df_lmpy.apply(lambda row: getEnd(row), axis=1)		## Make a column with the 'end value' - NOTE: 'BND's will not have an 'END' value
df_lmpy['ci_vals_start']=df_lmpy['info_dict'].apply(lambda x: x['CIPOS'] if "CIPOS" in x else 'na')	# Make a column with the 'start' adjustment values
df_lmpy['ci_vals_end']=df_lmpy['info_dict'].apply(lambda x: x['CIEND'] if "CIEND" in x else 'na')	# Make a column with the 'end' adjustment values
df_lmpy['total_evid']=df_lmpy['sample_dict'].apply(lambda x: x['SU'])	# Create columns for soft-clip and read pair evidence
df_lmpy['paired_end']=df_lmpy['sample_dict'].apply(lambda x: x['PE'])
df_lmpy['split_read']=df_lmpy['sample_dict'].apply(lambda x: x['SR'])

### SEPARATE BREAKPOINTS INTO OWN DF'S THEN CONCAT TO ONE DF

df_lmpy_1 = df_lmpy[['id','chr','pos','sv_type','total_evid','paired_end','split_read','ci_vals_start']] 	## Get start positions
df_lmpy_2 = df_lmpy[['id','chr','end','sv_type','total_evid','paired_end','split_read','ci_vals_end']]		## Get end positions
df_lmpy_2 = df_lmpy_2.loc[df_lmpy['end']!='na']																	## Remove rows without an end position
df_lmpy_colnames = ['id','chr','pos','sv_type','total_evid','paired_end','split_read','ci_vals']			## Rename columns
df_lmpy_1.columns = df_lmpy_colnames
df_lmpy_2.columns = df_lmpy_colnames
df_lmpy_both=pd.concat([df_lmpy_1,df_lmpy_2])																	## Then concatenate

### GENERATE COORDINATES OF CONFIDENCE INTERVALS AROUND POSITIONS

df_lmpy_both['ci_val_1']=df_lmpy_both['ci_vals'].apply(lambda x: x.split(",")[0])								## Get first value
df_lmpy_both['ci_val_2']=df_lmpy_both['ci_vals'].apply(lambda x: x.split(",")[1])								## Get second value
df_lmpy_both[['pos','ci_val_1','ci_val_2']]=df_lmpy_both[['pos','ci_val_1','ci_val_2']].astype(int)				
df_lmpy_both['ci_start'] = df_lmpy_both['pos'] + df_lmpy_both['ci_val_1']												## Get start position for SV window
df_lmpy_both['ci_end'] = df_lmpy_both['pos'] + df_lmpy_both['ci_val_2']													## Get end position for SV window

## Pad the start and stop positions

df_lmpy_both['ci_start'] = df_lmpy_both['ci_start'].apply(lambda x: 0 if x-int(pad_bp)<0 else x-int(pad_bp))
df_lmpy_both['ci_stop'] = df_lmpy_both['ci_end'] + int(pad_bp)

### FORMAT THE LUMPY DATA IN PREP FOR INTERSECTION

df_lmpy_both[['chr','id','pos','sv_type','paired_end','split_read']] = df_lmpy_both[['chr','id','pos','sv_type','paired_end','split_read']].astype(str)	## Change cols to strings
df_lmpy_both['info_lumpy(id,pos,sv_type,paired_end,split_read)'] = df_lmpy_both[['id','pos','sv_type','paired_end','split_read']].apply(lambda x: ','.join(x), axis=1)	## Merge columns
df_lmpy_both = df_lmpy_both[["chr","ci_start","ci_end","info_lumpy(id,pos,sv_type,paired_end,split_read)"]]	## Extract columns
df_lmpy_both.sort_values(by=['chr','ci_start'], inplace=True)
df_lmpy_both=df_lmpy_both.rename(columns = {'chr':'#chr'})

### TURN LUMPY DF INTO BEDTOOLS OBJECT

lmpy_str = df_lmpy_both.to_string(index=False) # convert df to string (seems to be only way to coerce pandas df into bedtool object)
lmpy_bed = BedTool(lmpy_str, from_string=True) # convert df as str to bedtool object

### INTERSECT WITH SV DATA FRAME

lmpy_isect = df_tum_bed.intersect(lmpy_bed, wa=True, wb=True) # intersect files
lmpy_header = list(df_tum_both.columns) + list(df_lmpy_both.columns) # create header
lmpy_header = [s.strip('#') for s in lmpy_header] # add header
lmpy_isect = lmpy_isect.to_dataframe(names=lmpy_header) # convert result to pandas data frame


### FORMAT OUTPUT

lmpy_isect = lmpy_isect[['name','bkpt','info_lumpy(id,pos,sv_type,paired_end,split_read)']]
lmpy_isect[['bkpt']] = lmpy_isect[['bkpt']].astype(str)
lmpy_isect_summ = lmpy_isect.groupby('name').agg({'bkpt': lambda x: ', '.join(x), 'info_lumpy(id,pos,sv_type,paired_end,split_read)': lambda x: '; '.join(x)}).reset_index()
lmpy_isect_summ = lmpy_isect_summ.rename(columns = {'bkpt':'bkpt_lumpy'})

## Annotate longranger calls with presence in lumpy output

lmpy_merge = pd.merge(som_merge, lmpy_isect_summ, how='left', on='name')
lmpy_merge['lumpy'] = lmpy_merge['bkpt_lumpy'].apply(lambda x: '' if str(x)=='nan' else 'X')

### FINAL OUTPUT
#print lmpy_merge
lmpy_merge.to_csv(out_prefix + ".som.lmpy.bedpe", sep="\t", index=False)

#################################################################
################                                 ################
################  DETERMINE SVs CALLED BY BICSEQ ################
################                                 ################
#################################################################

### PARSE BICSEQ OUTPUT

df_bic = pd.read_table(bic_file, sep="\t")
df_bic['cov_ratio'] = df_bic.apply(lambda row: (2**row['log2.copyRatio']), axis=1)
df_bic_sub = df_bic.loc[(df_bic['cov_ratio']<=0.95) | (df_bic['cov_ratio']>=1.25)]
df_bic_sub['id'] = df_bic_sub.index

df_bic_1 = df_bic_sub[['chrom','start','id','cov_ratio']]
df_bic_2 = df_bic_sub[['chrom','end','id','cov_ratio']]
df_bic_colnames = ['chrom','pos','id','cov_ratio']
df_bic_1.columns = df_bic_colnames
df_bic_2.columns = df_bic_colnames
df_bic_both = pd.concat([df_bic_1,df_bic_2], ignore_index=True)

df_bic_both['start_pad'] = df_bic_both['pos'].apply(lambda x: 0 if x-int(pad_bp)<0 else x-int(pad_bp))
df_bic_both['stop_pad'] = df_bic_both['pos'] + int(pad_bp)
df_bic_both = df_bic_both[['chrom','start_pad','stop_pad','id','pos','cov_ratio']]
df_bic_both[['chrom','start_pad','stop_pad','id','pos','cov_ratio']] = df_bic_both[['chrom','start_pad','stop_pad','id','pos','cov_ratio']].astype(str)
df_bic_both['info_bic(id,pos,cov_ratio)'] = df_bic_both[['id','pos','cov_ratio']].apply(lambda x: ','.join(x), axis=1)

df_bic_both = df_bic_both[['chrom','start_pad','stop_pad','info_bic(id,pos,cov_ratio)']]
df_bic_both = df_bic_both.rename(columns = {'chrom':'#chrom'})

### TURN BICSEQ DF INTO BEDTOOLS OBJECT

bic_str = df_bic_both.to_string(index=False) # convert df to string (seems to be only way to coerce pandas df into bedtool object)
bic_bed = BedTool(bic_str, from_string=True) # convert df as str to bedtool object

### INTERSECT

bic_isect = df_tum_bed.intersect(bic_bed, wa=True, wb=True) # intersect files
bic_header = list(df_tum_both.columns) + list(df_bic_both.columns) # create header
bic_header = [s.strip('#') for s in bic_header] # add header
bic_isect = bic_isect.to_dataframe(names=bic_header) # convert result to pandas data frame
#print bic_isect

### FORMAT OUTPUT

bic_isect = bic_isect[['name','bkpt','info_bic(id,pos,cov_ratio)']]
bic_isect[['bkpt']] = bic_isect[['bkpt']].astype(str)
bic_isect_summ = bic_isect.groupby('name').agg({'bkpt': lambda x: ', '.join(x), 'info_bic(id,pos,cov_ratio)': lambda x: '; '.join(x)}).reset_index()
bic_isect_summ = bic_isect_summ.rename(columns = {'bkpt':'bkpt_bic'})

## Annotate longranger calls with presence in lumpy output

bic_merge = pd.merge(lmpy_merge, bic_isect_summ, how='left', on='name')
bic_merge['bicseq'] = bic_merge['bkpt_bic'].apply(lambda x: '' if str(x)=='nan' else 'X')

### FINAL OUTPUT
#print bic_merge
bic_merge.to_csv(out_prefix + ".som.lmpy.bic.bedpe", sep="\t", index=False)


#################################################################
################                                 ################
################  DETERMINE NEARBY DRIVER GENES  ################
################                                 ################
#################################################################

#genome_file = "/home/sgreer2/tcga_stad_merge/A02_make_bed/goi.pad.bed"

pad_amt = gpad_bp

## Split the data frame into interchromosomal and intrachromosomal events
df_tum['inter'] = df_tum.apply(lambda row: "True" if row['#chrom1']!=row['chrom2'] else "False", axis=1)

df_inter = df_tum.loc[df_tum['inter'] == "True"]
df_intra = df_tum.loc[df_tum['inter'] == "False"]

## Create bed-style df for intrachromosomal events

df_intra['#chr'] = df_intra['#chrom1']
df_intra['start'] = df_intra['start1'].apply(lambda x: 0 if x-int(gpad_bp)<0 else x-int(gpad_bp))
df_intra['end'] = df_intra['stop2'] + int(gpad_bp)

df_intra = df_intra[['#chr','start','end','#chrom1','start1','stop1','chrom2','start2','stop2','name', 'info']]

## Create bed-style df for interchromosomal events

df_inter2 = df_inter.copy()

df_inter['#chr'] = df_inter['#chrom1']
df_inter['start'] = df_inter['start1'].apply(lambda x: 0 if x-int(gpad_bp)<0 else x-int(gpad_bp))
df_inter['end'] = df_inter['stop1'] + int(gpad_bp)
df_inter = df_inter[['#chr','start','end','#chrom1','start1','stop1','chrom2','start2','stop2','name', 'info']]

df_inter2['#chr']	= df_inter['chrom2']
df_inter2['start'] = df_inter['start2'].apply(lambda x: 0 if x-int(gpad_bp)<0 else x-int(gpad_bp))
df_inter2['end']	= df_inter['stop2'] + int(gpad_bp)
df_inter2 = df_inter2[['#chr','start','end','#chrom1','start1','stop1','chrom2','start2','stop2','name','info']]

df_sv = pd.concat([df_intra, df_inter, df_inter2])

sv_str = df_sv.to_string(index=False) # convert df to string (seems to be only way to coerce pandas df into bedtool object)
sv_bed = BedTool(sv_str, from_string=True) # convert df as str to bedtool object

genes_bed = BedTool(gene_file) # convert other file to bedtool object

genes_isect = sv_bed.intersect(genes_bed, wa=True, wb=True) # intersect files
genes_header = list(df_sv.columns) + ['gene_chr','gene_start','gene_end','gene']
genes_isect = genes_isect.to_dataframe(names=genes_header) # convert result to pandas data frame

genes_isect = genes_isect[['name','gene']]
genes_isect_summ = genes_isect.groupby('name')['gene'].apply(lambda x: ','.join(x)).reset_index()

## Annotate longranger calls with presence in lumpy output

genes_merge = pd.merge(bic_merge, genes_isect_summ, how='left', on='name')

### FINAL OUTPUT
genes_merge.to_csv(out_prefix + ".som.lmpy.bic.genes.bedpe", sep="\t", index=False)
