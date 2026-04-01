import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
#from sklearn.decomposition import PCA
import os
import numpy as np
import re
import subprocess
from upsetplot import UpSet, plot
#help(UpSet)
from matplotlib import pyplot
from collections import Counter
import ast
#import pyclusterprofiler
#import pyupset as pyu
import statistics
from itertools import combinations
import matplotlib.patches as mpatches

import goatools
from goatools.obo_parser import GODag
# from goatools.base import download_go_basic_obo
# from goatools.go_enrichment import GOEnrichmentStudy
# #help(GOEnrichmentStudy)
# from goatools import obo_parser

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
# import gffutils


#--------#
# Inputs #
#--------#
work_dir = "../../output"

#from panTEannot
dir_panTEannot_pop_structure = f"{work_dir}/panTEannot_5_Malus_TE_structure/"
dir_panTEannot_non_pop_structure = f"{work_dir}/panTEannot_5_Malus_TE_residuals/"

#from panREPET
dir_pop_structure = "panREPET_Malus_5accs_structure_95"
dir_non_pop_structure = "panREPET_Malus_5accs_residuals_95"
f_pop_structure = f"{work_dir}/{dir_pop_structure}/clique/filter80/cliques_copy_filter80.tsv"
f_non_pop_structure = f"{work_dir}/{dir_non_pop_structure}/clique/filter80/cliques_copy_filter80.tsv"

f_unique_popStruc = f"{work_dir}/{dir_pop_structure}/uniqueCopies/filter80/concat_uniqueCopies_filter80.tsv"
f_unique_nonpopStruc = f"{work_dir}/{dir_non_pop_structure}/uniqueCopies/filter80/concat_uniqueCopies_filter80.tsv"

f_TElibrary_popStruc = "../../data/TElists/250923_TE_structure.txt"
f_TElibrary_nonpopStruc = "../../data/TElists/250923_TE_residuals.txt"

# Aurélie
f_genes = "../../data/genes/concat.emapper.decorated.reformatChr.gff" # sed -E 's/Chr/chr/g; s/chr0+([0-9])/\1/g' concat.emapper.decorated.gff > concat.emapper.decorated.reformatChr.gff
#f_genes_Mdo = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/genes/Mdo.emapper.decorated.gff.gene"
f_GOs = "../../data/genes/concat.emapper.annotations.reformatChr" # sed -E 's/Chr([0-9]+)/chr\1/g; s/chr0+([0-9])([^0-9])?/chr\1\2/g' concat.emapper.annotations > concat.emapper.annotations.reformatChr

#Xilong
# f_genes_Xilong = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/gene_models_20170612.gene.sort.01.flank2k.bed"
# f_genes_Xilong_info = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/Mdomestica_491_v1.1.annotation_info.txt"
# f_genes_Xilong_floraison = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/flowering_GDDH_florid.txt"
#f_genes_GDDH13 = "../../gene_models_20170612.gff3.gene" #From Neo /data/atipe-workspace/xchen/assembly_GDDH13_v1_1
# f_positive_genes = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/positive_gene_table.bed" # on M.domestica/GDDH13_1-1 genomes
# f_balancing_genes = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/balancing_gene_table.bed" # on M.domestica/GDDH13_1-1 genomes

# color dictionnaries
d_order = {'LTR': 'blue', 'TIR': 'pink', 'noCat': 'white', 'LINE': 'yellow', 'SINE': 'orange', 'pararetrovirus': 'black', 'Helitron': 'red'}
d_order_color = {
	'LINE': 'yellow', #'#f2e914',
	'LTR_Other': 'orange', #'#faa11c',
	'LTR_Copia': 'lightgreen',
	'LTR_Gypsy': 'lightblue',
	'TIR': 'orchid', #'purple',
	'noCat': 'lightgrey',
	'Helitron': 'pink', #update
	'MITE': 'blueviolet', #update
	#'SINE': 'olive', #update
	'SSR': 'brown', #update
	'pararetrovirus': 'cyan', #update
}

d_gene_color = {
	'1. No genes at 2000 bp on either side of the TE': 'grey',
	'4. Response to stimulus': 'red',
	'3. No response to stimulus': 'green',
	'2. GO unavailable': 'lightblue'
}
d_ass_names_old = {"M.domestica": "MdoGDDH", "M.sieversii":"MsiKAZ", "M.sylvestris":"MsyGERA", "M.orientalis":"MorLC", "MdoR106A": "MdoR106A"}
d_ass_names = {"M.domestica": "Mdo", "M.sieversii":"Msie", "M.sylvestris":"Msyl", "M.orientalis":"Mori", "MdoR106A": "MdoR"}
d_species_color = {"green" : "Wild", "red": "Domesticated", "black": "Both", "brown": "Shared by all"}
d_species_color2 = {"3. GDDH13-specific": "green", "4. Rootstock-specific": "green", "3. M.sieversii":"green", "4. M.sylvestris":"green", "5. M.orientalis":"green", "1. Wild-specific": "red", "2. Domesticated-specific": "green"}
d_species_shape = {"3. GDDH13-specific": '<', "4. Rootstock-specific": '>', "3. M.sieversii": 'X', "4. M.sylvestris": 's', "5. M.orientalis": 'o',  "1. Wild-specific": "^", "2. Domesticated-specific": "v"}


# f_all_TEs_pop = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/panTEannot/workflow/results_Anthony_Malus_5ass_TE_structure/M.domestica/run_blaster-S0.align.clean_match.path.total.out.coordchr.merged.bed"
# f_all_TEs_non_pop = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/panTEannot/workflow/results_Anthony_Malus_5ass_TE_residuals_11March2025/M.domestica/run_blaster-S0.align.clean_match.path.total.out.coordchr.merged.bed"
# genome_GDDH13_pop = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/panTEannot/workflow/results_Anthony_Malus_5ass_TE_structure/M.domestica/M.domestica.fa.formated"
# genome_GDDH13_non_pop = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/panTEannot/workflow/results_Anthony_Malus_5ass_TE_residuals_11March2025/M.domestica/M.domestica.fa.formated"


# consensus with high TE copies
# f_highCopiesCount = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/highCopiesCount.tsv"

#PASTEC classif
# f_classif = "/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/Database_pourContig.xlsx"

# f_genes_GDDH13_flowering = "genes_GDDH13_flowering.tsv"


#---------#
# Outputs #
#---------#
outdir = "../../output"

if not os.path.exists("../../output/TE_structure"):
    os.mkdir("../../output/TE_structure")
if not os.path.exists("../../output/TE_residuals"):
    os.mkdir("../../output/TE_residuals")





def reformat( f, f_unique, dir_panTEannot, extension ):
    
    df = pd.read_csv(f, sep="\t")
    df = df.rename(columns={"id": "TE insertions"})
    
    out_file = f"../../output/{extension}/panREPET_identities_{extension}.tsv"
    
    #if os.path.exists(out_file):
    get_identity( df, dir_panTEannot, extension )
    
    df_identity = pd.read_csv(out_file, sep="\t")
    #print(f"df before merging with copy-consensus identity:\n{df}")
    df_identity = df_identity[["accession", "copy_name", "identity"]]
    df_identity = df_identity.drop_duplicates(subset=["accession", "copy_name"])
    df = df.merge(df_identity , on=["accession", "copy_name"], how="left")
    #print(f"df after merging with copy-consensus identity:\n{df}")
    
    df_unique = pd.read_csv( f_unique, sep="\t", names=['chromosome', 'tool', 'type', 'start', 'end', 'length', 'strand', '.', 'features', 'accession'] )
    df_unique['clique_size'] = 1
    df_unique.insert(0, 'TE insertions', range(1, 1 + len(df_unique)))
    df_unique['TE insertions'] = df_unique['TE insertions'].astype(str) + '_unique'
    df_unique['copy_name'] = df_unique['features'].apply(lambda x: x.split(';')[0].split('=')[1] + '_' + x.split(';')[1].split('=')[1])
    df_unique['consensus'] = df_unique['copy_name'].apply(lambda x: x.replace(x.split('_')[0] + '_', ''))
    df_unique['identity'] = df_unique['features'].apply(lambda x: x.split('identity:')[1]).astype('float')
    df_unique = df_unique[["TE insertions", "accession", "chromosome", "start", "end", "copy_name", "consensus", "clique_size", "identity"]]
    
    df = pd.concat([df, df_unique])
    
    #df["order"] = df["consensus"].str.split("_").str[2]
    df["order"] = df["consensus"].str.split("_").apply(lambda x: "_".join(x[2:4]) if len(x) > 3 else x[2])
    df["order"] = df["order"].replace("LTR", "LTR_Other")
    
    df["length"] = df["end"] - df["start"]

    out_file_concatUniques = f"../../output/{extension}/panREPET_identities_{extension}_concatUniques.tsv"
    df.to_csv(out_file_concatUniques, sep="\t")

    return df








def get_identity( df, dir_panTEannot, extension ):


	df_asss = pd.DataFrame()

	cols_bed = ["chr", "start", "end", "copy_name", "score", "phase", "start2", "end2", "identity", "nb_hits", "none"]

	for ass in os.listdir(dir_panTEannot):
		if ass != "log":

			print(ass)

			ass_bed = f"{dir_panTEannot}{ass}/run_blaster-S1.align.clean_match.path.total.out.coordchr.merged.bed"

			df_ass = pd.read_csv(ass_bed, sep="\t", names = cols_bed)
			df_ass = df_ass[["copy_name", "identity"]]
			df_ass["accession"] = ass
			df_ass["copy_name"] = df_ass["copy_name"].str.split(",").str[0]

			df_asss = pd.concat([df_asss, df_ass])
			#print(df_asss)


	df = df.merge(df_asss, on=["accession", "copy_name"], how="left")
	df.to_csv(f"../../output/{extension}/panREPET_identities_{extension}.tsv", sep="\t", index=False)








def clustermap( df, title, ax ):
	'''
		by default: method='average', metric='euclidean'
	'''

	size = 35
	title = title
	df = df[df['clique_size'] != 1]

	df_crosstab = pd.crosstab(df['TE insertions'], df['accession'])
	df_crosstab = df_crosstab[0:1000] # 'presence-absence' matrix
	#print("\nCrosstab:\n", df_crosstab)


	# Rows color
	index_to_color = df.set_index('TE insertions')['order'].map(d_order).to_dict()
	row_colors = [index_to_color.get(index, 'white') for index in df_crosstab.index]
	df_row_colors = pd.DataFrame({"Order": row_colors}, index=df_crosstab.index)


	cm = sns.clustermap(df_crosstab,
				figsize=(23,20),
				#metric=metric, 
				method="ward",
				xticklabels=True, 
				yticklabels=False,
				row_colors=df_row_colors,
				#dendrogram_ratio = 0.1,
				cmap=['black'],
				mask= df_crosstab == 0,
				cbar_pos=None)
	
	# text size
	cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(), fontsize=size-10, rotation=85)
	cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=size)
	cm.ax_heatmap.set_xlabel("Assembly", fontsize=size, labelpad=25)
	cm.ax_heatmap.set_ylabel(cm.ax_heatmap.get_ylabel(), fontsize=size, labelpad=25)
	plt.title(title, fontsize=50, loc='center', fontweight='bold', x=0.4, y=1.5) # x=-0.5, y=0.9) 

	## Legends

	# TE Presence-Absence
	if title == "TE non population structure":
		legend_labels = ['Presence', 'Absence']
		legend_colors = ['black', 'white']
		legend_patches = [Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]
		legend = cm.ax_heatmap.legend(handles=legend_patches, title='TE insertions', loc='lower left', bbox_to_anchor=(-0.5, -0.35), prop={'size': size+4})
		legend.get_title().set_fontsize(size+4)
		for patch in legend.get_patches():
			patch.set_edgecolor('black')

	if title == "TE population structure":
		# Order
		legend_patches3 = [Patch(color=color, label=label) for label, color in d_order.items()]
		legend3 = cm.ax_col_dendrogram.legend(handles=legend_patches3, title='Order', loc='lower left', bbox_to_anchor=(-0.6, 0.5), prop={'size': size}, ncol=2)
		legend3.get_title().set_fontsize(size)


	fig_name = f'clustermap_panREPET_Anthony_{title}.png'
	cm.savefig( fig_name, dpi=500, bbox_inches='tight', pad_inches=0.1 )
	fig = plt.gcf()
	plt.close(fig)
	img1 = plt.imread(f'{fig_name}')
	ax.imshow(img1)
	ax.axis('off')
	os.remove(fig_name)

	return df_crosstab









def TElibray_stat( f_TElibrary, ax, title ):

	print("\n\n****** TE library statistics starting... ******")

	df = pd.read_csv( f_TElibrary, names=["consensus"] )
	#print(df)
	df["order"] = df["consensus"].str.split("_").str[2]
	df = df.groupby(["order"]).count()
	#print(df)

	colors = [d_order.get(consensus, 'white') for consensus in df.index]
	df.plot.pie( y='consensus', colors=colors, labels=None, autopct=lambda p: '{:.1f}'.format(p) if p > 0.5 else '', pctdistance=1.27, labeldistance=1.05, textprops={'fontsize': 20, 'rotation': 0}, ax=ax, wedgeprops={'edgecolor': 'black'} ) #autopct=lambda p: '{:.1f}'.format(p) if p > 0.5 else '', pctdistance=1.27,
	ax.legend( labels=df.index, loc='lower left', fontsize=25, ncols=8, bbox_to_anchor=(0.05, -0.6)) #-0.3))
	ax.set_ylabel('')
	ax.set_title(title, fontsize=22, pad=30)
	if title == "TE-residuals":
		ax.get_legend().set_visible(False)

	print("****** TE library statistics end******\n\n")









def piechart_panCompVSTEfamily( df, f_TElibrary, axes, nrow, sub_title ):
	'''
		type: class or order
	'''

	#title = "TE insertions {type} proportions in terms of pangenomic compartments"

	# remove TE copie redundancy from a clique
	df = df[['TE insertions', 'consensus', 'order', 'clique_size']]
	df = df.sort_values(by=['TE insertions', 'consensus'])
	df = df.drop_duplicates()
	print(f"\n Pangenomic Compartment VS TE\n", df, '\n')

	# Plot pie charts
	#title = fig.suptitle("a", fontweight='bold', x=0.1, y=0.95) 
	#title.set_fontsize(48)

	## TE library: proportion of orderss
	TElibray_stat( f_TElibrary, axes[nrow, 0], sub_title)

	df['clique_size'] = df['clique_size'].replace(1, "1 species (non shared)")
	df['clique_size'] = df['clique_size'].replace(2, "2 species")
	df['clique_size'] = df['clique_size'].replace(3, "3 species")
	df['clique_size'] = df['clique_size'].replace(4, "4 species")
	df['clique_size'] = df['clique_size'].replace(5, "5 species (core)")

	grouped = df.groupby(['clique_size', 'order']).size().unstack(fill_value=0)
	#print("\n", grouped)

	all_axes = [axes[nrow, col] for col in range(6) if not (col == 0)] #axes.flatten() #

	for i, (size_batched, ax) in enumerate(zip(grouped.index, all_axes)):
		sizes = grouped.loc[size_batched]
		labels = sizes.index
		colors = [d_order[order] for order in labels]
		wedges, texts, autotexts = ax.pie( sizes, colors=colors, pctdistance=1.27, labeldistance=1.05, labels=None, textprops={'fontsize': 15, 'rotation': 0}, autopct=lambda p: '{:.1f}%\n({:.0f})'.format(p, p * sizes.sum() / 100) if p > 0.01 else '', wedgeprops={'edgecolor': 'black'} ) #autopct=lambda p: '{:.1f}'.format(p) if p > 0.5 else '', 
		ax.set_title(f'{size_batched}', fontsize=22, pad=30)





def determine_color(intersection,sets):

	if set(intersection) == {"Mdo", "Msie", "Msyl", "Mori", "MdoR"}: # {'MdoR106A', 'M.domestica', 'M.sieversii', 'M.orientalis', 'M.sylvestris'}:
		return 'brown'
	elif set(intersection).issubset({'Mdo','MdoR'}):
		return 'red'
	elif set(intersection).issubset({'Msie', 'Mori', 'Msyl'}):
		return 'green'
	elif set(intersection) == set(sets): 
		return 'black'
	else:
		return 'black'







def format_chromosome( chrom, accession ):

	if pd.isna(chrom):
		return "chr0"

	#rename chr
	chrom = re.sub(r'chromosome(\d+)', r'chr\1', chrom, flags=re.IGNORECASE)
	chrom = re.search(r'(chr\d+[a-zA-Z]?)', chrom, re.IGNORECASE)
	chrom = chrom.group(0) if chrom else "chr0"
	chrom = re.sub(r'chr(\d+)[a-zA-Z]?', r'chr\1', chrom, flags=re.IGNORECASE)
	chrom = re.sub(r'chr0+(\d+)', r'chr\1', chrom, flags=re.IGNORECASE)

	if accession == "M.domestica" or accession == "MdoR106A":
		chrom = re.sub(r'chr(\d+)', lambda m: f'Chr{int(m.group(1)):02}', chrom)

	return chrom








def get_bedfile( df_panREPET, outDir ):

	#df_panREPET['chromosome'] = df_panREPET.apply(lambda row: format_chromosome(row['chromosome'], row['accession']), axis=1)

	df_panREPET = df_panREPET[["chromosome", "start", "end", "TE insertions", "accession", "copy_name", "consensus", "clique_size", "order"]]
	df_panREPET.to_csv(f"{outDir}/panREPET.bed", sep="\t", index=False, header=False)






def intersect_against_genes( f_TEs, f_genes, windows ):

	# Bedtools window
	f_output_bedtools = f_TEs + ".intersectVSgenes"
	cmd = f"ml all && ml gencore/3 && ml bedtools/2.31.1 && bedtools window -b {f_genes} -a {f_TEs} -w {windows} > {f_output_bedtools}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)



def parse_TEvsgenes( f_intersection, d_ass_names ):


	#chr0    36939635        36946744        1       M.domestica     112826_TE_003334_LTR_Copia      TE_003334_LTR_Copia     3       chr0    Helixer exon    36935551        36940185        .       +       .       ID=Malus_chr0_003853.1.exon.1;Parent=Malus_chr0_003853.1;Assembly_name=MdoFUJ1  550
	#chr6    12106010        12113067        1       M.sylvestris    361993_TE_003334_LTR_Copia      TE_003334_LTR_Copia     3       chr6    Helixer three_prime_UTR 12112942        12113267        .       -       .       ID=Malus_chr6_001395.1.three_prime_UTR.1;Parent=Malus_chr6_001395.1;Assembly_name=MdoHAN
    #update
    # chr4    12219623        12221000        1       Mdo     181816_TE_103530_LTR    TE_103530_LTR   2       LTR_Other
    # chr4    11958227        11959602        1       MdoR    263215_TE_103530_LTR    TE_103530_LTR   2       LTR_Other

	header = "#chr\tstart\tend\tclique_id\taccession\tcopy_name\tconsensus\tclique_size\torder\tchr_gene\tsource_gene\ttype_gene\tstart_gene\tend_gene\tscore_gene\tstrand_gene\tphase_gene\tfeatures_gene\n"

	with open(f_intersection, "r") as f, open(f_intersection+".sameAss", "w") as f_sameAss:
		lines = f.readlines()
		f_sameAss.write(header)
		for line in lines:
			cutting = line.strip().split('\t')

			#TE
			assName_TE = cutting[4] #d_ass_names[cutting[4]] --> torm
			length_TE = abs(int(cutting[2]) - int(cutting[1]))
			chr_TE = cutting[0]

			#gene
			assName_gene = cutting[17].split(";Assembly_name=")[1]
			gene_type = cutting[11]

			#nb_bp_overlapping = int(cutting[17])
			#percent_overlapping = nb_bp_overlapping / length_TE * 100

			if assName_TE == assName_gene and chr_TE != "chr0" and chr_TE != "Chr00" and gene_type == "gene": #nb_bp_overlapping >=80 and 
				f_sameAss.write(line)

	print(f"\nFiltering TE x genes intersection {f_intersection+'.sameAss'}")






def merge_GO( f_intersection ):

    print('\n\n\n******************* Start - Merge GO:\n')
    
    df = pd.read_csv( f_intersection, sep="\t" )
    df["query"] = df["features_gene"].str.split(";Assembly_name").str[0].str.replace("ID=", "")
    df["Assembly_name"] = df["accession"] #.map(d_ass_names) #df["features_gene"].str.split("Assembly_name;").str[1]
    df['inGene'] = (df['start'] >= df['start_gene']) & (df['end'] <= df['end_gene'])
    print(df)
    
    df_GO =  pd.read_csv( f_GOs, sep="\t" )
    df_GO["query"] = df_GO["query"].str.split(".").str[0]
    df_GO = df_GO[["query", "Assembly_name", "GOs", "Description"]]
    print(df_GO)
    
    df_merged = df.merge(df_GO, on=["query", "Assembly_name"])
    #print(df_merged)
    output = f_intersection + ".GO"
    df_merged.to_csv(output, sep="\t", index=False)
    print(f"\nMerge with GO {output}")
    
    print('\n******************* End - Merge GO:\n\n\n')
    
    




def intersect_GO( df, outDir, window ):

	if not os.path.exists( outDir ):
		os.mkdir( outDir )

	# get bed file format
	get_bedfile( df, outDir )
	
	# intersect panREPET against genes
	intersect_against_genes( f"{outDir}/panREPET.bed", f_genes, window )

	# keep intersect from same accessions
	parse_TEvsgenes( f"{outDir}/panREPET.bed.intersectVSgenes", d_ass_names )

	# merge with GO
	merge_GO( f"{outDir}/panREPET.bed.intersectVSgenes.sameAss" )








def filter_clique1(group):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement 'M.sieversii', 'M.orientalis', 'M.sylvestris'
	return all(group['accession'].isin(['M.sieversii', 'M.orientalis', 'M.sylvestris']))

def filter_clique2(group):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement "M.domestica" ou "MdoR106A"
	return all(group['accession'].isin(['M.domestica', 'MdoR106A']))

def filter_clique3(group, sp):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement 'M.sieversii', 'M.orientalis', 'M.sylvestris'
	return all(group['accession'].isin([sp]))


def GOEA( f, inTE, ax, title, title2, species, order ):


	print("\n\n\n")
	print("***************************************************")
	print("************        Start GOEA        *************")
	print("***************************************************")
	print("\n")

	print("\ngoatools version : ", goatools.__version__)

	df = pd.read_csv(f, sep="\t")
	print(f"\nInitial df:\n{df}")

	if order != "all":
		df = df[df["order"]==order]
		print(f"\nSelecting TE order {order}:\n{df}")

	if inTE!= "both":
		df = df[df["inGene"] == inTE]
	print(f"\nFiltering ... params 'genes in TE' == {inTE}:\n{df}")

	if df.empty:
		return f"Empty df for {f}"

	species_name = species

	# select one species to avoid redundancy and/or select wild or cultivar species
	if species == "wild":
		print(df)
		#df_filtered = df[df['accession'].isin(['M.sieversii', 'M.orientalis', 'M.sylvestris'])]
		df = df.groupby('clique_id').filter(filter_clique1)
		print(df)
		species = "M.sieversii"
		df = df[df["accession"] == species]
		print(df)

	elif species == "cultivar":
		#df_filtered = df[df['accession'].isin(['M.domestica', 'MdoR106A'])]
		#clique_ids_concernes = df_filtered['clique_id'].unique()
		#df = df[df['clique_id'].isin(clique_ids_concernes)]
		df = df.groupby('clique_id').filter(filter_clique2)
		species = "M.domestica"
		df = df[df["accession"] == species]

	else:
		#df = df[df["accession"] == species] #to rm because I want unique TE insertions
		#df = df.groupby('clique_id').filter(filter_clique3, species)
		df = df.groupby('clique_id').filter(lambda group: filter_clique3(group, species))

	df = df[df['GOs'] != '-']
	#print(df)
	df = df[["clique_id", "accession", "query", "GOs"]]
	#print(df)
	df.to_csv(f + f".{inTE}.{species_name}.{order}", sep="\t", index=False)

	### Stats
	# Nombre de clique_id uniques
	nb_cliques = df[["clique_id", "GOs"]].nunique()
	# Calculer la médiane et la moyenne des nombres de GO
	df["GO_count"] = df["GOs"].apply(lambda x: len(x.split(",")))
	median_go = df["GO_count"].median()
	average_go = df["GO_count"].mean()
	# Affichage des résultats
	print(f"{nb_cliques} cliques having GOs")
	print(f"{median_go} GOs per gene in median")
	print(f"{average_go} GOs per gene in average")


	# ========> 1. Charger le fichier OBO GO
	go_obo_file = "go-basic.obo"
	go_dag = GODag(go_obo_file)
	go = obo_parser.GODag(go_obo_file)

	# ========> 2. Population de gènes (ensemble de gènes d'intérêt)
	df_annotation = pd.read_csv(f_GOs, sep="\t") #la population c'est ton annotation en gène de départ # example {'GENE1', 'GENE2', 'GENE3', 'GENE4'}
	print(f"\ndf_genes:\n{df_annotation}\n")
	df_annotation = df_annotation[df_annotation["Assembly_name"] == d_ass_names[species]]
	print(f"\ndf_genes ({species}):\n{df_annotation}\n")
	print(f"\nnb of query:\n{len(df_annotation['query'].unique())}\n")
	df_annotation["query"] = df_annotation["query"].str.split(".").str[0]
	df_valid = df_annotation[df_annotation["GOs"].notna() & (df_annotation["GOs"] != "-")]
	print(f"\nnb of query having GO:\n{len(df_valid['query'].unique())}\n")
	print(f"\nreformated df_genes:\n{df_valid}\n")
	pop = set(df_valid["query"])
	#print(pop)
	print(f"\nPopulation de référence (pop) : {len(pop)} gènes")


	# ========> 3. Associations GO pour les gènes
	#assoc = {
	#	'GENE1': {'GO:0008150', 'GO:0007049'},
	#	'GENE2': {'GO:0007049', 'GO:0051301'},
	#}
	assoc = {}
	for index, row in df.iterrows():
		gene = row["query"] #['clique_id']
		if row['GOs'] != '-':
			go_terms = row['GOs'].split(',')
			assoc[gene] = set(go_terms)
		else:
			assoc[gene] = set()
	#print(assoc)

	# ========> 3.1. Filtrer les termes GO de niveau x
	filtered_go_terms = {term: go_dag[term] for term in go_dag if go_dag[term].level == 1}

	# ========> 3.2. Filtrer les associations GO pour ne garder que les termes GO de niveau x
	filtered_assoc = {}
	for gene, terms in assoc.items():
		filtered_assoc[gene] = {term for term in terms if term in filtered_go_terms}


	if not filtered_assoc:
		print(f"Error: Associations file is empty or does not exist.")
	else:
		# ========> 4. Créer l'objet GOEnrichmentStudy
		goeaobj = GOEnrichmentStudy(
			pop=pop,
			assoc=filtered_assoc, #assoc,
			obo_dag=go_dag,
			propagate_counts=True,
			alpha=0.05,
			methods=['fdr', 'bonferroni']
		)


		# ========> 5. Lancer l'enrichissement GO
		study = df["query"].to_list()
		goea_results = goeaobj.run_study(study)
		#for result in goea_results:
		#	print(result)
		#	print(record.__dict__) 

		if not goea_results:
			print(f"No GO terms found for {f} {inTE} {species} {order}.")
		else:

			output = f+f'.goea_results_{inTE}.{species_name}.{order}.tsv'
			with open( output, 'w' ) as fout_tsv:
				goeaobj.prt_tsv(fout_tsv, goea_results)

			results_df = pd.read_csv( output, sep="\t" )
			
			results_df = results_df[results_df["NS"] == "BP"]
			#print(results_df)

			# Trier par p-value ajustée (fdr_bh)
			results_df = results_df.sort_values(by=['p_fdr', 'study_count'], ascending=[True, False])
			#print(results_df)
			# Filter par p-value ajustée (fdr_bh)
			#results_df = results_df[results_df['p_fdr'] < 0.001]
			results_df = results_df.head(35)
			print(f"\nGOEA df result:\n{results_df}\n")

			# Visualiser les résultats avec un dotplot
			sns.scatterplot(
				ax=ax,
				data=results_df,
				y='name',  # Le nom du terme GO
				x='p_fdr',  # P-value ajustée p_bonferroni
				size='study_count',  # Taille du terme GO
				hue='study_count',  # Nombre d'occurrences dans les gènes d'intérêt
				palette='coolwarm',  # Palette de couleurs
				#sizes=(50, 500)  # Plage des tailles des points
			)

			if ax != None:
				ax.set_title(title, fontweight='bold', loc='left', fontsize=15)
				ax.set_title(title2, fontsize=11) #'GO Enrichment Study'
				ax.set_ylabel('Significant GO Terms') #('Top 30 most represented GO Terms description')
				#ax.axvline(x=0.05, color='red', linestyle='--', linewidth=1)
				#ax.axvline(x=0.01, color='red', linestyle='--', linewidth=1)
				plt.tight_layout()
				#plt.show()

				labels_to_color = ["response to stimulus", "cellular response to stimulus", "response to stress", "response to hormone", "response to external stimulus"]
				for label in ax.get_yticklabels():
					text = label.get_text()
					label.set_size(9)
					#if text in labels_to_color:
					if "response to" in text or "stimulus" in text:
						label.set_color("red")
					


	print("\n")
	print("***************************************************")
	print("************        End GOEA          *************")
	print("***************************************************")
	print("\n\n\n")




# Function to get level 1 GO term names
def get_level_1_go_names( go_terms, go_dag ):

	if go_terms == "1. No genes at 2000 bp on either side of the TE":
		return go_terms
	else:
		go_terms_list = go_terms.split(',')
		#print(go_terms_list)
		go_terms_level0_list = []
		for go_term in go_terms_list:
			try:
				term_info = go_dag[go_term]
				if term_info.level == 1: #'biological_process'
					go_terms_level0_list.append(term_info.name)
			except KeyError:
				""
        			#print(f"Clé {go_term} non trouvée dans go_dag")

		#print(go_terms_level0_list)

	if len(go_terms_level0_list) == 0:
		return "2. GO unavailable"
	if "response to stimulus" in go_terms_level0_list:
		return "4. Response to stimulus"
	else:
		return "3. No response to stimulus" #str(set(go_terms_level0_list))




def add_genes_presence_absence( df, outDir ):
	"""
		Add for each TE, a new column "gene"
		-> True: at least one gene around the TE
		-> False: no gene
	"""

	df["TE insertions"] = df["TE insertions"].astype(str)

	f_TEs_vs_genes = f"{outDir}/panREPET.bed.intersectVSgenes.sameAss.GO"
	#chr    start   end     clique_id       accession       copy_name       consensus       clique_size     order   chr_gene        source_gene     type_gene       start_gene      end_gene        score_gene      strand_gene     phase_gene      features_gene   query   Assembly_name   inGene  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description     Preferred_name  GOs     EC      KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE    KEGG_TC        CAZy    BiGG_Reaction   PFAMs
	#Chr17   25145519        25155060        16      MdoR106A        103582_TE_104914_LTR    TE_104914_LTR   4       LTR     Chr17   Helixer gene    25142109        25145261        .       -       .       ID=MdoR_Chr17_001578;Assembly_name=MdoR MdoR_Chr17_001578       MdoR    False   225117.XP_009336308.1   0.0     1324.0  28JRK@1|root,2QS51@2759|Eukaryota,37MB5@33090|Viridiplantae,3G9W3@35493|Streptophyta35493|Streptophyta      S       Protein FAR1-RELATED SEQUENCE   -       -       -       ko:K17604       -       -       -       -       ko00000,ko01009 -       -       -       FAR1,MULE,SWIM

	# Fetch TE that have genes aroung them
	df_TEs_with_genes = pd.read_csv(f_TEs_vs_genes, sep="\t")
	#print(df_TEs_with_genes)

	"""
	clique_ids_with_genes = df_TEs_with_genes["clique_id"].unique().tolist()

	df["genes"] = df["TE insertions"].isin(clique_ids_with_genes)
	df["Description"] = df["TE insertions"].isin(clique_ids_with_genes)
	"""

	# Create a mapping from clique_id to Description
	#clique_to_description = df_TEs_with_genes.set_index("clique_id")["GOs"].to_dict()
	#print(df_TEs_with_genes.set_index("clique_id")["GOs"])
	clique_to_description = df_TEs_with_genes.groupby('clique_id')['GOs'].apply(lambda x: ','.join(x)).to_dict()
	#print(clique_to_description)

	# Get the list of clique_ids that have genes around them
	clique_ids_with_genes = df_TEs_with_genes["clique_id"].unique().tolist()

	# Add new columns to the original DataFrame
	df["genes"] = df["TE insertions"].isin(clique_ids_with_genes)
	df["GOs"] = df["TE insertions"].map(clique_to_description).fillna("1. No genes at 2000 bp on either side of the TE")
	#df["GOs"] = df["TE insertions"].map(clique_to_description).notna().map({True: "ok", False: "NA"})

	# Load the OBO GO database: wget https://current.geneontology.org/ontology/go-basic.obo
	go_obo_file = "go-basic.obo"
	goDag = GODag(go_obo_file)

	df['level_1_go_names'] = df['GOs'].apply( get_level_1_go_names, go_dag=goDag )


	return df





def check_both_annotation( ):
	"""
		Compare Aurelie's annotation VS Daccord et al's annotation (GDDH13)
	"""

	# Bedtools intersect Aurelie's annotation VS Daccord et al's annotation (GDDH13)
	f_output_bedtools = "compare_Aurelie_vs_Daccord_gene_annotation.tsv"
	cmd = f"bedtools intersect -wo -a {f_genes_Mdo} -b {f_genes_GDDH13} > {f_output_bedtools}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)





def add_gene_selection( df, outDir ):



	f_TEs = f"{outDir}/panREPET.bed.intersectVSgenes.sameAss"
	df_TEs = pd.read_csv(f_TEs, sep="\t")
	#select GDHH13
	df_TEs = df_TEs[df_TEs["accession"] == "M.domestica"]
	# reformat chromosomes
	df_TEs = df_TEs[df_TEs["accession"] == "M.domestica"]
	#df_TEs['chr'] = df_TEs.apply(lambda row: format_chromosome(row['chr'], row['accession']), axis=1) # to rm ?
	f_TEs_GDDH13 = f_TEs + ".Mdomestica"
	df_TEs.to_csv(f_TEs_GDDH13, sep="\t", header=False, index=False)

	# Bedtools windows POSITIVE
	f_output_bedtools = f_TEs_GDDH13 + ".intersectVSpositiveselection"
	cmd = f"bedtools window -b {f_positive_genes} -a {f_TEs_GDDH13} -w 10000 > {f_output_bedtools}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)


	# Bedtools intersect POSITIVE
	f_output_bedtools = f_TEs_GDDH13 + ".intersectVSbalancingselection"
	cmd = f"bedtools window -b {f_balancing_genes} -a {f_TEs_GDDH13} -w 10000 > {f_output_bedtools}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)




def upset_plot( df, ax, title ):


    print("\n\n****** Start - Upset plot ... ******")
    
    #sets = ['MdoR106A', 'M.domestica', 'M.sieversii', 'M.orientalis', 'M.sylvestris']
    sets = ["Mdo", "Msie", "Msyl", "Mori", "MdoR"]
    df = df[["accession", "TE insertions", "identity", "genes", "order", "GOs", "level_1_go_names", "length"]]
    print(df)
    print(df[["accession", "TE insertions", "genes"]])
    
    #df = df[df['level_1_go_names'] != "1. No genes at 2000 bp on either side of the TE"]
    #df = df[df['order'] == "LTR_Gypsy"]
    
    # Create binary table
    df_bool = df[["accession", "TE insertions"]].pivot_table(index='TE insertions', columns='accession', aggfunc=lambda x: True).fillna(False)
    print(df_bool)
    
    # merge with initial df for fetching other informations like genes presence
    df = pd.merge(df_bool, df, on='TE insertions', how='left')
    print(df)
    
    df = df.drop('accession', axis=1)
    df = df.drop_duplicates()
    print(df)

    #cols = ['TE insertions', 'M.domestica', 'M.orientalis', 'M.sieversii', 'M.sylvestris', 'MdoR106A']
    cols = ['TE insertions', 'Mdo', 'Mori', 'Msie', 'Msyl', 'MdoR']
    #df = df.groupby(cols).agg(identity=('identity', 'median'), length=('length', 'median'), genes=('genes', lambda x: 1 if x.any() else 0), order=('order', 'first'), level_1_go_names=('level_1_go_names', 'first')).reset_index()
    df = df.groupby(cols).agg(
        order=('order', 'first'),
        identity=('identity', 'median'),
        length=('length', 'median'),
        genes=('genes', lambda x: 1 if x.any() else 0),
        level_1_go_names=('level_1_go_names', lambda x: x.dropna().iloc[0] if x.dropna().size > 0 else np.nan)
    ).reset_index()
    print(df)

    df.set_index(sets, inplace=True)
    print(df)
    df.sort_index(inplace=True)
    print(df)



    # ===========================================================
    # 🔹 Création du DataFrame df_orderPerc : proportion d'ordre par intersection
    # ===========================================================

    # On remet à plat le DataFrame avec les colonnes booléennes
    df_reset = df.reset_index()

    # Compter le nombre unique de TE insertions par intersection et par ordre
    count_order = (
        df_reset.groupby(sets + ['order'])
                .agg(count=('TE insertions', 'nunique'))
                .reset_index()
    )
    count_GO = (
        df_reset.groupby(sets + ['level_1_go_names'])
                .agg(count=('TE insertions', 'nunique'))
                .reset_index()
    )
    
    # Calcul du total par intersection
    count_order['total'] = count_order.groupby(sets)['count'].transform('sum')
    count_GO['total'] = count_GO.groupby(sets)['count'].transform('sum')

    # Calcul du pourcentage (0–100)
    count_order['order_percent'] = 100 * count_order['count'] / count_order['total']
    count_GO['GO_percent'] = 100 * count_GO['count'] / count_GO['total']

    # Renommage pour plus de clarté
    df_orderPerc = count_order[sets + ['order', 'order_percent']].copy()
    df_GOPerc = count_GO[sets + ['level_1_go_names', 'GO_percent']].copy()

    # Optionnel : trier par intersection et ordre
    df_orderPerc = df_orderPerc.sort_values(by=['order_percent', 'order']).reset_index(drop=True)
    df_GOPerc = df_GOPerc.sort_values(by=['level_1_go_names']).reset_index(drop=True)

    print("\n--- Aperçu de df_orderPerc ---")
    print(df_orderPerc.head(15))
    df_orderPerc.to_csv(f"order_percent_{title}.csv", sep=',')

    # Fusion des pourcentages dans df
    df_reset = df.reset_index()
    df_reset = df_reset.merge(
        df_orderPerc,
        on=["Mdo", "Msie", "Msyl", "Mori", "MdoR", "order"],
        how="left"
    )
    df_reset = df_reset.merge(
        df_GOPerc,
        on=["Mdo", "Msie", "Msyl", "Mori", "MdoR", "level_1_go_names"],
        how="left"
    )
    df = df_reset.set_index(["Mdo", "Msie", "Msyl", "Mori", "MdoR"]).sort_index()
    print(df)
    #print(df.head(15))
    df = df.sort_values(by=['level_1_go_names', 'order_percent', 'order'])
    print(df)
  
  
    upset = UpSet(df, sort_by='cardinality', show_counts=True, include_empty_subsets=True, show_percentages=False) #show_counts='%0.2f') #show_counts='%d')
    #upset.add_stacked_bars(by="level_1_go_names", title="GO level 1", elements=10, colors=d_gene_color, sum_over="genes")
    upset.add_catplot(
        kind="bar",
        value="GO_percent",
        hue="level_1_go_names",           # pour colorer selon l’ordre TE
        dodge=True,            # sépare les barres par couleur
        palette=d_gene_color, # ta palette de couleurs
        elements=8,            # taille visuelle
        width=0.8,        # <— plus large (valeurs entre 0 et 1, par défaut 0.8)
        linewidth=0.1,    # <— contour des barres plus visible
        edgecolor="black" # <— optionnel : ajoute un contour noir
    )
    
    #Proportion of shared TE insertions with adjacent genes (within 2 kb), an insertion considered gene-associated if one copy is near a gene in any assembly.
    upset.add_catplot(value="genes", kind="bar", color="lightblue")
    
    upset.add_catplot(value="identity", kind="violin", color="lightblue", elements=6, scale="width") #, jitter=True, dodge=True)
    upset.add_catplot(value="length", kind="violin", color="lightblue", elements=6, scale="width", cut=0) #, jitter=True, dodge=True) inner="box", 
    #upset.add_stacked_bars(by='order', title="Intersection size and TE order", elements=10, colors=d_order_color, linewidth=0.3)#, sum_over='prop')
    upset.add_catplot(
        kind="bar",
        value="order_percent",
        hue="order",           # pour colorer selon l’ordre TE
        dodge=True,            # sépare les barres par couleur
        palette=d_order_color, # ta palette de couleurs
        elements=8,            # taille visuelle
        width=1,        # <— plus large (valeurs entre 0 et 1, par défaut 0.8)
        linewidth=0.3,    # <— contour des barres plus visible
        edgecolor="black" # <— optionnel : ajoute un contour noir
    )

    for i in range(1, len(sets) + 1):  # Boucle sur la taille des combinaisons
        for comb in combinations(sets, i):
            color = determine_color(comb,sets)  # Déterminer la couleur pour l'intersection
            upset.style_subsets(
                present=list(comb),  # Intersection spécifique
                facecolor=color,
                edgecolor=color,
                linewidth=2,
                label=d_species_color[color]
                #label=",".join(comb)  # Label indiquant les ensembles dans l'intersection
    )
    
    axes = upset.plot()
    # Exemple : récupérer le subplot principal (barres d’intersection)
    print("\n--- Clés disponibles dans axes ---")
    print(axes.keys())#dict_keys(['matrix', 'shading', 'totals', 'intersections', 'extra1', 'extra2', 'extra3', 'extra4', 'extra5'])
    
    # Déplacer la légende associée
    leg = axes['extra1'].get_legend()
    if leg:
        leg.set_bbox_to_anchor((1.05, 1))
        leg.set_loc('upper left')
    leg = axes['extra5'].get_legend()
    if leg:
        leg.set_bbox_to_anchor((1.05, 1))
        leg.set_loc('upper left')
        
    fig_name = f"../../Figures/upset_plot_Anthony_{title}.png"
    plt.savefig( fig_name, dpi=600, bbox_inches='tight', pad_inches=0.1 )
    plt.savefig( fig_name + ".svg", dpi=600, bbox_inches='tight', format="svg" )
    plt.savefig( fig_name + ".pdf", dpi=600, bbox_inches='tight', format="pdf" )
    fig = plt.gcf()
    plt.close(fig)
    img1 = plt.imread(f'{fig_name}')
    ax.imshow(img1)
    ax.axis('off')
    os.remove(fig_name)
    
    print("\n\n****** End - Upset plot ... ******")
    


def check_Xilong_genes( TEdir ):
	"""
		Xilong's genes: genes under a stress
	"""

	#merge Xilong genes and their annotation
	

	# Bedtools window 1
	Input = f"{TEdir}/panREPET.bed"
	Output = f"{TEdir}/panREPET.bed.Xilong"
	cmd = f"bedtools window -a {Input} -b {f_genes_Xilong} -w 2000 > {Output}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)

	# Bedtools window 2
	Input = f"{TEdir}/panREPET.bed.intersectVSgenes.sameAss"
	Output = f"{TEdir}/panREPET.bed.intersectVSgenes.sameAss.Xilong"
	#cmd = f"bedtools intersect -wo -a {Input} -b {f_genes_Xilong} > {Output}" # -wo  -nonamecheck
	cmd = f"bedtools window -a {Input} -b {f_genes_Xilong} -w 2000 > {Output}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)

	
	df_floraison = pd.read_csv(f_genes_Xilong_floraison, sep="\t")
	print(df_floraison)

	df_intersect = pd.read_csv(f"{TEdir}/panREPET.bed.Xilong", sep="\t", names=["chr", "start", "end", "clique_size", "acc", "copy_name", "consensus", "na", "order", "gene_chr", "start_gene", "end_gene", "GD_gene"])
	df_intersect = df_intersect[df_intersect["acc"] == "M.domestica"]
	print(df_intersect)


	df_merge = df_intersect.merge(df_floraison , on=["GD_gene"], how='inner')
	print(df_merge)
	count_intersection = df_merge.shape[0]
	print(f"Nombre d'intersections : {count_intersection}")


	df_merge2 = df_floraison.merge(df_intersect, on=["GD_gene"], how='inner')
	print(df_merge2)
	count_intersection = df_merge2.shape[0]
	print(f"Nombre d'intersections : {count_intersection}")





def GOEA_plot( ax, num, TE_type, order ): #"TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO.goea_results_both.M.domestica.tsv" )


	df_allSp = pd.DataFrame()

	species_list = ["wild", "cultivar", "M.domestica", "MdoR106A"] #[, 'M.sieversii', 'M.sylvestris', 'M.orientalis']:
	d_species_list = {"wild": "1. Wild-specific", "cultivar": "2. Domesticated-specific", "M.domestica": "3. GDDH13-specific", "MdoR106A": "4. Rootstock-specific"}

	for sp in species_list:

		f = f"{TE_type}/panREPET.bed.intersectVSgenes.sameAss.GO.goea_results_both.{sp}.{order}.tsv"
	
		if os.path.exists(f) and os.path.getsize(f) != 0:
			df = pd.read_csv(f, sep="\t")
			print(f"\nGOEA df {sp}:\n{df}\n")

			df = df[df["NS"] == "BP"]

			df = df[df['p_fdr'] < 0.05]

			df = df[["name", "p_fdr", 'study_count']]

			df["species"] = d_species_list[sp]

			# Merge data for comparison between species
			df_allSp = pd.concat([df_allSp, df])

		else:
			(f"{f} does not exist.")

	if not df_allSp.empty:

		print("2", f, df_allSp)
		df_allSp = df_allSp.sort_values(by=['species', 'p_fdr', 'study_count'], ascending=[True, True, False])
		print("3", f, df_allSp)

		df_allSp["name"] = df_allSp["name"].replace("biological process involved in interspecies interaction between organisms", "biological process involved in\ninterspecies interaction between organisms")
		df_allSp["name"] = df_allSp["name"].replace("biological process involved in intraspecies interaction between organisms", "biological process involved in\nintraspecies interaction between organisms")

		# Plot the differences between p_bonferroni values for each species
		sns.scatterplot(ax=ax,
				y='name', 
				x=f'p_fdr', 
				hue="species", 
				palette=d_species_color2, 
				data=df_allSp, 
				style="species", 
				markers=d_species_shape,
				alpha=0.6,
				s=200) #, s=100

		# Get the current y-axis tick labels
		y_labels = ax.get_yticklabels(True, "both")
		print(y_labels)
		# Loop through y-axis tick labels and change "response to" to red
		for label in y_labels:
			print(label.get_text())
			if 'response to' in label.get_text():
				label.set_color('red')

		# Customize the plot
		#ax.set_xticklabels(ax.get_xticks(), rotation=45)
		#ax.set_xlim(df["p_fdr"].min(), df["p_fdr"].max())
		ax.set_ylabel('Significant GO (p_fdr < 0.05)')
		ax.tick_params(axis='y', labelsize=10)
		ax.tick_params(axis='x', labelsize=9)
		if order == "LTR":
			order = "LTR_Other"
		title = f"{num}. {TE_type} ; {order}"
		ax.set_title( title ) #'GOEA : Significant GO and their p-value across Species')
		ax.legend(title='Shared TE insertions type' )
		ax.figure.tight_layout()
		#plt.savefig( "GOEA.png", dpi=500, bbox_inches='tight') #, pad_inches=0.1 )
		



def add_flowering( ):

	#genes : add floraison info
	#Chr15   IRHS-2017       gene    5301    6082    .       +       .       ID=gene:MD15G1000100;Name=MD15G1000100
	df_genes = pd.read_csv(f_genes_GDDH13, sep="\t", names=["chr", "gene_name", "start", "end", "score", "strand", "phase", "features"])
	df_genes["GD_gene"] = df_genes["features"].str.split("Name=").str[1]
	df_genes["GD_gene"] = df_genes["GD_gene"].str.split(";Note=").str[0]
	df_floraison = pd.read_csv(f_genes_Xilong_floraison, sep="\t")
	df_merge = df_genes.merge(df_floraison , on=["GD_gene"], how='inner')
	df_merge.to_csv(f_genes_GDDH13_flowering, sep="\t", index=False, header=False)



def add_GOs( ):
	"""
		genes annotation on GDDH13 by Xilong
		1. merge 2 files: genes position and the associated GOs
	"""

	df_genes = pd.read_csv(f_genes_GDDH13, sep="\t", names=["chr", "source", "gene_type", "start", "end", "score", "strand", "phase", "features"])
	df_genes["locusName"] = df_genes["features"].str.split("Name=").str[1]
	df_genes["locusName"] = df_genes["locusName"].str.split(";Note=").str[0]
	df_GOs = pd.read_csv(f_genes_Xilong_info, sep="\t")
	df_merge = df_genes.merge(df_GOs , on=["locusName"], how='left')
	df_merge = df_merge.fillna("None")
	df_merge.to_csv(f_genes_GDDH13 + ".GOs", sep="\t", index=False, header=False)




def get_fasta_for_AME( genome_GDDH13, outdir ) :

	f_TEs = f"{outdir}/panREPET.tsv"

	output_sequences = f_TEs.replace("tsv", "fasta")
	cmd = f"sort -u {f_TEs} > {f_TEs}.sortu"
	print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = f"bedtools getfasta -fi {genome_GDDH13} -bed {f_TEs}.sortu -fo {output_sequences} -name"
	print(cmd)
	subprocess.call(cmd, shell=True)




def calculate_gc_content(sequence):
	"""Calculate the GC content percentage of a sequence."""
	gc_count = sum(1 for base in sequence if base in 'GCgc')
	return (gc_count / len(sequence)) * 100

def is_valid_sequence(sequence, max_n_fraction=0.1):
	"""
	Check if a sequence is valid, i.e., it doesn't contain too many 'N's.
	:param sequence: The sequence to check.
	:param max_n_fraction: The maximum fraction of 'N's allowed in the sequence.
	:return: True if the sequence is valid, False otherwise.
	"""
	n_count = sequence.count('N') + sequence.count('n')
	n_fraction = n_count / len(sequence)
	return n_fraction <= max_n_fraction  # Allow up to max_n_fraction of 'N's

def load_gene_annotations(gff_file):
	"""Parse a GFF file and return a dictionary of gene coordinates indexed by chromosome."""
	db = gffutils.create_db(gff_file, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
	gene_annotations = {}

	# Extract gene annotations per chromosome
	for feature in db.features_of_type("gene"):
		chromosome = feature.chrom
		start = feature.start
		end = feature.end
		if chromosome not in gene_annotations:
			gene_annotations[chromosome] = []
		gene_annotations[chromosome].append((start, end))
	    
		return gene_annotations

def check_overlap(seq_start, seq_end, gene_annotations):
	"""Check if the sequence overlaps with any gene."""
	for gene_start, gene_end in gene_annotations:
		if seq_end >= gene_start and seq_start <= gene_end:  # Check for overlap
			return True
	return False


def extract_random_sequences(
	TE_origin,
	genome_origin,
	output_file,
	gff_file=None,
	gc_tolerance=10,
	max_n_fraction=0.1,
	max_attempts=5000
	):
	"""
	Extracts one random genome sequence per TE from anywhere in the genome,
	matching length and GC content, avoiding gene regions and duplicates.

	Parameters:
	TE_origin (str): FASTA with multiple TE sequences.
	genome_origin (str): Genome FASTA (multi-chromosome).
	output_file (str): Output FASTA path.
	gff_file (str): Optional GFF file for gene annotation.
	gc_tolerance (float): Max allowed GC% difference.
	max_n_fraction (float): Max N fraction allowed.
	max_attempts (int): Attempts per TE.
	"""
	# Load TE sequences
	te_records = list(SeqIO.parse(TE_origin, "fasta"))
	te_lengths = [len(r.seq) for r in te_records]
	te_gcs = [calculate_gc_content(r.seq) for r in te_records]

	# Load genome
	genome_records = list(SeqIO.parse(genome_origin, "fasta"))

	# Load GFF if provided
	gene_annotations = {}
	if gff_file:
		gene_annotations = load_gene_annotations(gff_file)

	extracted = []
	used_seqs = set()
	used_positions = set()

	for idx, (length, target_gc) in enumerate(zip(te_lengths, te_gcs)):
		found = False

		for attempt in range(max_attempts):
			chrom = random.choice(genome_records)
			chrom_seq = chrom.seq
			chrom_name = chrom.id

			if len(chrom_seq) < length:
				continue

			start = random.randint(0, len(chrom_seq) - length)
			end = start + length
			candidate = chrom_seq[start:end]
			candidate_str = str(candidate)

			if (chrom_name, start, end) in used_positions:
				continue

			if candidate_str in used_seqs:
				continue

			if not is_valid_sequence(candidate, max_n_fraction):
				continue

			if gff_file and chrom_name in gene_annotations:
				if check_overlap(start + 1, end, gene_annotations[chrom_name]):
					continue

			gc = calculate_gc_content(candidate)
			if abs(gc - target_gc) > gc_tolerance:
				continue

			# Match found
			seq_id = f"extracted_{length}_{target_gc:.2f}_{abs(gc - target_gc):.3f}_{chrom_name}:{start+1}-{end}"
			new_record = SeqRecord(candidate, id=seq_id, description="") #f"{chrom_name}:{start+1}-{end}")
			extracted.append(new_record)
			used_seqs.add(candidate_str)
			used_positions.add((chrom_name, start, end))
			found = True
			break

		if not found:
			raise ValueError(f"No match found for TE #{idx+1} with length {length} and GC {target_gc:.2f}")

	SeqIO.write(extracted, output_file, "fasta")
	print(f"✅ {len(extracted)} sequences written to {output_file}")












def get_level2_terms(go_list_str, go_dag):

	if pd.isna(go_list_str) or go_list_str=="None":  # Gère les NaN ou None
		return "None"

	go_ids = [go.strip() for go in go_list_str.split(',')]
	level2_terms = set()

	for go_id in go_ids:
		try:
			term = go_dag[go_id]
			ancestors = term.get_all_parents()
			for ancestor_id in ancestors:
				ancestor_term = go_dag[ancestor_id]
				if ancestor_term.depth == 2:
					level2_terms.add(ancestor_term.name)
		except KeyError:
			continue

	return ",".join(sorted(level2_terms)) if level2_terms else "None"






def add_GOs_from_bedtools( outdir, f ):


	#allTEs
	#Chr11   9555    10174   977626_TE_057616_TIR    565.0   +       2       653     91.159  1               TE_057616_TIR   Chr11   IRHS-2017       gene    1       35794   .       -       .       ID=gene:MD11G1000100;Name=MD11G1000100;Note=DNA binding
	#panREPET
	#Chr03   30130230        30130953        Chr03   IRHS-2017       gene    30116176        30121173        .       -       .       ID=gene:MD03G1217000;Name=MD03G1217000;Note=Met-10+ like family protein

	l_cols_allTEs = ["query_name", "query_start", "query_end", "subject_name", "score", "strand", "subject_start", "subject_end", "identity", "nb_of_subjects", "na", "consensus"]
	l_cols_panREPET = ["chr_TE", "start_TE", "end_TE"]
	l_cols_genes = ["gene_chr", "gene_source", "gene_type", "gene_start", "gene_end", "gene_score", "gene_strand", "gene_phase", "gene_features"]
	#l_cols_GOs = ["pacId", "locusName", "transcriptName", "peptideName", "Pfam", "Panther", "KOG", "ec", "KO", "GOs", "Best-hit-arabi-name", "arabi-symbol", "arabi-defline"]

	if "allTEs" in f :
		df = pd.read_csv(f, sep="\t", names=l_cols_allTEs+l_cols_genes)
		print(df)
	elif "panREPET" in f:
		df = pd.read_csv(f, sep="\t", names=l_cols_panREPET+l_cols_genes)
		#print(df)

	df["locusName"] = df["gene_features"].str.split("Name=").str[1]
	#df["locusName"] = df["locusName"].str.split(";Note=").str[0]
	df["locusName"] = df["locusName"].str.split(";").str[0]
	df_GOs = pd.read_csv(f_genes_Xilong_info, sep="\t")
	df = df.merge(df_GOs , on=["locusName"], how='left')
	df = df.fillna("None")
	print(df["locusName"])

	# Charger le fichier OBO GO
	go_obo_file = "go-basic.obo"
	go_dag = GODag(go_obo_file)
	df['GO_level2'] = df['GO'].apply(lambda gos: get_level2_terms(gos, go_dag))
	df.to_csv(f + ".GOs", sep="\t", index=False)






def bedtools_windows( f1, f2, f_output, flank_size ):

	cmd = f"bedtools window -b {f1} -a {f2} -w {flank_size} > {f_output}" # -wo  -nonamecheck
	print(cmd)
	subprocess.call(cmd, shell=True)





def launch_bedtools_window( outdir, flank_size ):

	f_GDDH13_allTE = f"{outdir}/GDDH13_allTEs.tsv"
	f_TEs = f"{outdir}/panREPET.tsv"

	#-------------------------#
	#         All TEs         #
	#-------------------------#
	# GDDH13 genes vs all TEs
	bedtools_windows( f_genes_GDDH13, f_GDDH13_allTE, f"{outdir}/GDDH13_allTEs_intersectVSgenes.tsv", flank_size )

	# GDDH13 flowering genes vs all TEs
	#bedtools_windows( f_genes_GDDH13_flowering, f_GDDH13_allTE, f"{outdir}/GDDH13_allTEs_intersectVSgenes_flowering.tsv", flank_size )


	#-------------------------#
	#     TEs on panREPET     #
	#-------------------------#
	# GDDH13 genes vs panREPET TEs
	bedtools_windows( f_genes_GDDH13, f_TEs, f"{outdir}/panREPET_intersectVSgenes.tsv", flank_size )

	# GDDH13 flowering genes vs panREPET TEs
	#bedtools_windows( f_genes_GDDH13_flowering, f_TEs, f"{outdir}/panREPET_intersectVSgenes_flowering.tsv", flank_size )

	# GDDH13 positive genes vs panREPET TEs
	#bedtools_windows( f_positive_genes, f_TEs, f"{outdir}/panREPET_intersectVSpositivegenes.tsv", flank_size )

	# GDDH13 balancing genes vs panREPET TEs
	#bedtools_windows( f_balancing_genes, f_TEs, f"{outdir}/panREPET_intersectVSbalancinggenes.tsv", flank_size )






def filter_TEs( df, outdir, f_all_TEs, TE_names_values ):

	#-----------------------#
	#  All TEs from GDDH13  #
	#-----------------------#
	#Chr11   9555    10174   132369_TE_144861_LINE   564.0   +       3       646     90.895  1       
	df_allTEs = pd.read_csv(f_all_TEs, sep="\t", names=["#chr", "start", "end", "copy_name", "score", "strand", "start_align", "size", "identity", "nb_hits", "na"])
	df_allTEs["consensus"] = df_allTEs["copy_name"].astype(str).str.split(",").str[0]
	df_allTEs["consensus"] = df_allTEs["consensus"].apply( lambda x: "_".join(x.split("_")[1:]) )
	f_GDDH13_allTE = f"{outdir}/GDDH13_allTEs.tsv"
	df_allTEs = df_allTEs[df_allTEs['consensus'].isin(TE_names_values)]
	df_allTEs.to_csv(f_GDDH13_allTE, sep="\t", index=False)#, header=False)

	#----------------------------#
	#  panREPET TEs from GDDH13  #
	#----------------------------#
	f_TEs = f"{outdir}/panREPET.tsv"
	f_TEs2 = f"{outdir}/panREPET.tsv2"
	print(TE_names_values)
	df= df[df['consensus'].isin(TE_names_values)]
	df = df[["chromosome", "start", "end"]] #"consensus", "copy_name", "TE insertions"]]
	df = df.drop_duplicates()
	df.to_csv(f_TEs, sep="\t", index=False, header=False)

	df["length"] = df["end"] - df["start"]
	median_length = np.median(df["length"])
	mean_length = np.mean(df["length"])
	print("Median Length:", median_length)
	print("Mean Length:", mean_length)





def motif_enrichment( df, outdir, genome_GDDH13, f_all_TEs, highTE_list, flank_size ):

	filter_TEs( df, outdir, f_all_TEs, highTE_list )
	launch_bedtools_window( outdir, flank_size )
	
	add_GOs_from_bedtools( outdir, f"{outdir}/GDDH13_allTEs_intersectVSgenes.tsv")
	add_GOs_from_bedtools( outdir, f"{outdir}/panREPET_intersectVSgenes.tsv")

	"""
	get_fasta_for_AME( genome_GDDH13, outdir)
	# Extract random sequences for TE residuals / M.domestica
	extract_random_sequences(f"{outdir}/panREPET.fasta", 
					genome_GDDH13_pop, 
					f"{outdir}/panREPET_random_seqs.fasta", 
					f_genes_GDDH13)
	"""





def get_highTE( VS ):
	"""
		Select 36 consensus potentially involved in climate adaptation
		x consensus, depends on VS value
	"""

	if VS == "noVS":
		f_highTE = f_TElibrary_nonpopStruc
	else:
		f_highTE = f"/mnt/Volume_Somia/Malus_PLEASURE_ANR/data_Anthony/TE_list_VS_gt_{VS}.txt"

	df_highTE = pd.read_csv(f_highTE, sep="\t", names=['TE_name'])
	#print(df_highTE)
	return df_highTE['TE_name']






def filter_clique_1(group):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement "M.domestica" ou "M.sylvestris"
	return all(group['accession'].isin(['M.domestica', 'M.sylvestris']))

def filter_clique_2(group):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement "M.domestica" ou "MdoR106A"
	return all(group['accession'].isin(['M.domestica', 'MdoR106A']))

#def filter_clique_3(group):
	# Vérifier si toutes les valeurs dans 'accession' sont strictement "M.domestica" ou "M.orientalis"
	#return all(group['accession'].isin(['M.domestica', 'M.orientalis']))







def AME_main( df_pop_structure, df_non_pop_structure, VS, flank_size ):

	out_dir = f"motif_enrichment_{VS}_{flank_size}"

	os.makedirs(out_dir, exist_ok=True)

	os.makedirs(f"{out_dir}/pop", exist_ok=True)
	os.makedirs(f"{out_dir}/pop/GDDH13_Msyl_specific", exist_ok=True)

	os.makedirs(f"{out_dir}/TEresiduals", exist_ok=True)
	os.makedirs(f"{out_dir}/TEresiduals/GDDH13specific", exist_ok=True)
	os.makedirs(f"{out_dir}/TEresiduals/nonGDDH13specific", exist_ok=True)
	os.makedirs(f"{out_dir}/TEresiduals/GDDH13_Msyl_specific", exist_ok=True)
	os.makedirs(f"{out_dir}/TEresiduals/GDDH13_Rootstock_specific", exist_ok=True)

		
	TE_names_values = set(get_highTE(VS))

	"""
	# 1. TE pop structure / M.domestica
	df_pop_structure_Mdo = df_pop_structure[df_pop_structure["accession"] == "M.domestica"]
	motif_enrichment( df_pop_structure_Mdo, f"{out_dir}/pop", genome_GDDH13_pop, f_all_TEs_pop, None )

	# 1.1 TE pop structure / M.domestica/M.Sylvestris specific
	df_pop_structure_Mdo_Msyl = df_pop_structure.groupby('TE insertions').filter(filter_clique_1)
	motif_enrichment( df_pop_structure_Mdo_Msyl, f"{out_dir}/pop/GDDH13_Msyl_specific", genome_GDDH13_pop, f_all_TEs_pop, None )
	"""

	# 2. TE residuals / M.domestica : x consensus (high TE) ==> depends on VS value
	df_non_pop_struct_GDDH13 = df_non_pop_structure[df_non_pop_structure["accession"] == "M.domestica"]
	motif_enrichment( df_non_pop_struct_GDDH13, f"{out_dir}/TEresiduals", genome_GDDH13_non_pop, f_all_TEs_non_pop, TE_names_values, flank_size )

	"""
	# 4. TE residuals / M.domestica (unique copies)
	df_non_pop_struct_GDDH13specific = df_non_pop_struct_GDDH13[df_non_pop_struct_GDDH13["clique_size"] == 1]
	motif_enrichment( df_non_pop_struct_GDDH13specific, f"{out_dir}/TEresiduals/GDDH13specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, None )

	# 4.1 TE residuals / M.domestica (unique copies)
	df_non_pop_struct_GDDH13specific = df_non_pop_struct_GDDH13[df_non_pop_struct_GDDH13["clique_size"] == 1]
	motif_enrichment( df_non_pop_struct_GDDH13specific, f"{out_dir}/TEresiduals/high_te/GDDH13specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, TE_names_values )

	# 5. TE residuals / M.domestica (shared TE copies but not unique copies)
	df_non_pop_struct_nonGDDH13specific = df_non_pop_struct_GDDH13[df_non_pop_struct_GDDH13["clique_size"] != 1]
	motif_enrichment( df_non_pop_struct_nonGDDH13specific, f"{out_dir}/TEresiduals/nonGDDH13specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, None )

	# 5.1 TE residuals / M.domestica (shared TE copies but not unique copies)
	df_non_pop_struct_nonGDDH13specific = df_non_pop_struct_GDDH13[df_non_pop_struct_GDDH13["clique_size"] != 1]
	motif_enrichment( df_non_pop_struct_nonGDDH13specific, f"{out_dir}/TEresiduals/high_te/nonGDDH13specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, TE_names_values )
	
	# 6. TE residuals / M.domestica/M.Sylvestris specific
	df_non_pop_struct_Mdo_Msyl = df_non_pop_structure.groupby('TE insertions').filter(filter_clique_1)
	df_non_pop_struct_Mdo_Msyl = df_non_pop_struct_Mdo_Msyl[df_non_pop_struct_Mdo_Msyl["accession"] == "M.domestica"]
	motif_enrichment( df_non_pop_struct_Mdo_Msyl, f"{out_dir}/TEresiduals/GDDH13_Msyl_specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, None )

	# 6.1 TE residuals / M.domestica/M.Sylvestris specific
	df_non_pop_struct_Mdo_Msyl = df_non_pop_structure.groupby('TE insertions').filter(filter_clique_1)
	df_non_pop_struct_Mdo_Msyl = df_non_pop_struct_Mdo_Msyl[df_non_pop_struct_Mdo_Msyl["accession"] == "M.domestica"]
	motif_enrichment( df_non_pop_struct_Mdo_Msyl, f"{out_dir}/TEresiduals/high_te/GDDH13_Msyl_specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, TE_names_values )

	# 7. TE residuals / M.domestica/Rootstock specific
	df_non_pop_struct_Mdo_Rootstock = df_non_pop_structure.groupby('TE insertions').filter(filter_clique_2)
	df_non_pop_struct_Mdo_Rootstock = df_non_pop_struct_Mdo_Rootstock[df_non_pop_struct_Mdo_Rootstock["accession"] == "M.domestica"]
	motif_enrichment( df_non_pop_struct_Mdo_Rootstock, f"{out_dir}/TEresiduals/GDDH13_Rootstock_specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, None )

	# 7.1 TE residuals / M.domestica/Rootstock specific
	df_non_pop_struct_Mdo_Rootstock = df_non_pop_structure.groupby('TE insertions').filter(filter_clique_2)
	df_non_pop_struct_Mdo_Rootstock = df_non_pop_struct_Mdo_Rootstock[df_non_pop_struct_Mdo_Rootstock["accession"] == "M.domestica"]
	motif_enrichment( df_non_pop_struct_Mdo_Rootstock, f"{out_dir}/TEresiduals/high_te/GDDH13_Rootstock_specific", genome_GDDH13_non_pop, f_all_TEs_non_pop, TE_names_values )
	"""





if __name__ == "__main__":


	df_pop_structure = reformat( f_pop_structure, f_unique_popStruc, dir_panTEannot_pop_structure, "TE_structure" )
	print(f"TE_pop_structure:\n{df_pop_structure}\n")

	df_non_pop_structure = reformat( f_non_pop_structure, f_unique_nonpopStruc, dir_panTEannot_non_pop_structure, "TE_residuals" )
	print(f"\nTE_residuals:\n{df_non_pop_structure}\n")
	print(df_non_pop_structure.columns)

#------ TORM
# 	df = pd.concat([df_pop_structure, df_non_pop_structure])
# 	df["seq_ID"] = "::" + df["chromosome"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
# 	df = df[["seq_ID", "consensus", "copy_name", "TE insertions"]]
# 	df_TE_highCopies = pd.read_csv( f_highCopiesCount, sep="\t" )
# 	df = df.merge(df_TE_highCopies, on="consensus", how="left")
# 	df.to_csv("panREPET.tsv", sep="\t", index=False)

# 	df_classif = pd.read_csv(f_classif, sep="\t", names=["consensus", "length", "classif"])
# 	print(f"\ndf_classif:\n{df_classif}\n")
# 	df_classif["consensus"] = df_classif["consensus"].str.replace(">", "")
# 	print(f"\ndf_classif:\n{df_classif}\n")
# 	df_highCopiesCount = pd.read_csv(f_highCopiesCount, sep="\t")
# 	print(df_highCopiesCount)
# 	df_highCopiesCount = pd.merge(df_highCopiesCount, df_classif, on="consensus", how="left")
# 	print(df_highCopiesCount)
# 	df_highCopiesCount.to_csv(f_highCopiesCount+".classif", sep="\t", index=False)
#------ TORM

# 	add_flowering()
# 	add_GOs()

	
# 	#--------------------------------------#
# 	#       orders_per_panCompartment      #
# 	#--------------------------------------#
# 	orders_filename = "orders_per_panCompartment_Anthony.png"
# 	if not os.path.exists( orders_filename ):

# 		fig, axes = plt.subplots(nrows=2, ncols=6, figsize=(30, 15))
# 		plt.subplots_adjust(wspace=0.2, hspace=0.15)
# 		piechart_panCompVSTEfamily( df_pop_structure, f_TElibrary_popStruc, axes, 0, "TE-core-structure" )
# 		piechart_panCompVSTEfamily( df_non_pop_structure, f_TElibrary_nonpopStruc, axes, 1, "TE-residuals" )
# 		plt.savefig( orders_filename, dpi=400, bbox_inches='tight' )
	

	#-----------------------#
	#       clustermap      #
	#-----------------------#
# 	clustermap_filename = "clustermap_panREPET_Anthony.png"
# 	if not os.path.exists( clustermap_filename ):

# 		fig, axs = plt.subplots(1,2)
# 		plt.subplots_adjust(wspace=0, hfspace=0)
# 		clustermap( df_pop_structure, "TE-core-structure", axs[0] )
# 		clustermap( df_non_pop_structure, "TE-residuals", axs[1] )
# 		fig.savefig( clustermap_filename, dpi=500, bbox_inches='tight' )
		


	#-----------------------#
	#   intersect with GO   #
	#-----------------------#

	window = 2000
	intersect_GO( df_pop_structure, "../../output/TE_structure",  window )
	intersect_GO( df_non_pop_structure, "../../output/TE_residuals", window )


# 	goea1_filename = "goea_core_vs_residuals.png"
# 	if not os.path.exists( goea1_filename ):
# 		fig, axes = plt.subplots(2, 2, figsize=(20, 10))

# 		# Case 1 : genes inside TEs
# 		GOEA( "TE_core_structure/panREPET.bed.intersectVSgenes.sameAss.GO", True, axes[0,0], "a", "TE_pop_structure - Genes in TE", "M.domestica" )
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", True, axes[1,0], "c", "TE_residuals - Genes in TE", "M.domestica" )

# 		# Case 2 : genes around TEs
# 		GOEA( "TE_core_structure/panREPET.bed.intersectVSgenes.sameAss.GO", False, axes[0,1], "b", "TE_pop_structure - 2kbp Upstream and Downstream TE-Associated Genes", "M.domestica" )
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", False, axes[1,1], "d", "TE_residuals - 2kbp Upstream and Downstream TE-Associated Genes", "M.domestica" )

# 		fig.savefig( goea1_filename )


# 	goea2_filename = "goea_wild_vs_cultivar.png"
# 	if not os.path.exists( goea2_filename ):
# 		fig, axes = plt.subplots(2, 2, figsize=(20, 15))

# 		#\textcolor{green}{\textbf{Wild}}

# 		# genes in TE
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", True, axes[0,0], "a", "Genes in TE \n Wild (M. sieversii)", "wild" ) #wild (M.sieversii/M.orientalis/M.sylvestris) I choose M.sieversii
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", True, axes[1,0], "c", "Genes in TE \n Cultivar (M. domestica)", "cultivar" ) #cultivar (M.domestica/MdoR106A) I choose M.domestica
# 		# 2kbp Upstream and Downstream TE-Associated Genes
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", False, axes[0,1], "b", "2kbp Upstream and Downstream TE-Associated Genes \n Wild (M. sieversii)", "wild" )
# 		GOEA( "TE_residuals/panREPET.bed.intersectVSgenes.sameAss.GO", False, axes[1,1], "d", "2kbp Upstream and Downstream TE-Associated Genes \n Cultivar (M. domestica)", "cultivar" )

# 		fig.savefig( goea2_filename )



	#-----------------------#
	#       upset plot      #
	#-----------------------#

	#check_both_annotation( )

	upset_filename = "../../Figures/upset_plot_Anthony.png"
	if os.path.exists( upset_filename ):

# 		add_gene_selection( df_pop_structure, "../../output/TE_structure" )
# 		add_gene_selection( df_non_pop_structure, "../../output/TE_residuals" )
	
		df_pop_structure = add_genes_presence_absence( df_pop_structure, "../../output/TE_structure" )
		print(f"\nAdd genes, TE_pop_structure:\n{df_pop_structure}\n")
		df_non_pop_structure = add_genes_presence_absence( df_non_pop_structure, "../../output/TE_residuals" )
		print(f"\nAdd genes, TE_residuals:\n{df_non_pop_structure}\n")

		fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,35))
		plt.subplots_adjust(wspace=0, hspace=0.05)
		upset_plot( df_pop_structure,  axes[0], "a - TE-structure" ) #"TE population structure" )
		upset_plot( df_non_pop_structure, axes[1], "b - TE-residuals" ) #TE non population structure" )
		plt.savefig( upset_filename, dpi=400, bbox_inches='tight' )
		#plt.savefig( upset_filename + ".svg", dpi=400, bbox_inches='tight', format="svg" )

	#---------------------------------------------------#
	#       check genes: Xilong's stress gene list      #
	#---------------------------------------------------#

		
	#check_Xilong_genes( "TE_core_structure" )
	#check_Xilong_genes( "TE_residuals" )


	#---------------------------------------#
	#       GOEA plot per species pair      #
	#---------------------------------------#

	
	#for TE_type in ["TE_core_structure", "TE_residuals"]:
	#	for sp in ['wild', 'cultivar', "M.domestica", "MdoR106A"]: #  'M.sieversii', 'M.orientalis', 'M.sylvestris', 
	#		for order in ["LTR_Gypsy", "all", "LTR", "LTR_Copia", "LTR_Gypsy", "TIR", "noCat", "LINE", "Helitron"]:
	#			f = f"{TE_type}/panREPET.bed.intersectVSgenes.sameAss.GO"
	#			GOEA( f"{TE_type}/panREPET.bed.intersectVSgenes.sameAss.GO", "both", None, None, None, sp, order )
	
	

# 	goea3_filename = "GOEA.png"
# 	if not os.path.exists( goea3_filename ):

# 		fig, axes = plt.subplots(2, 3, figsize=(25, 13)) #40

# 		GOEA_plot( axes[0,0], "a", "TE_core_structure", "all" )
# 		GOEA_plot( axes[0,1], "b", "TE_residuals", "all" )
# 		GOEA_plot( axes[0,2], "c", "TE_residuals", "LTR_Copia" )
# 		GOEA_plot( axes[1,0], "d", "TE_residuals", "LTR_Gypsy" )
# 		GOEA_plot( axes[1,1], "e", "TE_residuals", "LTR" )
# 		GOEA_plot( axes[1,2], "f", "TE_residuals", "TIR" )
# 		# and noCat ?

# 		fig.savefig( goea3_filename, bbox_inches='tight' )
	


	#---------------------------------------#
	#            motif_enrichment           #
	#---------------------------------------#

	#flank_size = 2000
	#for VS in ['noVS', '0.05', '0.20', '0.35', '0.50', '0.80']: 
	#	AME_main( df_pop_structure, df_non_pop_structure, VS, flank_size)
	



