
import sys
sys.path.append('/root/geneticexplorer/')

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as npy
import gene_analyzer
from sys import argv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def listpopulator(bed_file):
    """creates list of genes for use in clustering
    """
    gene_list = gene_analyzer.bed_analyzer(bed_file)
    genes = gene_list[1]
    names = gene_list[0]
    return (names, genes)


def gene_val_select(exclusion_dict, gene):
    """creates a list of values that were not excluded
    """
    gene_vals = []
    if exclusion_dict['GCcont'] != 1:
        gene_vals.append(gene.GC_Content_entgene())
    if exclusion_dict['totleng'] != 1:
        gene_vals.append(gene.totallength)
    if exclusion_dict['track'] != 1:
        gene.get_Track_average()
        gene_vals.append(gene.NucleotideTrackLength)
        gene_vals.append(gene.standard_deviation_tracklength)
    if exclusion_dict['18mer'] != 1:
        gene_vals.append(gene.get_eighteen())
    if exclusion_dict['percex'] != 1:
        gene_vals.append(gene.get_percentexonic())
    if exclusion_dict['percint'] != 1:
        gene_vals.append(gene.get_percentintronic())
    if exclusion_dict['lower'] != 1:
        gene_vals.append(gene.Lower_case_count_entgene())
    return gene_vals


def clustermethod(gene_list, exclusion_dict):
    """clusters a list of genes using scipy linkage
    """
    array_val = 0
    array_list = []
    for gene in gene_list:
        array_val += 1
        gene_vals = gene_val_select(exclusion_dict, gene)
        array_val = npy.fromiter(gene_vals, float)
        array_list.append(array_val)
    heirclust = linkage(array_list, "ward")
    return heirclust


def workflow(bed, exclusion_dict):
    """workflow given a bed file and dict of what to exclude which produces a clustered ndarray
    """
    gene_list = listpopulator(bed)
    names = gene_list[0]
    clustered_array = clustermethod(gene_list[1], exclusion_dict)
    return (names, clustered_array)


def visualization(clustered_data, gene_names):
    gene_names = gene_names
    plt.figure(figsize=(15, 20))
    plt.title('Gene Clusters')
    plt.ylabel('Distance')
    plt.xlabel('Genes')
    dendrogram(clustered_data, labels=gene_names, show_leaf_counts=False)
    plt.show()
    plt.savefig('rad_clusters.png')


def exclusion_builder(tracknum, totlengnum, GCcontnum, eighteenmernum,
                      percexnum, percintnum, lowernum):
    """creates a dictionary using command line arguments which dictates what data will be excluded
    """
    exclusion_dict = {}
    exclusion_dict['track'] = tracknum
    exclusion_dict['totleng'] = totlengnum
    exclusion_dict['GCcont'] = GCcontnum
    exclusion_dict['18mer'] = eighteenmernum
    exclusion_dict['percex'] = percexnum
    exclusion_dict['percint'] = percintnum
    exclusion_dict['lower'] = lowernum
    return exclusion_dict


def commandLine(arg_list=argv):
    """defines command line arguments
    """
    options_parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    options_parser.add_argument('-b', type=str, help='bedfile with genes to be clustered', action="store", dest='bed')
    options_parser.add_argument('-track', type=int, help='exclude track length', default=0, dest='track')
    options_parser.add_argument('-totleng', type=int, help='exclude total length', default=0, dest='totleng')
    options_parser.add_argument('-GCcont', type=int, help='exclude GCcontent', default=0, dest='GCcont')
    options_parser.add_argument('-18mer', type=int, help='exclude 18mer count', default=0, dest='eighteenmer')
    options_parser.add_argument('-percex', type=int, help='exclude percent exonic', default=0, dest='percex')
    options_parser.add_argument('-percint', type=int, help='exclude percent intronic', default=0, dest='percint')
    options_parser.add_argument('-lower', type=int, help='exclude percent lowercase', default=0, dest='lower')

    opts, unknown = options_parser.parse_known_args(args=arg_list)

    return opts


def main(args, args_parsed=None):
    """function which uses command line arguments and runs the workflow based off the given parameters
    """
    if args_parsed is not None:
        opts = args_parsed

    else:
        opts = commandLine(args)

    exclusion_dict = exclusion_builder(opts.track, opts.totleng, opts.GCcont, opts.eighteenmer,
                                       opts.percex, opts.percint, opts.lower)

    if opts.bed is not None:
        bed_file = opts.bed
        clustered_data = workflow(bed_file, exclusion_dict)
        visualization(clustered_data[1], clustered_data[0])

    else:
        print "please give a bed file"


if __name__ == "__main__":
    main(argv)
