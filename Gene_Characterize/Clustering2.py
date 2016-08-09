
import sys
sys.path.append('/root/gene_clusters/geneticexplorer/')

from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as npy
import gene_analyzer
# gene_list = []
# var_num = 0
from sys import argv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def listpopulator(bed_file):
    """creates list of genes for use in clustering
    """
    gene_list = gene_analyzer.bed_analyzer(bed_file)
    return gene_list


def gene_val_select(exclusion_dict, gene):
    """creates a list of values that were not excluded
    """
    gene_vals = []
    if exclusion_dict['GCcont'] != 1:
        gene_vals.append(gene.get_GCperc())
    if exclusion_dict['totleng'] != 1:
        gene_vals.append(gene.get_totalleng())
    if exclusion_dict['track'] != 1:
        gene_vals.append(gene.get_Track_average())
    if exclusion_dict['18mer'] != 1:
        gene_vals.append(gene.get_eighteen())
    if exclusion_dict['percex'] != 1:
        gene_vals.append(gene.get_percentexonic())
    if exclusion_dict['percint'] != 1:
        gene_vals.append(gene.get_percentintronic())
    if exclusion_dict['lower'] != 1:
        gene_vals.append(gene.get_Lower())
    print gene_vals
    return gene_vals


def clustermethod(gene_list, exclusion_dict):
    """clusters a list of genes using scipy linkage
    """
    array_val = 0
    array_list = []
    past_point = False
    for gene in gene_list:
        array_val += 1
        gene_vals = gene_val_select(exclusion_dict, gene)
        array_val = npy.fromiter(gene_vals, float)
        array_list.append(array_val)
        # loop produces a ndarray for every gene, containing it's data, for use in concatenation then clustering
        #for point in array_list:
        # loop concatenates all the data
        #current_point = point
        #if past_point is False:
        #    past_point = True
        #    past_data = point
        #else:
        #    data = npy.concatenate([past_data, current_point])
        #    past_data = data
        # this code is most likely uneeded ^^^^
    print array_list
    heirclust = linkage(array_list, "ward")
    print heirclust[0]
    return heirclust


def workflow(bed, exclusion_dict):
    """workflow given a bed file and dict of what to exclude which produces a clustered ndarray
    """
    clustered_array = clustermethod(listpopulator(bed), exclusion_dict)
    print clustered_array
    return clustered_array


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
        workflow(bed_file, exclusion_dict)

    else:
        print "please give a bed file"


if __name__ == "__main__":
    main(argv)
