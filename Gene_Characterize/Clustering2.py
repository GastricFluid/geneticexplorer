
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as npy
import gene_analyzer
gene_list = []
var_num = 0
from sys import argv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def listpopulator():
    """creates list of genes for use in clustering
    """
    #return gene_list
    pass


def clustermethod(gene_list):
    """clusters a list of genes using scipy linkage
    """
    array_val = 0
    array_list = []
    past_point = False
    for gene in gene_list:
        array_val += 1
        gene_vals = (gene.get_GCperc(), gene.get_Lower(), gene.get_totalleng())
        array_val = npy.fromiter(gene_vals, float)
        array_list.append(array_val)
        #loop produces a ndarray for every gene, containing it's data, for use in concatenation then clustering
    for point in array_list:
        #loop concatenates all the data
        current_point = point
        if past_point is False:
            past_point = True
            past_data = point
            pass
        else:
            data = npy.concatenate(past_data, current_point)
            past_data = data
    heirclust = linkage(data, "ward")
    return heirclust


def workflow(bed):
    clustermethod(listpopulator(bed))


def commandLine(arg_list=argv):
    """defines command line arguments
    """
    options_parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    options_parser.add_argument('-b', help='bedfile with genes to be clustered', action="store", dest='bed')

    opts, unkown = options_parser.parse_known_args(args=arg_list)


def main(args, args_parsed="None"):
    if args_parsed is not None:
        opts = args_parsed

    else:
        opts = commandLine(args)

    if opts.bed is True:
        bed_file = opts.bed

    if bed_file is True:
        workflow(bed_file)

    else:
        print "please give a bed file"

if __name__ == "__main__":
    main(argv)