__author__ = 'theep_000'

import math
import logging
import os
import tempfile
import sys
sys.path.append('/root/hg19_exons_introns.gtf')
import __builtin__

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

logger = logging.getLogger(__file__)


import archer_core.app.bedtools as bedtools


class Gene:
    """format for gene object that contains collected data
    """
    list = __builtin__.list
    percent_lower = 0
    base_a_l = []
    base_c_l = []
    base_g_l = []
    base_t_l = []
    base_a_track = 0
    base_c_track = 0
    base_g_track = 0
    base_t_track = 0
    total_base_l = []
    NucleotideTrackLength = 0
    stndard_deviation_tracklength = 0
    ExonTrackLength = 0
    exonstandard_deviation_tracklength = 0
    IntronTrackLength = 0
    intronstandard_deviation_tracklength = 0
    exon_l = []
    exonlength = 0
    intron_l = []
    intronlength = 0
    percent_intronic = 0
    percent_exonic = 0
    average_exon = 0
    standard_deviation_exon = 0
    average_intron = 0
    standard_deviation_intron = 0

    LowerCasePercent = 0

    def __init__(self, sequence=None, GCpercent=None, LowerCasePercent=None, NucleotideTracklength=None,
                 EighteenMer_Count=None, ExonLength=None, IntronLength=None, TotalLength=None, PercentIntronic=None, PercentExonic=None):
        self.base_a_l = []
        self.base_c_l = []
        self.base_g_l = []
        self.base_t_l = []
        self.base_a_track = 0
        self.base_c_track = 0
        self.base_g_track = 0
        self.base_t_track = 0

        self.sequence = sequence
        self.GC_percent = GCpercent
        self.LowerCasePercent = LowerCasePercent
        self.NucleotideTrackLength = NucleotideTracklength
        self.eighteenMers = {}
        self.exonlength = ExonLength
        self.intronlength = IntronLength
        self.totallength = TotalLength
        self.percent_intronic = PercentIntronic
        self.percent_exonic = PercentExonic

    def __str__(self):
        return "\t".join([str(x) for x in [self.sequence, self.GC_percent,
                                           self.LowerCasePercent, self.NucleotideTrackLength,
                                           self.eighteenmer_count, self.exonlength,
                                           self.intronlength, self.totallength,
                                           self.percent_intronic, self.percent_exonic]])

    def GC_Content_entgene(self):
        """determines percentage of bases which are G or C in the entire gene
        """
        gc_counter = 0
        total_base = 0
        for base in (self.sequence + " "):
            if '>' not in self.sequence:
                if base != '\n':
                    total_base += 1
                if base == 'g':
                    gc_counter += 1
                if base == 'c':
                    gc_counter += 1
                if base == 'G':
                    gc_counter += 1
                if base == 'C':
                    gc_counter += 1
        self.GC_percent = float(gc_counter) / float(total_base)
        self.totallength = total_base
        return self.GC_percent

    def TotalLength_entgene(self):
        """gives the amount of base pairs in the gene
        """
        self.total_count = 0
        for base in (self.sequence + " "):
            if '>' not in self.sequence:
                for base in self.sequence:
                    self.total_count += 1

    def Lower_case_count_entgene(self):
        """gives amount of base pairs which are lower case in the gene
        """
        lower_count = 0
        total = 0
        for base in (self.sequence + " "):
            if '>' not in self.sequence:
                for base in self.sequence:
                    total += 1
                    if base.islower():
                        lower_count += 1
        self.percent_lower = float(lower_count) / float(total)
        return self.percent_lower

    def get_TrackLength_a(self):
        """gives the length of all polynucleotide tracks for Base A
        """
        prev_character = "A" or "a"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_a_l.append(track_count)
                track_count = 0
            while 0 in self.base_a_l:
                    self.base_a_l.remove(0)
        self.base_a_track = sum(self.base_a_l) / float(len(self.base_a_l))
        prev_character = base.upper()

    def get_TrackLength_g(self):
        """gives the length of all polynucleotide tracks for Base G
        """
        prev_character = "G" or "g"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_g_l.append(track_count)
                track_count = 0
            while 0 in self.base_g_l:
                    self.base_g_l.remove(0)
        self.base_g_track = sum(self.base_g_l) / float(len(self.base_g_l))
        prev_character = base.upper()
        return self.base_g_l

    def get_TrackLength_t(self):
        """gives the length of all polynucleotide tracks for Base T
        """
        prev_character = "T" or "t"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_t_l.append(track_count)
                track_count = 0
            while 0 in self.base_t_l:
                    self.base_t_l.remove(0)
        self.base_t_track = sum(self.base_t_l) / float(len(self.base_t_l))
        prev_character = base.upper()
        return self.base_a_l

    def get_TrackLength_c(self):
        """gives the length of all polynucleotide tracks for Base T
        """
        prev_character = "C" or "c"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_c_l.append(track_count)
                track_count = 0
            while 0 in self.base_c_l:
                    self.base_c_l.remove(0)
        self.base_c_track = sum(self.base_c_l) / float(len(self.base_c_l))
        prev_character = base.upper()

    def get_Track_average(self):
        """gives the lengths of the polynucleotide tracks for all bases and standard deviation
        """
        prev_character = "G" or "g"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_g_l.append(track_count)
                track_count = 0
            while 0 in self.base_g_l:
                    self.base_g_l.remove(0)
        prevC_character = "C" or "c"
        trackC_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevC_character:
                trackC_count += 1
            elif base.upper() != prevC_character:
                self.base_c_l.append(trackC_count)
                trackC_count = 0
            while 0 in self.base_c_l:
                    self.base_c_l.remove(0)
        prevA_character = "A" or "a"
        trackA_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevA_character:
                trackA_count += 1
            elif base.upper() != prevA_character:
                self.base_a_l.append(trackA_count)
                trackA_count = 0
            while 0 in self.base_a_l:
                    self.base_a_l.remove(0)
        prevT_character = "T" or "t"
        trackT_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevT_character:
                trackT_count += 1
            elif base.upper() != prevT_character:
                self.base_t_l.append(trackT_count)
                trackT_count = 0
            while 0 in self.base_t_l:
                    self.base_t_l.remove(0)
        self.total_base_l = self.base_a_l + self.base_c_l + self.base_g_l + self.base_t_l
        self.NucleotideTrackLength = sum(self.total_base_l) / (len(self.total_base_l))
        variance_tracklength = sum(map(lambda x: (x - self.NucleotideTrackLength) ** 2, self.total_base_l))
        self.standard_deviation_tracklength = math.sqrt((variance_tracklength) / (len(self.total_base_l) - 1))

    def get_TrackLength_exon(self):
        """gives the lengths of the polynucleotide tracks for exons and standard deviation
        """
        for character in exon:
            prev_character = "G" or "g"
            track_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prev_character:
                    track_count += 1
                elif base.upper() != prev_character:
                    self.base_g_l.append(track_count)
                    track_count = 0
                while 0 in self.base_g_l:
                    self.base_g_l.remove(0)
            prevC_character = "C" or "c"
            trackC_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevC_character:
                    trackC_count += 1
                elif base.upper() != prevC_character:
                    self.base_c_l.append(trackC_count)
                    trackC_count = 0
                while 0 in self.base_c_l:
                    self.base_c_l.remove(0)
            prevA_character = "A" or "a"
            trackA_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevA_character:
                    trackA_count += 1
                elif base.upper() != prevA_character:
                    self.base_a_l.append(trackA_count)
                    trackA_count = 0
                while 0 in self.base_a_l:
                    self.base_a_l.remove(0)
            prevT_character = "T" or "t"
            trackT_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevT_character:
                    trackT_count += 1
                elif base.upper() != prevT_character:
                    self.base_t_l.append(trackT_count)
                    trackT_count = 0
                while 0 in self.base_t_l:
                    self.base_t_l.remove(0)
        exon_base_l = self.base_a_l + self.base_c_l + self.base_g_l + self.base_t_l
        self.ExonTrackLength = sum(exon_base_l) / (len(self.exon_base_l))
        exonvariance_tracklength = sum(map(lambda x: (x - self.ExonTrackLength) ** 2, exon_base_l))
        self.exonstandard_deviation_tracklength = math.sqrt((exonvariance_tracklength) / (len(self.exon_base_l) - 1))

    def get_TrackLength_intron(self):
        """gives the lengths of the polynucleotide tracks for introns and standard deviation
        """
        for character in intron:
            prev_character = "G" or "g"
            track_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prev_character:
                    track_count += 1
                elif base.upper() != prev_character:
                    self.base_g_l.append(track_count)
                    track_count = 0
                while 0 in self.base_g_l:
                    self.base_g_l.remove(0)
            prevC_character = "C" or "c"
            trackC_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevC_character:
                    trackC_count += 1
                elif base.upper() != prevC_character:
                    self.base_c_l.append(trackC_count)
                    trackC_count = 0
                while 0 in self.base_c_l:
                    self.base_c_l.remove(0)
            prevA_character = "A" or "a"
            trackA_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevA_character:
                    trackA_count += 1
                elif base.upper() != prevA_character:
                    self.base_a_l.append(trackA_count)
                    trackA_count = 0
                while 0 in self.base_a_l:
                    self.base_a_l.remove(0)
            prevT_character = "T" or "t"
            trackT_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevT_character:
                    trackT_count += 1
                elif base.upper() != prevT_character:
                    self.base_t_l.append(trackT_count)
                    trackT_count = 0
                while 0 in self.base_t_l:
                    self.base_t_l.remove(0)
        intron_base_l = self.base_a_l + self.base_c_l + self.base_g_l + self.base_t_l
        self.IntronTrackLength = sum(intron_base_l) / (len(self.intron_base_l))
        intronvariance_tracklength = sum(map(lambda x: (x - self.IntronTrackLength) ** 2, intron_base_l))
        self.intronstandard_deviation_tracklength = math.sqrt((intronvariance_tracklength) /
                                                              (len(self.intron_base_l) - 1))

    def get_eighteen(self):
        k = 18
        for i in range(len(self.sequence) - k + 1):
            EighteenMer_Count = self.sequence[i:i + k]
            if EighteenMer_Count in self.eighteenMers:
                self.eighteenMers[EighteenMer_Count] += 1
            else:
                self.eighteenMers[EighteenMer_Count] = 1
        return self.eighteenMers[EighteenMer_Count]

    """def get_exleng(self):
        exon_count = 0
        for base in (self.sequence):
            if base == exon:
                exon_count += 1
            if base != exon:
                self.exon_l.append(exon_count)
                exon_count = 0
        self.exonlength = sum(int(exon_count) for exon_count in self.exon_l)

    def get_inleng(self):
        intron_count = 0
        for base in (self.sequence):
            if base == intron:
                intron_count += 1
            if base != intron:
                self.intron_l.append(intron_count)
                intron_count = 0
        self.intronlength = sum(int(intron_count) for intron_count in self.intron_l)"""

    def get_totalleng(self):
        """Find the exon/intron lengths, percents, averages, and standard deviations
        """
        length_count = 0
        for base in (self.sequence):
            if base:
                length_count += 1
        exon_count = 0
        for base in (self.sequence):
            if base == exon:
                exon_count += 1
            if base != exon:
                self.exon_l.append(exon_count)
                exon_count = 0
        self.exonlength = sum(int(exon_count) for exon_count in self.exon_l)
        intron_count = 0
        for base in (self.sequence):
            if base == intron:
                intron_count += 1
            if base != intron:
                self.intron_l.append(intron_count)
                intron_count = 0
        self.intronlength = sum(int(intron_count) for intron_count in self.intron_l)
        self.percent_exonic = ((self.exonlength) / (length_count)) * 100
        self.percent_intronic = ((self.intronlength) / (length_count)) * 100
        self.average_exon = sum(self.exon_l) * 1.0 / len(self.exon_l)
        variance_exon = sum(map(lambda x: (x - self.average_exon) ** 2, self.exon_l))
        self.standard_deviation_exon = math.sqrt(self.average_exon(variance_exon))
        self.average_intron = sum(self.intron_l) * 1.0 / len(self.intron_l)
        variance_intron = sum(map(lambda x: (x - self.average_intron) ** 2, self.intron_l))
        self.standard_deviation_intron = math.sqrt(variance_intron / (len(self.intron_l) - 1))


def extract_sequence_from_genome(chromosome, start, stop,
                                human_genome = os.environ["ARCHER_GENOME"]):

    """Extracts a sequence from the human genome given a chromosom, start and stop

    :param str chromosome: Chromsome to extract from
    :param int start: Bed coordinate start of the sequence
    :parm int stop: Bed coordinate stop of the sequence
    """

    # Create a dummy bed line to use fastaFromBed of bedtools
    bed_lines = ["\t".join([chromosome, str(start), str(stop), "Name", "0", "+"])]

    temp_output_fasta = tempfile.NamedTemporaryFile()
    result = bedtools.extract_fasta_using_bed(bed_lines, human_genome, temp_output_fasta.name)

    # Grab the resulting output file
    output_file = result['output_file']
    output_fasta = [x.strip() for x in output_file.readlines()]

    #print output_fasta
    # The sequence will be the second line in the result
    return \
        "".join(output_fasta[1:]).replace("\n", "").strip()

def only_exon(bed):
    f1 = open('hg19 file', 'r')
    f2 = open('exon', 'w')
    f4 = open('hg19.exons_introns.gtf', 'r')

    for line in f4:
        if gene_coordinate_start in f1 <= exon_end in f4:
            f2.write(line)

    f1.close()
    f2.close()
    f4.close()

"""def GC_Content_entgene(extracted_seq):
    ex_seq = extracted_seq
    gc_counter = 0
    total_base = 0
    if '>' not in ex_seq:
        for base in ex_seq:
            if base != '\n':
                total_base += 1
            if base == 'g':
                gc_counter += 1
            if base == 'c':
                gc_counter += 1
            if base == 'G':
                gc_counter += 1
            if base == 'C':
                gc_counter += 1
    GC_percent = float(gc_counter) / float(total_base)
    return GC_percent


def TotalLength_entgene(extracted_seq):
    ex_seq = extracted_seq
    total_count = 0
    if '>' not in ex_seq:
        for base in ex_seq:
            total_count += 1
    return total_count


def Lower_case_count_entgene(extraced_seq):
    ex_seq = extraced_seq
    lower_count = 0
    total = 0
    if '>' not in ex_seq:
        for base in ex_seq:
            total += 1
            if base.islower():
                lower_count += 1
    percent_lower = float(lower_count) / float(total)
    return percent_lower"""


def bed_analyzer(bed):
    """the main process of Gene_Characterize
    """
    gene_list = []
    gene_num = 0
    for line in open(bed, 'r'):
        gene_num += 1
        gene_name = gene_num
        bed_data = line.split()
        single_gene = extract_sequence_from_genome(bed_data[0], int(bed_data[1]), int(bed_data[2]))
        single_gene_clone = single_gene
        #gene_lower = Lower_case_count_entgene(single_gene)
        #gene_gc = GC_Content_entgene(single_gene_clone)
        #totalleng = TotalLength_entgene(single_gene)
        gene_name = Gene(single_gene_clone)
        gene_list.append(gene_name)

    return gene_list


def parseCmdlineParams(arg_list=sys.argv):
    """Parses commandline arguments.

    :param list arg_list: Arguments to parse. Default is argv when called from the
    command-line.
    """
    # Create instance of ArgumentParser
    options_parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    options_parser.add_argument('--bed', '-b', type=str,
                                help='Bed file to analyze', required=True)

    # Parse options
    opts = options_parser.parse_args(args=arg_list[1:])

    return opts


def main(argv):
    opts = parseCmdlineParams(argv)
    gene = bed_analyzer(opts.bed)
    for curr in gene:
        print curr

if __name__ == "__main__":
    main(sys.argv)
