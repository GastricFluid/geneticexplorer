__author__ = 'theep_000'

import math
import logging
import os
import tempfile
import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

logger = logging.getLogger(__file__)


import archer_core.app.bedtools as bedtools


class Gene:
    """format for gene object that contains collected data
    """

    base_a_list = []
    base_c_list = []
    base_g_list = []
    base_t_list = []
    base_a_track = 0
    base_c_track = 0
    base_g_track = 0
    base_t_track = 0
    total_base_list = []
    NucleotideTrackLength = 0
    stndard_deviation_tracklength = 0
    ExonTrackLength = 0
    exonstandard_deviation_tracklength = 0
    IntronTrackLength = 0
    intronstandard_deviation_tracklength = 0
    eighteenmer_count = {}
    exon_list = []
    exonlength = 0
    intron_list = []
    intronlength = 0
    totallength = 0
    percent_intronic = 0
    percent_exonic = 0
    average_exon = 0
    standard_deviation_exon = 0
    average_intron = 0
    standard_deviation_intron = 0

    GC_percent = 0
    LowerCasePercent = 0

    def __init__(self, sequence, GCpercent, LowerCasePercent, NucleotideTracklength,
                 EighteenMer_Count, ExonLength, IntronLength, TotalLength, PercentIntronic, PercentExonic):
        self.base_a_list = []
        self.base_c_list = []
        self.base_g_list = []
        self.base_t_list = []
        self.base_a_track = 0
        self.base_c_track = 0
        self.base_g_track = 0
        self.base_t_track = 0

        self.sequence = sequence
        self.GC_percent = GCpercent
        self.LowerCasePercent = LowerCasePercent
        self.NucleotideTrackLength = NucleotideTracklength
        self.eighteenmer_count = EighteenMer_Count
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

    def get_GCperc(self):
        return self.GC_percent

    def get_Lower(self):
        return self.LowerCasePercent

    def get_TrackLength_a(self):
        prev_character = "A" or "a"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_a_list.append(track_count)
                track_count = 0
            while 0 in self.base_a_list:
                    self.base_a_list.remove(0)
        self.base_a_track = sum(self.base_a_list) / float(len(self.base_a_list))
        prev_character = base.upper()

    def get_TrackLength_g(self):
        prev_character = "G" or "g"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_g_list.append(track_count)
                track_count = 0
            while 0 in self.base_g_list:
                    self.base_g_list.remove(0)
        self.base_g_track = sum(self.base_g_list) / float(len(self.base_g_list))
        prev_character = base.upper()
        return self.base_g_list

    def get_TrackLength_t(self):
        prev_character = "T" or "t"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_t_list.append(track_count)
                track_count = 0
            while 0 in self.base_t_list:
                    self.base_t_list.remove(0)
        self.base_t_track = sum(self.base_t_list) / float(len(self.base_t_list))
        prev_character = base.upper()

    def get_TrackLength_c(self):
        prev_character = "C" or "c"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_c_list.append(track_count)
                track_count = 0
            while 0 in self.base_c_list:
                    self.base_c_list.remove(0)
        self.base_c_track = sum(self.base_c_list) / float(len(self.base_c_list))
        prev_character = base.upper()

    def get_Track_average(self):
        """prev_character = "G" or "g"
        track_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prev_character:
                track_count += 1
            elif base.upper() != prev_character:
                self.base_g_list.append(track_count)
                track_count = 0
            while 0 in self.base_g_list:
                    self.base_g_list.remove(0)
        prevC_character = "C" or "c"
        trackC_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevC_character:
                trackC_count += 1
            elif base.upper() != prevC_character:
                self.base_c_list.append(trackC_count)
                trackC_count = 0
            while 0 in self.base_c_list:
                    self.base_c_list.remove(0)
        prevA_character = "A" or "a"
        trackA_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevA_character:
                trackA_count += 1
            elif base.upper() != prevA_character:
                self.base_a_list.append(trackA_count)
                trackA_count = 0
            while 0 in self.base_a_list:
                    self.base_a_list.remove(0)
        prevT_character = "T" or "t"
        trackT_count = 0
        for base in (self.sequence + " "):
            if base.upper() == prevT_character:
                trackT_count += 1
            elif base.upper() != prevT_character:
                self.base_t_list.append(trackT_count)
                trackT_count = 0
            while 0 in self.base_t_list:
                    self.base_t_list.remove(0)"""
        self.total_base_list = self.base_a_list + self.base_c_list + self.base_g_list + self.base_t_list
        """self.NucleotideTrackLength = sum(self.total_base_list) / (len(self.total_base_list))
        variance_tracklength = sum(map(lambda x: (x - self.NucleotideTrackLength) ** 2, self.total_base_list))
        self.standard_deviation_tracklength = math.sqrt((variance_tracklength) / (len(self.total_base_list) - 1))"""

    def get_TrackLength_exon(self):
        for character in exon:
            prev_character = "G" or "g"
            track_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prev_character:
                    track_count += 1
                elif base.upper() != prev_character:
                    self.base_g_list.append(track_count)
                    track_count = 0
                while 0 in self.base_g_list:
                    self.base_g_list.remove(0)
            prevC_character = "C" or "c"
            trackC_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevC_character:
                    trackC_count += 1
                elif base.upper() != prevC_character:
                    self.base_c_list.append(trackC_count)
                    trackC_count = 0
                while 0 in self.base_c_list:
                    self.base_c_list.remove(0)
            prevA_character = "A" or "a"
            trackA_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevA_character:
                    trackA_count += 1
                elif base.upper() != prevA_character:
                    self.base_a_list.append(trackA_count)
                    trackA_count = 0
                while 0 in self.base_a_list:
                    self.base_a_list.remove(0)
            prevT_character = "T" or "t"
            trackT_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevT_character:
                    trackT_count += 1
                elif base.upper() != prevT_character:
                    self.base_t_list.append(trackT_count)
                    trackT_count = 0
                while 0 in self.base_t_list:
                    self.base_t_list.remove(0)
        exon_base_list = self.base_a_list + self.base_c_list + self.base_g_list + self.base_t_list
        self.ExonTrackLength = sum(exon_base_list) / (len(self.exon_base_list))
        exonvariance_tracklength = sum(map(lambda x: (x - self.ExonTrackLength) ** 2, exon_base_list))
        self.exonstandard_deviation_tracklength = math.sqrt((exonvariance_tracklength) / (len(self.exon_base_list) - 1))

    def get_TrackLength_intron(self):
        for character in intron:
            prev_character = "G" or "g"
            track_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prev_character:
                    track_count += 1
                elif base.upper() != prev_character:
                    self.base_g_list.append(track_count)
                    track_count = 0
                while 0 in self.base_g_list:
                    self.base_g_list.remove(0)
            prevC_character = "C" or "c"
            trackC_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevC_character:
                    trackC_count += 1
                elif base.upper() != prevC_character:
                    self.base_c_list.append(trackC_count)
                    trackC_count = 0
                while 0 in self.base_c_list:
                    self.base_c_list.remove(0)
            prevA_character = "A" or "a"
            trackA_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevA_character:
                    trackA_count += 1
                elif base.upper() != prevA_character:
                    self.base_a_list.append(trackA_count)
                    trackA_count = 0
                while 0 in self.base_a_list:
                    self.base_a_list.remove(0)
            prevT_character = "T" or "t"
            trackT_count = 0
            for base in (self.sequence + " "):
                if base.upper() == prevT_character:
                    trackT_count += 1
                elif base.upper() != prevT_character:
                    self.base_t_list.append(trackT_count)
                    trackT_count = 0
                while 0 in self.base_t_list:
                    self.base_t_list.remove(0)
        intron_base_list = self.base_a_list + self.base_c_list + self.base_g_list + self.base_t_list
        self.IntronTrackLength = sum(intron_base_list) / (len(self.intron_base_list))
        intronvariance_tracklength = sum(map(lambda x: (x - self.IntronTrackLength) ** 2, intron_base_list))
        self.intronstandard_deviation_tracklength = math.sqrt((intronvariance_tracklength) /
                                                              (len(self.intron_base_list) - 1))

    def get_eighteen(self):
        k = 18
        for i in range(len(self.sequence) - k + 1):
            EighteenMer_Count = self.sequence[i:i + k]
            if EighteenMer_Count.haskey(EighteenMer_Count):
                EighteenMer_Count[EighteenMer_Count] += 1
            else:
                EighteenMer_Count[EighteenMer_Count] = 1
            for EighteenMer_Count, count in EighteenMer_Count.items():
                return self.eighteenmer_count

    def get_exleng(self):
        exon_count = 0
        for base in (self.sequence):
            if base == exon:
                exon_count += 1
            if base != exon:
                self.exon_list.append(exon_count)
                exon_count = 0
        self.exonlength = sum(int(exon_count) for exon_count in self.exon_list)

    def get_inleng(self):
        intron_count = 0
        for base in (self.sequence):
            if base == intron:
                intron_count += 1
            if base != intron:
                self.intron_list.append(intron_count)
                intron_count = 0
        self.intronlength = sum(int(intron_count) for intron_count in self.intron_list)

    def get_totalleng(self):
        length_count = 0
        for base in (self.sequence):
            if base:
                length_count += 1
        exon_count = 0
        for base in (self.sequence):
            if base == exon:
                exon_count += 1
            if base != exon:
                self.exon_list.append(exon_count)
                exon_count = 0
        self.exonlength = sum(int(exon_count) for exon_count in self.exon_list)
        intron_count = 0
        for base in (self.sequence):
            if base == intron:
                intron_count += 1
            if base != intron:
                self.intron_list.append(intron_count)
                intron_count = 0
        self.intronlength = sum(int(intron_count) for intron_count in self.intron_list)
        self.percent_exonic = ((self.exonlength) / (length_count)) * 100
        self.percent_intronic = ((self.intronlength) / (length_count)) * 100
        self.average_exon = sum(self.exon_list) * 1.0 / len(self.exon_list)
        variance_exon = sum(map(lambda x: (x - self.average_exon) ** 2, self.exon_list))
        self.standard_deviation_exon = math.sqrt(self.average_exon(variance_exon))
        self.average_intron = sum(self.intron_list) * 1.0 / len(self.intron_list)
        variance_intron = sum(map(lambda x: (x - self.average_intron) ** 2, self.intron_list))
        self.standard_deviation_intron = math.sqrt(variance_intron / (len(self.intron_list) - 1))


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


def GC_Content_entgene(extracted_seq):
    """determines percentage of bases which are G or C in the entire gene
    """
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
    """gives the amount of base pairs in the gene
    """
    ex_seq = extracted_seq
    total_count = 0
    if '>' not in ex_seq:
        for base in ex_seq:
            total_count += 1
    return total_count


def Lower_case_count_entgene(extraced_seq):
    """gives amount of base pairs which are lower case in the gene
    """
    ex_seq = extraced_seq
    lower_count = 0
    total = 0
    if '>' not in ex_seq:
        for base in ex_seq:
            total += 1
            if base.islower():
                lower_count += 1
    percent_lower = float(lower_count) / float(total)
    return percent_lower


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
        gene_lower = Lower_case_count_entgene(single_gene)
        gene_gc = GC_Content_entgene(single_gene_clone)
        totalleng = TotalLength_entgene(single_gene)
        gene_name = Gene("", gene_gc, gene_lower, 0, 0, 0, 0, totalleng, 0, 0)
        gene_list.append(gene_name)
        logger.info(gene_name.get_GCperc)
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
