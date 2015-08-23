#!/usr/bin/env python

# Bedtools coverage processor for exon capture data!

from __future__ import print_function
from os import path
# The fasta parser.
import screed
import docopt
import sys
from pybedtools import BedTool
import itertools

__author__ = "lteasnail"

CLI_ARGS = """
USAGE:
bedtools_cover.py [-r REF] <BAM_FILES> ...

OPTIONS:
 -r REF         The reference fasta file you used as a reference in bfast. This
                code is designed to be used for single reference pipeline runs
                i.e. where you have used the same fasta file for each
                reference.

This bedtools coveage processor will analyse bam files to ascertain average
coverage, the number of exons with coverage and the proportion of exons with
coverage using bedtools.

All it needs to run is the bam files and the reference fasta file used in
bfast.

You need to have the python packages docopt, pybedtools and screed, hence
while you can get these to run on windows it's more straight forward to run
this script on linux.

e.g.

Single ref:

python bedtools_cover.py -r reference.fasta *.bam > output_table.txt

Species specific ref:

python bedtools_cover.py *.bam > output_table.txt

For this option you don't need to specify the reference fasta files but you do
need to have them in the same directory each respective BAM is in. Both the
reference fasta file and the bam file need to have the same name for
the species specific mode to work.

e.g Sphareospira_fraseri.bam and Sphareospira_fraseri.fasta


"""

# Creates a genome file from the reference fasta sequences the reads were
# mapped against.


def genome_file(reference):
    contig_lens = {}
    file_basename = path.splitext(reference)[0]
    new_filename = file_basename + '_genome.txt'
    with open(new_filename, 'w') as gen_fh:
        for seq in screed.open(reference):
            seq_name = seq.name
            length = len(seq.sequence)
            contig_lens[seq_name] = length
            print(seq_name, length, sep='\t', file=gen_fh)
    return contig_lens

# Converts the output of bedtools into a table with the average coverage
# for each contig.


def parse_bed_output(cov):
    def iter_cov(filename):
        with open(filename) as bedtools_output:
            for line in bedtools_output:
                cells = line.strip('\r\n').split('\t')
                yield cells
    for contig, cov_lines in itertools.groupby(iter_cov(cov.fn),
                                               lambda x: x[0]):
        bp_covs = []
        for line in cov_lines:
            bp_covs.append(int(line[2]))
        yield (contig, bp_covs)

# Calculates the average total coverage, the number of contigs with coverage,
# and the proportion of contigs with coverage.


def calc_stats(cov_table, genome_lens, output_name, bamname):
    def mean(lst):
        return sum(lst) / float(len(lst))
    stats = {}  # dict of contg: (mean, med, sd)
    for contig, covs in parse_bed_output(cov_table):
        contig_len = genome_lens[contig]
        for x in range(len(covs), contig_len):
            covs.append(0)
        stats[contig] = mean(covs)

    with open(output_name, 'w') as exon_means_file:
        sum_sums = 0
        bp_total = 0
        n_uncovered_exon = 0
        for contig, ctglen in sorted(genome_lens.items()):
            ctgmean = stats.get(contig, 0.0)
            if ctgmean == 0.0:
                n_uncovered_exon += 1
            sum_sums += ctgmean
            bp_total += 1
            print(contig, ctglen, ctgmean, sep='\t', file=exon_means_file)

    average_average_cov = sum_sums / float(bp_total)
    num_contigs_with_cov = len(genome_lens) - n_uncovered_exon
    prop_contigs_with_cov = 1 - (n_uncovered_exon / float(len(genome_lens)))
    print(bamname, average_average_cov, num_contigs_with_cov,
          prop_contigs_with_cov, sep='\t', file=sys.stdout)

# If I am being run as a script...
if __name__ == '__main__':
    opts = docopt.docopt(CLI_ARGS)
    reference = opts['-r']
    BAMfiles = opts['<BAM_FILES>']
# Run bedtools if Single reference run (i.e. all samples were mapped to the
# one reference)
    if reference is not None:
        print('producing reference genome file...', file=sys.stderr)
        genome_lens = genome_file(reference)
        print('Finished producing genome file', file=sys.stderr)
        print('Running bedtools for each BAM file', file=sys.stderr)
        ref_basename = path.splitext(reference)[0]
        ref_genome = '{}_genome.txt'.format(ref_basename)
        print('Bam_file_name', 'Average_exon_coverage',
              'Num_exons_with_coverage', 'Proportion_of_exons_with_coverage',
              sep='\t', file=sys.stdout)
        for bamfile in BAMfiles:
            print('Running bedtools for {}'.format(bamfile), file=sys.stderr)
            bed = BedTool(bamfile).bam_to_bed()
            cov_table = bed.genome_coverage(d=True, g=ref_genome)
            output_name = path.splitext(bamfile)[0]
            means_filename = output_name + '_mean_per_exon.txt'
            calc_stats(cov_table, genome_lens, means_filename, bamfile)

        print('Finished!', file=sys.stderr)

    else:
        # run bedtools if a reference per sample was used (i.e. Species
        # specific references)
        print('Bam_file_name', 'Average_exon_coverage',
              'Num_exons_with_coverage', 'Proportion_of_exons_with_coverage',
              sep='\t', file=sys.stdout)
        for bamfile in BAMfiles:
            bam_basename = path.splitext(bamfile)[0]
            ref_file = bam_basename + '.fasta'
            print('producing reference {} genome file...'.format(bam_basename),
                  file=sys.stderr)
            genome_lens = genome_file(ref_file)
            print('Finished producing genome file', file=sys.stderr)
            print('Running bedtools for {}.bam'.format(bam_basename),
                  file=sys.stderr)
            ref_genome = '{}_genome.txt'.format(bam_basename)
            print('Running bedtools for {}.bam'.format(bam_basename),
                  file=sys.stderr)
            bed = BedTool(bamfile).bam_to_bed()
            cov_table = bed.genome_coverage(d=True, g=ref_genome)
            means_filename = bam_basename + '_mean_per_exon.txt'
            calc_stats(cov_table, genome_lens, means_filename, bamfile)
        print('Finished!', file=sys.stderr)
