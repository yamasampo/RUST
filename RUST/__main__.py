'''
script to run RUST as a cli tool 

Usage:
    RUST <command> [<args>...]

Options:
    -h --help   Show this screen.
    --version   Show version.


'''

import RUST
# import sys
# import time
# import traceback
import argparse

def parser():
    parser = argparse.ArgumentParser(description='RUST - Ribo-seq Unit Step Transform')
    parser.usage = 'RUST <command> [<args>]'

    subparsers = parser.add_subparsers(help='mode of analysis')

    amino = subparsers.add_parser('amino', help='Run RUST on each amino acid')
    amino.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    amino.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    amino.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    amino.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    amino.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "amino")
    amino.set_defaults(mode='amino')
    # amino_args = amino.parse_args()

    codon = subparsers.add_parser('codon', help='Run RUST on each codon')
    codon.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    codon.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    codon.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    codon.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    codon.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "codon")
    codon.set_defaults(mode='codon')
    # codon_args = codon.parse_args()

    nucleotide = subparsers.add_parser('nucleotide', help='Run RUST on each nucleotide')
    nucleotide.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    nucleotide.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    nucleotide.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    nucleotide.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    nucleotide.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "nucleotide")
    nucleotide.set_defaults(mode='nucleotide')
    # nucleotide_args = nucleotide.parse_args()

    dipeptide = subparsers.add_parser('dipeptide', help='Run RUST on each dipeptide')
    dipeptide.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    dipeptide.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    dipeptide.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    dipeptide.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    dipeptide.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "dipeptide")
    dipeptide.set_defaults(mode='dipeptide')
    # dipeptide_args = dipeptide.parse_args()

    tripeptide = subparsers.add_parser('tripeptide', help='Run RUST on each tripeptide')
    tripeptide.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    tripeptide.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    tripeptide.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    tripeptide.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    tripeptide.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "tripeptide")
    tripeptide.set_defaults(mode='tripeptide')
    # tripeptide_args = tripeptide.parse_args()

    predict = subparsers.add_parser('predict', help='Correlation between observed and predicted profiles from CDS start + 120 to CDS stop - 60')
    predict.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    predict.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    predict.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    predict.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    predict.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "predict_profiles")
    predict.set_defaults(mode='predict')
    # predict_args = predict.parse_args()

    synergy = subparsers.add_parser('synergy', help='Identifies tripeptides that are candidates for synergistic interactions')
    synergy.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    synergy.add_argument('--aa', help='path to file produced from "rust_amino"', required=True)
    synergy.add_argument('--tri', help='path to file produced from "rust_tripeptide"', required=True)
    synergy.add_argument('-P', '--Path', help='path to outputfile, default is "synergy"',default = "synergy")
    synergy.set_defaults(mode='synergy')
    # synergy_args = synergy.parse_args()

    plot = subparsers.add_parser('plot', help='Plot observed and predicted ribosome profiles')
    plot.add_argument('-t', '--transcriptome', help='fasta file of transcripts, CDS start and end may be provided on description line using tab separation e.g. >NM_0001  10  5000, otherwise it searches for longest ORF'', required=True')
    plot.add_argument('-a', '--alignment', help='sorted bam file of transcriptome alignments', required=True)
    plot.add_argument('-o', '--offset', help='nucleotide offset to A-site',type=int )
    plot.add_argument('-l', '--lengths', help='lengths of footprints included, for example 28:32 is 28,29,30,31,32')
    plot.add_argument('-P', '--Path', help='path to outputfile, default is "amino"',default = "plot")
    plot.add_argument('-i', '--identifier', help='Specific transcript to plot (Use of unique identifier is sufficient for example "NM_031946"', required=True)
    plot.add_argument('-r', '--rustfile', help='path to file produced from "rust_codon"', required=True)
    plot.set_defaults(mode='plot')
    # plot_args = plot.parse_args()


    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.3.0')


    args = parser.parse_args()
    return args


def main():
    import sys
    args = parser()
    print(args)
    print(sys.argv[1])
    print(args.mode)
    # if args.subparser_name == 'amino':
    #     RUST.run_amino.main(args)
    #     print('hello from amino')
    # elif args.subparser_name == 'codon':
    #     print('hello from codon')
    # elif args.subparser_name == 'nucleotide':
    #     print('hello from nucleotide')
    # elif args.subparser_name == 'dipeptide':
    #     print('hello from dipeptide')
    # elif args.subparser_name == 'tripeptide':
    #     print('hello from tripeptide')
    # elif args.subparser_name == 'predict':
    #     print('hello from predict')
    # elif args.subparser_name == 'synergy':
    #     print('hello from synergy')
    # elif args.subparser_name == 'plot':
    #     print('hello from plot')
    # else:
    #     print('hello from else')



