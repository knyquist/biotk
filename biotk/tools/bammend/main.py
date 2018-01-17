#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import bammend as bm

def parse_args():
    """Parse command-line arguments"""
    summary = ('Remove pulses from reads in Pacbio Bam. Annotation indices \n'
               'are indexed from the beginning of the ZMW read (i.e. query \n'
               'indexing).')
    parser = argparse.ArgumentParser(prog='bammend',
                                     description=summary,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('subreads',
                        help='Path to subread bam')
    parser.add_argument('bammend_csv',
                        help=('Path to CSV with scheme \n'
                              '| ZMW | Annotation Start Index '
                              '| Annotation End Index |'))
    parser.add_argument('output_subreads',
                        help='Path to output bam')
    args = parser.parse_args()
    return args.subreads, args.bammend_csv, args.output_subreads

def main():
    """Bammend a subreadset."""
    read_bam_path, annotation_csv_path, out_bam_path = parse_args()
    bm.reject_basecalls(read_bam_path, annotation_csv_path, out_bam_path)


if __name__ == '__main__':
    main()
