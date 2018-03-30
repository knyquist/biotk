#!/usr/bin/env python
import argparse
from argparse import RawTextHelpFormatter
import taulysis as ta
import utils as ut
import numpy as np
from numpy.lib import recfunctions as recfuncs
import logging
logging.basicConfig(level=logging.INFO)

def parse_args():
    """Parse command-line arguments"""
    summary = ('Estimate first-pass, second-pass, and \n'
               'rolling-circle taus for fixed-template \n'
               'sequencing experiments.')
    parser = argparse.ArgumentParser(prog='taulysis',
                                     description=summary,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('alignments',
                        help='Path to alignmentset')
    parser.add_argument('output_directory',
                        help='Path to output directory')
    parser.add_argument('--smrtbellsize',
                        default=None,
                        help=('Set SMRTBell length (in bases). Optional flag.\n'
                              'Defaults to trying to determine from the \n'
                              'reference name by pattern-matching \n'
                              "unrolled_circular_\n"
                              "[number of times unrolled]x_l[unrolled length]"))
    parser.add_argument('--first-adapter-hit',
                        default=None,
                        help=('Set template position of first adapter hit.\n'
                              'Used for selecting data to fit first-pass tau.\n'
                              'Optional flag. Default value is \n'
                              'smrtbellsize / 2.'))
    parser.add_argument('--second-adapter-hit',
                        default=None,
                        help=('Set template position of second adapter hit.\n'
                              'Used for selecting data to fit second-pass\n'
                              'and rolling-circle taus. Optional flag.\n'
                              'Default value is smrtbellsize.'))
    parser.add_argument('--min-template-start',
                        default=0,
                        help=('Set minimum template-start position. Used for\n'
                              'rejecting alignments that did not start near \n'
                              'the beginning of the SMRTBell. Optional flag.\n'
                              'Default value is 0 bases.'))
    parser.add_argument('--max-template-start',
                        default=100,
                        help=('Set maximum template-start position. Used for\n'
                              'rejecting alignments that did not start near \n'
                              'the beginning of the SMRTBell. Optional flag.\n'
                              'Default value is 100 bases.'))
    parser.add_argument('--movie-length',
                        default=None,
                        help=('Set movie length (in minutes). Used for \n'
                              'flagging alignments that were movie-limited. \n'
                              'Optional flag. Will default to using all \n'
                              'alignments.'))
    parser.add_argument('--movie-limited-threshold',
                        default=0.5,
                        help=('Set threshold for determining when alignment\n'
                              'should be flagged as movie-limited. Value \n'
                              'sets duration, in minutes, before movie end.\n'
                              'Optional flag. Default value is 0.5 minutes.'))
    parser.add_argument('--base-file',
                        default=None,
                        help=('Provide path to bas.h5 file. To be used in \n'
                              'conjunction with legacy cmp.h5 format to \n'
                              'retrieve start-time and end-time of \n'
                              'alignments. Required for flagging alignments \n'
                              'that were movie-limited. Optional flag.\n'
                              'Will default to all alignments.'))
    parser.add_argument('--coarse-grain-binsize',
                        default=100,
                        help=('Choose binsize for coarse-graining survival\n'
                              'curve during termination rate calulation.\n'
                              'This does not affect plotted survival curve,\n'
                              'only the measurement of the termination rate.\n'
                              'Optional flag. Default value is 100 bases'))
    parser.add_argument('--subsample-to',
                        default=10000,
                        help=('Select number of alignments to subsample to.\n'
                              'Defaults to 10000'))
    parser.add_argument('--legacy',
                        action='store_true',
                        help='Use if running on cmp.h5 (legacy) format')
    args = parser.parse_args()
    return (args.alignments,
            args.base_file,
            args.output_directory,
            args.smrtbellsize,
            args.first_adapter_hit,
            args.second_adapter_hit,
            args.min_template_start,
            args.max_template_start,
            args.coarse_grain_binsize,
            args.movie_length,
            args.movie_limited_threshold,
            args.subsample_to,
            args.legacy)

def main():
    """Perform tau-analysis on an alignmentset."""
    # Parse inputs and determine alignmentset format
    (alignments_path,
     base_file_path,
     output_directory,
     smrtbellsize,
     first_adapter_hit,
     second_adapter_hit,
     template_min_start,
     template_max_start,
     coarse_grain_binsize,
     movie_length,
     movie_limited_threshold,
     subsample_to,
     is_legacy) = parse_args()

    if is_legacy:
        logging.info('Alignmentset expected to be cmp.h5 format.')
    else:
        logging.info('Alignmentset expected to be bam format.')

    # grab reference information and smrtbellsize
    alignment_set = ut.open_alignment_set(alignments_path, is_legacy)
    if ut.override_movie_limited_flag(alignment_set, is_legacy):
        movie_length = None # can't do movie-limited analysis

    if smrtbellsize is None: # try to get smrtbell size from reference name
        reference_information = ut.get_smrtbellsize(alignment_set)
    else:
        logging.info('SMRTbell size set by user. All templates are assumed\n'
                     'to have same length of ' + str(smrtbellsize) + ' '
                     'bases.')
        reference_information = alignment_set.referenceInfoTable[['ID',
                                                                  'FullName']]
        smrtbellsize = np.array([smrtbellsize for ref in reference_information],
                                dtype=int)
        reference_information = recfuncs.append_fields(reference_information.copy(),
                                                       'SMRTBellSize',
                                                       data=smrtbellsize,
                                                       asrecarray=True)
    # open basfile if provided and appropriate
    if not is_legacy and base_file_path is not None:
        raise RuntimeError('Bas.h5 path was provided, but alignment set '
                           'is bam format. Bas.h5 file is unnecessary.')
    elif not is_legacy:
        bas_reader = None
    elif is_legacy and (base_file_path and movie_length) is not None:
        bas_reader = ut.open_base_file(base_file_path)
    elif movie_length is None and base_file_path is not None:
        logging.info('Bas.h5 provided, but movie length was not. Analysis\n'
                     'will be done without tracking movie-limited '
                     'alignments.')
        bas_reader = None
    else:
        logging.info('Base file was not provided as input. Analysis will '
                     'be done without tracking movie-limited alignments, '
                     'even if movie length was inputted.')
        bas_reader = None

    ta.get_taus(alignment_set,
                bas_reader,
                reference_information,
                output_directory,
                movie_length,
                movie_limited_threshold,
                first_adapter_hit,
                second_adapter_hit,
                template_min_start,
                template_max_start,
                int(coarse_grain_binsize),
                int(subsample_to),
                is_legacy)


if __name__ == '__main__':
    main()
