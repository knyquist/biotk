from pbcore.io import SubreadSet
from pricompare import FastMetrics as fm
from biotk.libs.BurstMetrics import NoSplitSubreads
import argparse
import logging
import numpy as np
import csv

logging.basicConfig()
log = logging.getLogger(__name__)

def parseArgs():
    parser = argparse.ArgumentParser(description=('Create '
                                                  'input file '
                                                  'for bammend.'))
    parser.add_argument('subreadset',
                        help='Path to subreadset')
    parser.add_argument('output',
                        help='Path to output')
    parser.add_argument('-s',
                        '--reject-start',
                        default=0,
                        help=('Choose pulse-rejection start time '
                              '(minutes)'),
                        required=True)
    parser.add_argument('-e',
                        '--reject-end',
                        default=45*60*80,
                        help=('Choose pulse-rejection end time '
                              '(minutes)'),
                        required=True)
    args = parser.parse_args()
    return (args.subreadset,
            args.output,
            args.reject_start,
            args.reject_end)

def main():
    log.setLevel(logging.INFO)
    (subreadset,
     output,
     reject_start,
     reject_end) = parseArgs()
    sset = SubreadSet(subreadset)
    with open(output, 'wb') as csvfile:
        fieldnames = ['holeNumber',
                      'annotationStartIndex',
                      'annotationEndIndex']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        log.info('Iterating through subreads and writing rejection indices.')
        for subread in sset:
            framerate = subread.readGroupInfo['FrameRate']
            rs = float(reject_start) * 60. * framerate # reject start in frames
            re = float(reject_end) * 60. * framerate # reject end in frames
            sf = np.array(subread.peer.get_tag('sf'))
            pc = fm.s2npl(subread.peer.get_tag('pc'))
            b2p = fm.base2pls(pc)
            bsf = sf[b2p] # start frames of basecalls
            reject_indices = np.flatnonzero((bsf >= rs) &
                                            (bsf <= re))
            if not reject_indices.any():
                log.info(('Subread from ZMW {i} has no bases '
                          'to reject.'.format(i=subread.holeNumber)))
                continue
            annotation_start_index = (reject_indices[0] + subread.aStart)
            annotation_end_index = (reject_indices[-1] + subread.aStart)
            row = {'holeNumber': subread.holeNumber,
                   'annotationStartIndex': annotation_start_index,
                   'annotationEndIndex': annotation_end_index}
            writer.writerow(row)
    log.info('Done writing rejection indices to {i}'.format(i=output))

if __name__ == '__main__':
    main()