from pbcore.io import AlignmentSet
import pandas as pd
import numpy as np
import argparse
import logging

logging.basicConfig()
log = logging.getLogger(__name__)

def parseArgs():
    parser = argparse.ArgumentParser(description=('Generate AlignmentSet'
                                                  'summary metrics '
                                                  'from .pbi file.'))
    parser.add_argument('alignmentset',
                        help='Path to AlignmentSet')
    parser.add_argument('output',
                        help='Path to output CSV')
    args = parser.parse_args()
    return (args.alignmentset,
            args.output)

def refIdToName(referenceInfoTable):
    id_to_name = {}
    for refInfo in referenceInfoTable:
        id_to_name[refInfo['ID']] = refInfo['Name']
    return id_to_name

def generatePbiCsv(aset_path, output_path):
    """
    Given alignmentset, generate CSV summary of
    pbi index and save to output_path
    """
    aset = AlignmentSet(aset_path)
    log.info('Opened AlignmentSet for processing...')
    df = pd.DataFrame.from_records(aset.index)
    df['refName'] = df['tId'].map(refIdToName(aset.referenceInfoTable))
    df['accuracy'] = np.round(1. - np.divide(df['nDel'] + df['nMM'] + df['nIns'],
                                          df['tEnd'] - df['tStart'],
                                          dtype=float),
                              decimals=2)
    df.drop(['qId',
             'qStart',
             'qEnd',
             'tId',
             'nM',
             'nMM',
             'nIns',
             'nDel',
             'readQual',
             'contextFlag',
             'virtualFileOffset',
             'isReverseStrand',
             'mapQV'],
             axis=1,
             inplace=True)
    log.info('Saving info to CSV...')
    df.to_csv(output_path, index=False)
    log.info('Done.')

def main():
    log.setLevel(logging.INFO)
    (aset_path,
     output_path) = parseArgs()
    log.info('Parsed command line arguments...')
    generatePbiCsv(aset_path, output_path)


if __name__ == '__main__':
    main()