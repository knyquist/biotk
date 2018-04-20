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
    subparsers = parser.add_subparsers(help='Choose mode of operation',
                                       dest='command')
    summarize_parser = subparsers.add_parser('summarize')
    full_parser = subparsers.add_parser('full-mode')
    
    summarize_parser.add_argument('alignmentset',
                                  help='Path to AlignmentSet')
    summarize_parser.add_argument('output',
                                  help='Path to output CSV')
    summarize_parser.add_argument('--by-zmw',
                                  help=('Combine subreads from a particular ZMW \n'
                                        'when displaying readlength and template-span \n'
                                        'summary information. Flag does nothing when \n'
                                        'program is operated in full-mode.'),
                                  action="store_true")
    full_parser.add_argument('alignmentset',
                             help='Path to AlignmentSet')
    full_parser.add_argument('output',
                             help='Path to output CSV')    
    args = parser.parse_args()
    return args

def refIdToName(referenceInfoTable):
    id_to_name = {}
    for refInfo in referenceInfoTable:
        id_to_name[refInfo['ID']] = refInfo['Name']
    return id_to_name

def saveFullCsv(df, output_path):
    """
    Save full summary to CSV
    """
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

def saveSummaryCsv(df, output_path, by_zmw):
    """
    Save by-contig summary to CSV
    """
    df['readlength'] = df['aEnd'] - df['aStart']
    df['templatespan'] = df['tEnd'] - df['tStart']
    df.drop(['qId',
             'qStart',
             'qEnd',
             'tStart',
             'tEnd',
             'aEnd',
             'aStart',
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
    if not by_zmw:
        mean = df.groupby(by=['refName']).mean()
        mean['refName'] = mean.index
        median = df.groupby(by=['refName']).median()
        median['refName'] = median.index
        count = df.groupby(by=['refName']).count()
        count['nReads'] = count['accuracy'] # choose a column, all same
        sdf = pd.merge(mean,
                       median,
                       on='refName',
                       suffixes=['_mean', '_median'])
        sdf['nReads'] = pd.Series(count['nReads'].values)
        sdf['nZmws'] = pd.Series(df.groupby(by=['refName'])['holeNumber'].nunique().values)

    else: # group metrics by ZMW as well
        acc_mean = df.groupby(by=['refName'])['accuracy'].mean()
        acc_median = df.groupby(by=['refName'])['accuracy'].median()
        unique_count = df.groupby(by=['refName'])['holeNumber'].nunique()
        count = df.groupby(by=['refName'])['holeNumber'].count()
        sdf = pd.DataFrame({'refName': acc_mean.index,
                            'nReads': count,
                            'nZmws': unique_count,
                            'accuracy_mean': acc_mean,
                            'accuracy_median': acc_median})

        for length in ['readlength', 'templatespan']:
            lengths = df.groupby(by=['refName', 'holeNumber'])[length].sum()
            mean_length = lengths.groupby(by=['refName']).mean()
            median_length = lengths.groupby(by=['refName']).median()
            sdf[length + '_mean'] = mean_length
            sdf[length + '_median'] = median_length

    sdf = sdf[['refName', 'nReads', 'nZmws',
               'accuracy_mean', 'accuracy_median',
               'templatespan_mean', 'templatespan_median',
               'readlength_mean', 'readlength_median']]
    sdf.to_csv(output_path, index=False)
    log.info('Done.')

def generatePbiCsv(aset_path, output_path,
                   program_mode,
                   by_zmw):
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
    if program_mode == 'full-mode':
        log.info('Saving full-mode CSV')
        saveFullCsv(df, output_path)
    elif program_mode == 'summarize':
        log.info('Saving summarize CSV')
        saveSummaryCsv(df, output_path, by_zmw)

def main():
    log.setLevel(logging.INFO)
    args = parseArgs()
    log.info('Parsed command-line arguments...')
    program_mode = args.command

    # Based on command-line arguments, decide
    # whether the results should presented on
    # a per-subread bases or on a per-zmw basis.
    by_zmw = False
    if args.by_zmw:
        by_zmw = True

    generatePbiCsv(args.alignmentset, args.output,
                   program_mode,
                   by_zmw)


if __name__ == '__main__':
    main()