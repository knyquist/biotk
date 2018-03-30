import argparse
import pandas
import csv
import numpy as np
import random
random.seed(a=666)
import logging
import os
import plotly
from plotly.graph_objs import *
from plotly.offline import plot


logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def parseArgs():
    """
    CLI inputs
    """
    log.info('Parsing command-line arguments...')
    parser = argparse.ArgumentParser(description='Explore enzyme-screening variability historically')
    subparsers = parser.add_subparsers(help='Choose mode of operation')
    plot_parser = subparsers.add_parser('only-plots')
    full_parser = subparsers.add_parser('daisy-chain')

    full_parser.add_argument('stats',
                             help='Path to summary_stats.csv produced by sequin')
    full_parser.add_argument('output_old',
                             help='Path to previous output')
    full_parser.add_argument('output_new',
                             help='Path to new output')
    full_parser.add_argument('enzymes',
                             help='Name of enzyme(s) to retrieve data for, space-separated',
                             nargs='*')
    
    plot_parser.add_argument('output_old',
                             help='Path to previous output')
    plot_parser.add_argument('output_path',
                             help='Path for plots')
    
    args = parser.parse_args()
    return args
    #return args.stats, args.output_old, args.output_new, args.enzymes

def get_fieldnames(arr):
    """
    Get fieldnames from recarray
    """
    common_column_names = ['esPlate', 'chip', 'ESID', 'plateId', 'plateName', 'limsExptName', 'runCode',
                           'analogType', 'analogConcentration', 'movieTime',
                           'platform', 'chipLotId', 'chipType', 'esComment', 'laserPower']
    fn1 = []
    for val in arr.dtype.names:
        n = val.split('test')
        if len(n) == 1:
            fn1.append(n)
        else:
            if len(n[0]) > 0:
                n[0] = n[0][0:-1]
            fn1.append(''.join(n))    
    fn2 = []
    for val in arr.dtype.names:
        n = val.split('control')
        if len(n) == 1:
            fn2.append(n)
        else:
            if len(n[0]) > 0:
                n[0] = n[0][0:-1]
            fn2.append(''.join(n))
            
    ans = np.intersect1d(fn1, fn2)
    r = []
    for col in ans:
        if type(col) is list:
            r.append(col[0])
        else:
            r.append(col)
            
    f = []
    for col in r:
        s = col.split('.')
        if ('delta' in s) | ('ratio' in s):
            continue
        else:
            f.append(col)
    
    return f

def retrieve_new_data(stats, enzymes):
    """
    Retrieve summary stats about chosen enzymes
    """
    log.info('Retrieving new summary stats...')
    summary_stats = pandas.read_csv(stats, dtype=str)
    data = []
    places = ['test', 'control']
    for enzyme in enzymes:
        for place in places: # look at test and controls
            df = summary_stats[summary_stats[place + 'Enzyme'] == enzyme]
            df = df.drop_duplicates(subset=[place+'Enzyme', place+'Template'])
            if df.values.any():
                df = df.to_records()
                columns = df.dtype.names
                fieldnames = get_fieldnames(df)
                for row in df:
                    data_row = {}
                    for column in columns:
                        if column in fieldnames:
                            data_row[column] = row[column]
                        else:
                            new_column = column
                            m = column.split(place)
                            if len(m) > 1:
                                if len(m[0]) > 0:
                                    m[0] = m[0][0:-1]
                                    new_column = ''.join(m)
                                else:
                                    new_column = ''.join(m)
                                data_row[new_column] = row[column]
                    data.append(data_row)
    
    data = pandas.DataFrame.from_dict(data)
    data['plateName'] = pandas.Series(data['esPlate'].str.cat(data['chip'].values.astype(int).astype(str), sep='-'), index=data.index)
    data['plateName'] = data['plateName'].astype('category')
    return data
    
def read_previous_data(output):
    """
    Read the data in the previous CSV
    """
    log.info('Reading previous data...')
    try:
        old_data = pandas.read_csv(output, dtype=str)
        old_data['plateName'] = pandas.Series(old_data['esPlate'].str.cat(old_data['chip'].values.astype(int).astype(str), sep='-'), index=old_data.index)
        old_data['plateName'] = old_data['plateName'].astype('category')
        return old_data
    except ValueError:
        log.info('    {i} is empty'.format(i=output))
        return None
    except IOError:
        log.info('    {i} does not exist'.format(i=output))
        return None
        
def merge_data(old, new):
    """
    Merge the previous CSV with the new one,
    using a full outer join
    """
    log.info('Merging together datasets...')
    if old is not None:
        merged_df = pandas.merge(old, new, how='outer')
        #merged_df.sort_values(by=['esPlate', 'chip', 'ESID'], inplace=True)
        merged_df['plateName'] = merged_df['plateName'].astype('category')
        merged_df['plateId'] = merged_df['plateName'].cat.rename_categories(np.arange(1, len(np.unique(merged_df['plateName']))+1, 1))
        return merged_df
    else:
        return new

def write_output(output, df):
    """
    Write the merged df to output csv
    """
    log.info('Writing output...')
    df.to_csv(output, index=False)
    
def scatter_object(df, enzyme, y, symbol):
    """
    Make scatter plot object for particular enzyme
    """
    edf = df[df['Enzyme'] == enzyme]
    d = Scatter(
        x = edf['plateId'],
        y = edf[y],
        name = enzyme,
        mode = 'markers',
        marker = dict(
            symbol=symbol,
            size=20
        )
    )
    return d
    
def make_plots(merged_data, path):
    """
    Make performance plots
    """
    log.info('Generating plots...')
    if not merged_data.empty:
        features = ['median.snr.C', 'median.snr.T',
                    'median.pkmid.C', 'median.pkmid.T',
                    'median.baseline.sigma.C', 'median.baseline.sigma.T',
                    'median.baseline.level.C', 'median.baseline.level.T',
                    'median.acc',
                    'median.ts',
                    'polrate',
                    'tmean90.ipd.A', 'tmean90.ipd.C', 'tmean90.ipd.T', 'tmean90.ipd.G',
                    'tmean95.pw.A', 'tmean95.pw.C', 'tmean95.pw.T', 'tmean95.pw.G',
                    'del.DarkA', 'del.DarkC', 'del.DarkT', 'del.DarkG',
                    'laserPower']

        enzymes = np.unique(merged_data['Enzyme'].values)
        symbols = random.sample(np.arange(100, 144, 1), len(enzymes))
        ticktext = {}
        for index, entry in enumerate(merged_data['plateId'].values):
            if entry not in ticktext:
                ticktext[entry] = merged_data['plateName'].values[index]
        tickvals = np.unique(merged_data['plateId'].values)
        for feature in features:
            if feature in merged_data.columns:
                scats = []
                for index, enzyme in enumerate(enzymes):
                    symbol = symbols[index]
                    scat = scatter_object(merged_data, enzyme, feature, symbol)
                    scats.append(scat)
                layout=Layout(
                    showlegend=True,
                    font = dict(
                        size = 20
                    ),
                    xaxis = dict(
                        tickangle=45,
                        tickvals=tickvals,
                        ticktext=[ticktext[val] for val in tickvals]
                    ),
                    yaxis = dict(
                        title = feature
                    )
                )
                fig = Figure(data=scats, layout=layout)
                plot(fig, filename=path+feature+'.html', show_link=False, auto_open=False)
    
def main():
    """
    Extract metrics for particular enzymes and catenate to output file.
    """
    args = parseArgs()
    if hasattr(args, 'stats'): # daisy-chain mode
        log.info('Running in daisy-chain mode...')
        stats, output_old, output_new, enzymes = (args.stats,
                                                  args.output_old,
                                                  args.output_new,
                                                  args.enzymes)
        old_data = read_previous_data(output_old)
        new_data = retrieve_new_data(stats, enzymes)
        merged_data = merge_data(old_data, new_data)
        write_output(output_new, merged_data)
        path = os.path.split(output_new)[0] + os.path.sep
        make_plots(merged_data, path)
        log.info('Done.')
    else: # plot only mode
        log.info('Running in plot-only mode...')
        output_old, path = args.output_old, args.output_path
        old_data = read_previous_data(output_old)
        make_plots(old_data, path)
        log.info('Done.')
    
if __name__ == '__main__':
    main()
