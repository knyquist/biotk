from pbcore.io import (AlignmentSet,
                       SubreadSet)
import argparse
import random
import numpy as np
import plotly
from plotly import tools
from plotly.graph_objs import *
from plotly.tools import FigureFactory as FF
from plotly.offline import download_plotlyjs, plot

UINTMAX16 = 65536

def parseArgs():
    parser = argparse.ArgumentParser(description='Plot heatmap of subreadsets.')
    parser.add_argument('datasets',
                        help='List the datasets to combine for heatmap.',
                        nargs='*')
    parser.add_argument('-d',
                        '--dataSetType',
                        help=('Specify type of dataset, '
                              'AlignmentSet or SubreadSet.'),
                        dest='dtype')
    parser.add_argument('-s',
                        '--subSampleTo',
                        help=('Specify number of ZMWs to subsample each '
                              'subreadset to. Default is 1000.'),
                        default=1000,
                        dest='subsampleto'),
    parser.add_argument('-t',
                        '--title',
                        help='Set title of heatmap. Default is None.',
                        default=None,
                        dest='title')
    parser.add_argument('-o',
                        '--output',
                        help='Specify location of output.',
                        dest='output')
    args = parser.parse_args()
    return args.datasets, args.dtype, args.subsampleto, args.title, args.output

def main():
    datasets, dtype, subsampleto, title, output = parseArgs()
    d = []
    for dset in datasets:
        if dtype == 'AlignmentSet':
            f = AlignmentSet(dset)
        elif dtype == 'SubreadSet':
            f = SubreadSet(dset)
        else:
            raise ValueError('invalid dataSetType')
        x = f.index['holeNumber'] / UINTMAX16
        y = f.index['holeNumber'] - x * UINTMAX16
        if len(f) > subsampleto:
            x, y = zip(*random.sample(zip(x, y), subsampleto))
        h = Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                size=5,
                opacity=0.2
            ),
            showlegend=False
        )
        d.append(h)
    layout = Layout(
        title=title,
        height=600,
        width=600,
        xaxis=dict(
            title='X',
            range=[0, 1500]
        ),
        yaxis=dict(
            title='Y',
            range=[0, 1500]
        )
    )
    fig = Figure(data=d, layout=layout)
    plot(fig, show_link=False, auto_open=False, filename=output)

if __name__ == '__main__':
    main()
