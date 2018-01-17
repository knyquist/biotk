from pbcore.io import SubreadSet
import argparse
import itertools
import numpy as np

UINTMAX16 = 65536

def ranges(s):
    try:
        xmin, xmax = map(int, s.split(','))
        return xmin, xmax
    except:
        raise argparse.ArgumentTypeError('Coordinates must by x,y')

def parseArgs():
    parser = argparse.ArgumentParser(description=\
                                     'Create filtered subreadsets by x, y ' + \
                                     'coordinates')
    parser.add_argument('--mergeOutput',
                        help='If flagged, filter should produce one subreadset',
                        dest='one_output',
                        action='store_true')
    parser.add_argument('--subreadset',
                        help='Select subreadset for creating filtered ' + \
                             'subreadsets',
                        dest='sset')
    parser.add_argument('--xranges',
                        help='Select xranges for creating filtered ' + \
                             'alignmentsets xmin1,xmax1 xmin2, xmax2, ...',
                        dest='xranges',
                        type=ranges,
                        nargs='*')
    parser.add_argument('--yranges',
                        help='Select yranges for creating filtered ' + \
                             'alignmentsets ymin1,ymax1 ymin2,ymax2, ...',
                        dest='yranges',
                        type=ranges,
                        nargs='*')
    parser.add_argument('--outputDir',
                        help='Choose directory for output',
                        dest='output_dir')
    parser.add_argument('--name',
                        help='Choose a name for the filtered files',
                        dest='name')
    args = parser.parse_args()
    return args.one_output, args.sset, args.xranges, args.yranges, args.output_dir, args.name

def validateCoordinates(xranges, yranges):
    if len(xranges) == len(yranges):
        pass
    else:
        raise ValueError('xranges and yranges must have same length')

def openSubreadSet(subreadset):
    try:
        sset = SubreadSet(subreadset)
        return sset
    except:
        raise IOError('subreadset could not be opened')

def filterSubreadSet(sset, coordinate_ranges, one_output):
    filtered_ssets = []
    holenumbers = []
    for x_range, y_range in coordinate_ranges:
        xmin, xmax = x_range
        ymin, ymax = y_range
        x_vals = np.arange(xmin, xmax+1, 1)
        y_vals = np.arange(ymin, ymax+1, 1)
        xy = np.array(np.meshgrid(x_vals, y_vals)).T.reshape(-1, 2)
        if one_output:
            holenumbers.append([v[0]*UINTMAX16 + v[1] for v in xy])
        else:
            holenumbers = np.array([v[0]*UINTMAX16 + v[1] for v in xy])
            tmp_sset = openSubreadSet(sset)
            tmp_sset.filters.addRequirement(zm=[('=', holenumbers)])
            filtered_ssets.append(tmp_sset)
    if one_output:
        holenumbers = itertools.chain(*holenumbers)
        holenumbers = np.array(list(holenumbers))
        tmp_sset = openSubreadSet(sset)
        tmp_sset.filters.addRequirement(zm=[('=', holenumbers)])
        filtered_ssets.append(tmp_sset)
    return filtered_ssets

def writeFilteredSubreadSets(filtered_subreadsets,
                             coordinate_ranges,
                             one_output,
                             output_dir, name):
    if one_output: # write one file, named appropriately
        for sset in filtered_subreadsets:    
            cx = 'x'
            cy = 'y'
            for x, y in coordinate_ranges:
                cx += str(x[0]) + '-' + str(x[1]) + '_'
                cy += str(y[0]) + '-' + str(y[1]) + '_'
            sset.newUuid()
            sset.write(output_dir + \
                       '/filtered_' + \
                       name + '_' + \
                       cx + cy[0:-1] + \
                       '.subreadset.xml')
    else: # write one output for each window
        for index, sset in enumerate(filtered_subreadsets):
            x_range = coordinate_ranges[index][0]
            y_range = coordinate_ranges[index][1]
            sset.newUuid()
            sset.write(output_dir + \
                       '/filtered_' + \
                       name + \
                       '_x' + \
                       str(x_range[0]) + '-' + str(x_range[1]) + \
                       '_y' + \
                       str(y_range[0]) + '-' + str(y_range[1]) + \
                       '.subreadset.xml')

def main():
    one_output, sset, xranges, yranges, output_dir, name = parseArgs()
    validateCoordinates(xranges, yranges)
    filtered_subread_sets = filterSubreadSet(sset, zip(xranges, yranges), one_output)
    writeFilteredSubreadSets(filtered_subread_sets, zip(xranges, yranges), one_output,
                             output_dir,
                             name)

if __name__ == '__main__':
    main()
