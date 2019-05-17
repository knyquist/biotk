"""
Generate SMRTLink URL

Usage:
  view_trace.py [-hv]
        [--server_host server_host] [--server_port server_port] [--app_path app_path]
        [-t trc] [-s sset] [-a aset] [-z zmw] [-o browser]

Options:
  -h --help                             This message
  -v --version                          Print version and exit
     --server_host=<server_host>        Hostname of smrtserver [default: http://smrtserver]
     --server_port=<server_port>        Port of smrtserver [default: 8080]
     --app_path=<app_path>              Path to traceviewer app on server [default: smrtserver/store/apps/traceviewer]
  -t --trc=<trc>                        Path to trc.h5 file
  -s --sset=<sset>                      Path to subreadset.xml file
  -a --aset=<aset>                      Path to alignmentset.xml file
  -z --zmw=<zmw>                        ZMW holenumber
  -o --browser=<browser>                Auto open in browser [default: False]
"""


import os
import sys
from docopt import docopt
import webbrowser
import logging, coloredlogs
log = logging.getLogger('TRACEVIEWER')
log.addHandler(logging.StreamHandler())
log.setLevel(logging.DEBUG)
coloredlogs.install(level='DEBUG', logger=log)


def validate_args(args):
    """
    The trc file and zmw are required
    """
    cnt = 0
    if os.path.isfile(args['--trc']):
        cnt += 1
        log.info('Path to {trc} exists'.format(trc=args['--trc']))
    if args['--zmw'] != 'None':
        cnt += 1
        log.info('Validated --zmw is not None')
    if cnt < 2:
        log.error('Trace file and ZMW must be provided, exiting')
        sys.exit(1)

    if os.path.isfile(args['--sset']):
        log.info('Path to {sset} exists'.format(sset=args['--sset']))
    else:
        log.info('Sset path: {sset} not valid or not filled out'.format(sset=args['--sset']))
        args['--sset'] = None

    if os.path.isfile(args['--aset']):
        log.info('Path to {aset} exists'.format(aset=args['--aset']))
    else:
        log.info('Aset path: {aset} not valid or not filled out'.format(aset=args['--aset']))
        args['--aset'] = None

    return args

def construct_traceview_url(args):
    """
    Construct url using artifacts
    """
    url = "{host}:{port}/{app}/traceviewer.html?trace={trc}&alignmentset={aln}&subreads={sset}&holen={zmw}"
    url = url.format(host=args['--server_host'],
                     port=args['--server_port'],
                     app=args['--app_path'],
                     trc=args['--trc'],
                     aln=args['--aset'],
                     sset=args['--sset'],
                     zmw=args['--zmw'])
    # log.info("Opening up {url} in the browser".format(url))
    return url

if __name__ == '__main__':
    args = docopt(__doc__, version='View Trace 0.1')
    print(args)
    args = validate_args(args)
    url = construct_traceview_url(args)
    webbrowser.open_new_tab(url)