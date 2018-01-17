#!/pbi/dept/itg/modules/biochemistry-toolkit/0.4/bin/python

import numpy as np
import random
import string
import operator
from pbcore.io import (AlignmentSet,
                       SubreadSet)
from pricompare import FastMetrics as fm
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime

class BurstMetrics:
    """
    Class for retrieving burst metrics. Two flavors
        Alignment based
        HMM Classifier based
        
    If required information not available, return None.
    """
    def __init__(self, subread_set_path,
                       alignment_set_path=None,
                       reference_set_path=None,
                       subsampleto=None,
                       false_negative_burst_length=10):
        if alignment_set_path:
            self.alignment_set_path = alignment_set_path
            self.alignment_set = AlignmentSet(alignment_set_path,
                                              referenceFastaFname=reference_set_path)
        self.subread_set_path = subread_set_path
        self.subread_set = SubreadSet(subread_set_path)
        self.hqregions_constructor = NoSplitSubreads(self.subread_set_path)
        self.hqregion_set = self.hqregions_constructor.hqregion_reads
        self.subsampleto = subsampleto
        self.false_negative_burst_length = false_negative_burst_length

    def resize_array(self, arr, index, increase_by):
        """
        Resize NumPy array if necessary
        """
        if index >= len(arr): # extend array if needed
                new_size = tuple(map(operator.add, arr.shape, (increase_by, )))
                arr = np.resize(arr, new_size)
        return arr
        
    def retrieve_classifier_bursts(self, dset):
        """
        Retrieve information about the bursts detected by the classifier.
        Returns a recarray with the following columns:
            zmw
            queryStart
            queryEnd
            burstStart
            burstLength
            burstStartTime
            burstEndTime
            fractionC
            fractionA
            fractionT
            fractionG
        """
        if 'pe' not in [t[0] for t in dset[0].peer.tags]: # check for burst classification
            return None 
        if self.subsampleto is not None: # default behavior does not subsample
            if len(dset) > self.subsampleto: # subsample subreadset randomly
                subsample_zmws = np.unique(random.sample(dset.index['holeNumber'],
                                                         self.subsampleto))
                dset.filters.addRequirement(zm=[('=', subsample_zmws)])
        bursts = np.zeros((len(dset), ), dtype=[('zmw', int),
                                                ('qStart', int),
                                                ('qEnd', int),
                                                ('burstStart', int),
                                                ('burstLength', int),
                                                ('burstStartTime', int),
                                                ('burstEndTime', int),
                                                ('fractionC', float),
                                                ('fractionA', float),
                                                ('fractionT', float),
                                                ('fractionG', float)])
        bases = ['a', 'c', 'g', 't']
        burst_count = 0
        for read in dset:
            bursty_indices = np.flatnonzero(np.array(read.peer.get_tag('pe')) == 2)
            bursty_gaps = np.diff(bursty_indices)
            bursty_breaks = np.flatnonzero(bursty_gaps > self.false_negative_burst_length)
            if bursty_indices.any():
                bursts['zmw'][burst_count] = read.holeNumber
                bursts['qStart'][burst_count] = read.qStart
                bursts['qEnd'][burst_count] = read.qEnd
                start_frames = read.peer.get_tag('sf')
                bursts['burstStart'][burst_count] = bursty_indices[0]
                if bursty_breaks.any():
                    for bursty_break in bursty_breaks:
                        bursts = self.resize_array(bursts, burst_count, self.subsampleto)
                        index = bursty_indices[bursty_break]
                        bursts['burstLength'][burst_count] = index - bursts['burstStart'][burst_count] + 1
                        bursts['burstStartTime'][burst_count] = start_frames[bursts['burstStart']
                                                                                   [burst_count]]
                        bursts['burstEndTime'][burst_count] = start_frames[(bursts['burstStart'][burst_count] +
                                                                            bursts['burstLength'][burst_count])]
                        bs = bursts['burstStart'][burst_count]
                        be = bursts['burstStart'][burst_count] + bursts['burstLength'][burst_count]
                        burstcalls = list(read.peer.get_tag('pc')[bs:be])
                        for base in bases:
                            f1 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == base)),
                                           len(burstcalls),
                                           dtype=float) # include rejected bases (pulses)
                            f2 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == string.upper(base))),
                                           len(burstcalls),
                                           dtype=float) # include basecalls
                            bursts['fraction' + string.upper(base)][burst_count] = f1 + f2
                        burst_count += 1
                        bursts = self.resize_array(bursts, burst_count, self.subsampleto)
                        bursts['zmw'][burst_count] = read.holeNumber
                        bursts['qStart'][burst_count] = read.qStart
                        bursts['qEnd'][burst_count] = read.qEnd
                        next_index = bursty_indices[bursty_break + 1]
                        bursts['burstStart'][burst_count] = next_index
                bursts['burstLength'][burst_count] = bursty_indices[-1] - bursts['burstStart'][burst_count] + 1
                bursts['burstStartTime'][burst_count] = start_frames[bursts['burstStart'][burst_count]]
                bursts['burstEndTime'][burst_count] = start_frames[(bursts['burstStart'][burst_count] +
                                                                    bursts['burstLength'][burst_count])]
                bs = bursts['burstStart'][burst_count]
                be = bursts['burstStart'][burst_count] + bursts['burstLength'][burst_count]
                burstcalls = list(read.peer.get_tag('pc')[bs:be])
                for base in bases:
                    f1 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == base)),
                                   len(burstcalls),
                                   dtype=float) # include rejected bases (pulses)
                    f2 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == string.upper(base))),
                                   len(burstcalls),
                                   dtype=float) # include basecalls
                    bursts['fraction' + string.upper(base)][burst_count] = f1 + f2
                burst_count += 1
                bursts = self.resize_array(bursts, burst_count, self.subsampleto)
            
        bursts = bursts[bursts['zmw'] != 0]
        return bursts

    def retrieve_alignment_bursts(self, aset):
        """
        Use insertions from the alignmentset to tabulate bursts.
        Returns a recarray with the following columns:
            zmw
            queryStart
            queryEnd
            burstStart
            burstLength
            burstStartTime
            burstEndTime
            fractionC
            fractionA
            fractionT
            fractionG
        """
        if len(aset) > self.subsampleto: # subsample subreadset randomly
            subsample_zmws = np.unique(random.sample(aset.index['holeNumber'],
                                                     self.subsampleto))
            aset.filters.addRequirement(zm=[('=', subsample_zmws)])
        bursts = np.zeros((len(aset), ), dtype=[('zmw', int),
                                                ('qStart', int),
                                                ('qEnd', int),
                                                ('burstStart', int),
                                                ('burstLength', int),
                                                ('burstStartTime', int),
                                                ('burstEndTime', int),
                                                ('fractionC', float),
                                                ('fractionA', float),
                                                ('fractionT', float),
                                                ('fractionG', float)])
        bases = ['a', 'c', 'g', 't']
        burst_count = 0
        for read in aset:
            ref_sequence = np.array(list(read.reference()), dtype='S')
            read_sequence = np.array(list(read.read()), dtype='S')
            # remove deletions to keep base-space
            ref_sequence = ref_sequence[read_sequence != '-']
            read_sequence = read_sequence[read_sequence != '-']
            bursty_indices = np.flatnonzero(ref_sequence == '-')
            b2p = fm.base2pls(fm.s2npl(read.peer.get_tag('pc')))
            bursty_indices = b2p[bursty_indices] # in pulse indices
            bursty_gaps = np.diff(bursty_indices)
            bursty_breaks = np.flatnonzero(bursty_gaps > self.false_negative_burst_length)
            if bursty_indices.any():
                bursts['zmw'][burst_count] = read.holeNumber
                bursts['qStart'][burst_count] = read.qStart
                bursts['qEnd'][burst_count] = read.qEnd
                if 'sf' in [t[0] for t in read.peer.tags]:    
                    start_frames = read.peer.get_tag('sf')
                else:
                    start_frames = None # start frame info not available
                bursts['burstStart'][burst_count] = bursty_indices[0]
                if bursty_breaks.any():
                    for bursty_break in bursty_breaks:
                        bursts = self.resize_array(bursts, burst_count, self.subsampleto)
                        index = bursty_indices[bursty_break]
                        bursts['burstLength'][burst_count] = (index -
                                                              bursts['burstStart'][
                                                                  burst_count] + 1)
                        if start_frames:
                            bursts['burstStartTime'][burst_count] = start_frames[
                                                                        bursts['burstStart'][
                                                                            burst_count]]
                            bursts['burstEndTime'][burst_count] = start_frames[(bursts['burstStart'][burst_count] +
                                                                                bursts['burstLength'][burst_count])]
                        bs = bursts['burstStart'][burst_count]
                        be = bursts['burstStart'][burst_count] + bursts['burstLength'][burst_count]
                        burstcalls = list(read.read()[bs:be])
                        for base in bases:
                            f = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == base)),
                                          len(burstcalls),
                                          dtype=float)
                            bursts['fraction' + string.upper(base)][burst_count] = f
                        burst_count += 1
                        bursts = self.resize_array(bursts, burst_count, self.subsampleto)
                        bursts['zmw'][burst_count] = read.holeNumber
                        bursts['qStart'][burst_count] = read.qStart
                        bursts['qEnd'][burst_count] = read.qEnd
                        next_index = bursty_indices[bursty_break + 1]
                        bursts['burstStart'][burst_count] = next_index
                bursts['burstLength'][burst_count] = bursty_indices[-1] - bursts['burstStart'][burst_count] + 1
                if start_frames:
                    bursts['burstStartTime'][burst_count] = start_frames[bursts['burstStart'][burst_count]]
                    bursts['burstEndTime'][burst_count] = start_frames[(bursts['burstStart'][burst_count] +
                                                                        bursts['burstLength'][burst_count])]
                bs = bursts['burstStart'][burst_count]
                be = bursts['burstStart'][burst_count] + bursts['burstLength'][burst_count]
                burstcalls = list(read.peer.get_tag('pc')[bs:be])
                for base in bases:
                    f1 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == base)),
                                   len(burstcalls),
                                   dtype=float) # include rejected bases (pulses)
                    f2 = np.divide(len(np.flatnonzero(np.array(burstcalls, 'S') == string.upper(base))),
                                   len(burstcalls),
                                   dtype=float) # include basecalls
                    bursts['fraction' + string.upper(base)][burst_count] = f1+f2
                burst_count += 1
                bursts = self.resize_array(bursts, burst_count, self.subsampleto)
            
        bursts = bursts[bursts['zmw'] != 0]
        return bursts
        
class NoSplitSubreads():
    """
    Class for generating and opening a ZmwRead Dataset from a
    SubreadSet using bam2bam.
    """
    
    def __init__(self, sset_path,
                       output_folder='/tmp/',
                       no_split_mode='--hqregion'):
        self.sset_path = sset_path
        (self.hqregion_reads,
         self.bam2bam_stdout,
         self.bam2bam_stderr) = self.build_zmw_read()
    
    def build_zmw_read(self):
        """
        Use bam2bam to build ZMW read. Return DataSet object.
        """
        cmd = ['bam2bam',
               self.sset_path]
        output_path = (output_folder +
                       str(datetime.now()))
        cmd = cmd + ['-o', output_path]
        options = ['-j', '4', # number of processors
                   '-b', '4', # number of processors
                   no_split_mode]
        cmd = cmd + options
        
        stdout = None
        stderr = None
        subprocess.check_call(cmd,
                              stdout=stdout,
                              stderr=stderr)
        return (SubreadSet(output_path + '.subreadset.xml'),
                stdout,
                stderr)
