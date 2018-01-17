#!/usr/bin/env python

import numpy as np
import random
import string
import operator
from pbcore.io import (AlignmentSet,
                       SubreadSet,
                       IndexedBamReader)
from pricompare import FastMetrics as fm
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime

class PpaBurstMetrics:
    """
    Class for retrieving burst metrics. Two flavors
        Alignment based
        HMM Classifier based
        
    If required information not available, return None.
    """
    def __init__(self, subread_set_path,
                       subsampleto=None,
                       false_negative_burst_length=10):
        if alignment_set_path:
            self.alignment_set_path = alignment_set_path
            self.alignment_set = AlignmentSet(alignment_set_path,
                                              referenceFastaFname=reference_set_path)
        self.subsampleto = subsampleto
        self.subread_set_path = subread_set_path
        self.false_negative_burst_length = false_negative_burst_length
        self.subread_set = SubreadSet(subread_set_path)
        self.scraps = IndexedBamReader(self.subread_set.externalResources[0].scraps)
        self.subread_ppa_bursts = self.retrieve_classifier_bursts(self.subread_set)
        self.scraps_ppa_bursts = self.retrieve_classifier_bursts(self.scraps)
        self.ppa_bursts = np.hstack((self.subread_ppa_bursts,
                                     self.scraps_ppa_bursts))

    def resize_array(self, arr, index, increase_by):
        """
        Resize NumPy array if necessary
        """
        if index >= len(arr): # extend array if needed
            new_size = tuple(map(operator.add, arr.shape, (increase_by, )))
            arr = np.resize(arr, new_size)
        return arr

    def subsample_zmws(self, dset):
        """
        Subsample Zmws for measurement
        """
        zmws = np.unique(dset.index['holeNumber'])
        if self.subsampleto is not None:
            if len(dset) > self.subsampleto:
                zmws = np.unique(random.sample(dset.index['holeNumber'],
                                               self.subsampleto))
        return zmws

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
            previousBasecall
            previousBaseIndex
            fractionC
            fractionA
            fractionT
            fractionG
        """
        if 'pe' not in [t[0] for t in dset[0].peer.tags]: # check for burst classification
            return None

        bursts = np.zeros((len(dset), ), dtype=[('zmw', int),
                                                ('qStart', int),
                                                ('qEnd', int),
                                                ('burstStart', int),
                                                ('burstLength', int),
                                                ('burstStartTime', int),
                                                ('burstEndTime', int),
                                                ('previousBasecall', 'S1'),
                                                ('previousBaseIndex', int),
                                                ('fractionC', float),
                                                ('fractionA', float),
                                                ('fractionT', float),
                                                ('fractionG', float)])
        bases = ['a', 'c', 'g', 't']
        burst_count = 0

        zmws = self.subsample_zmws(dset)
        read_indices = np.flatnonzero(np.in1d(dset.index['holeNumber'], zmws))
        for index in read_indices:
            read = dset[index]
            bursty_indices = np.flatnonzero(np.array(read.peer.get_tag('pe')) == 2)
            bursty_gaps = np.diff(bursty_indices)
            bursty_breaks = np.flatnonzero(bursty_gaps > self.false_negative_burst_length)
            if bursty_indices.any():
                bursts['zmw'][burst_count] = read.holeNumber
                bursts['qStart'][burst_count] = read.qStart
                bursts['qEnd'][burst_count] = read.qEnd
                start_frames = read.peer.get_tag('sf')
                bursts['burstStart'][burst_count] = bursty_indices[0]
                p2b = fm.pls2base(fm.s2npl(read.peer.get_tag('pc')))
                index = bursty_indices[0] - 1
                previous_base_index = p2b[index]
                while (previous_base_index < 0) and (index >= 0):
                    index -= 1
                    previous_base_index = p2b[index]
                try:
                    bursts['previousBaseIndex'][burst_count] = previous_base_index
                    bursts['previousBasecall'][burst_count] = read.read(
                        aligned=False)[previous_base_index]
                except IndexError: # catch reads where there are no previous basecalls
                    bursts['previousBaseIndex'][burst_count] = -1
                    bursts['previousBasecall'][burst_count] = 'Z'
                if bursty_breaks.any():
                    for bursty_break in bursty_breaks:
                        bursts = self.resize_array(bursts, burst_count, self.subsampleto)
                        index = bursty_indices[bursty_break]
                        bursts['burstLength'][burst_count] = index - bursts[
                            'burstStart'][burst_count] + 1
                        bursts['burstStartTime'][burst_count] = start_frames[
                            bursts['burstStart'][burst_count]]
                        bursts['burstEndTime'][burst_count] = start_frames[
                            (bursts['burstStart'][burst_count] +
                             bursts['burstLength'][burst_count])]
                        bs = bursts['burstStart'][burst_count]
                        be = (bursts['burstStart'][burst_count] +
                              bursts['burstLength'][burst_count])
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
                        bursts = self.resize_array(bursts,
                                                   burst_count,
                                                   self.subsampleto)
                        bursts['zmw'][burst_count] = read.holeNumber
                        bursts['qStart'][burst_count] = read.qStart
                        bursts['qEnd'][burst_count] = read.qEnd
                        next_index = bursty_indices[bursty_break + 1]
                        bursts['burstStart'][burst_count] = next_index
                        index = next_index - 1
                        previous_base_index = p2b[index]
                        while (previous_base_index < 0) and (index >= 0):
                            index -= 1
                            previous_base_index = p2b[index]
                        bursts['previousBaseIndex'][burst_count] = previous_base_index
                        bursts['previousBasecall'][burst_count] = read.read(
                            aligned=False)[previous_base_index]
                bursts['burstLength'][burst_count] = (bursty_indices[-1] -
                                                      bursts['burstStart'][burst_count])
                bursts['burstStartTime'][burst_count] = start_frames[
                    bursts['burstStart'][burst_count]]
                bursts['burstEndTime'][burst_count] = start_frames[
                    (bursts['burstStart'][burst_count] +
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
