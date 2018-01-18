#!/usr/bin/env python

import numpy as np
import random
import string
import operator
from pbcore.io import (SubreadSet,
                       IndexedBamReader)
from pricompare import FastMetrics as fm

class PpaBurstMetrics:
    """
    Class for retrieving burst metrics. Two flavors
        Alignment based
        HMM Classifier based
        
    If required information not available, return None.
    """
    def __init__(self, subread_set_path,
                       subsampleto=None):
        self.subread_set_path = subread_set_path
        self.subread_set = SubreadSet(subread_set_path)
        self.scraps = IndexedBamReader(self.subread_set.externalResources[0].scraps)
        self.subsampleto = subsampleto
        self.zmws = self._subsample_zmws()
        (self.subread_ppa_bursts,
         self.subread_reads) = self.retrieve_classifier_bursts(self.subread_set,
                                                               'subreads')
        (self.scraps_ppa_bursts,
         self.scrap_reads) = self.retrieve_classifier_bursts(self.scraps,
                                                              'scraps')
        self.ppa_bursts = np.hstack((self.subread_ppa_bursts,
                                     self.scraps_ppa_bursts))
        self.reads = np.hstack((self.subread_reads,
                                self.scrap_reads))

    def _resize_array(self, arr, index, increase_by):
        """
        Resize NumPy array if necessary
        """
        if index >= len(arr): # extend array if needed
            new_size = tuple(map(operator.add, arr.shape, (increase_by, )))
            arr = np.resize(arr, new_size)
        return arr

    def _subsample_zmws(self):
        """
        Subsample Zmws for measurement
        """
        zmws = np.union1d(self.subread_set.index.holeNumber,
                          self.scraps.index.holeNumber) # scraps index bug should be fixed
        if self.subsampleto is not None:
            if len(zmws) > self.subsampleto:
                zmws = np.unique(random.sample(zmws,
                                               self.subsampleto))
        return zmws

    def retrieve_classifier_bursts(self, dset, dset_type):
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
                                                ('seqType', 'S1'), # seqType -> {H, L, A}
                                                ('burstStart', int),
                                                ('burstLength', int),
                                                ('numShorties', int),
                                                ('burstStartTime', int),
                                                ('burstEndTime', int),
                                                ('previousBasecall', 'S1'),
                                                ('previousBaseIndex', int),
                                                ('fractionC', float),
                                                ('fractionA', float),
                                                ('fractionT', float),
                                                ('fractionG', float)])
        burst_count = 0

        reads = np.zeros((len(dset), ), dtype=[('zmw', int),
                                               ('seqType', 'S1'),
                                               ('qStart', int),
                                               ('qEnd', int),
                                               ('startTime', int),
                                               ('endTime', int)])
        read_count = 0
        
        bases = ['a', 'c', 'g', 't']

        read_indices = np.flatnonzero(np.in1d(dset.index.holeNumber, self.zmws))
        for index in read_indices:
            read = dset[index]
            
            # Store information about the read being considered
            # Keep info even if read doesn't contain a burst 
            reads['zmw'][read_count] = read.holeNumber
            reads['qStart'][read_count] = read.qStart
            reads['qEnd'][read_count] = read.qEnd
            p2b = fm.pls2base(fm.s2npl(read.peer.get_tag('pc')))
            p2b = np.flatnonzero(
                fm.pls2base(fm.s2npl(read.peer.get_tag('pc'))) >= 0)
            start_frames = read.peer.get_tag('sf')
            reads['startTime'][read_count] = start_frames[p2b[0]]
            reads['endTime'][read_count] = start_frames[p2b[-1]]
            if dset_type == 'subreads':
                reads['seqType'] = 'H'
            elif dset_type == 'scraps':
                reads['seqType'] = read.scrapType
            read_count += 1

            # Consider read for bursts and record burst
            # information if they exist
            pe_reason = np.array(read.peer.get_tag('pe'))
            """convert short-frame exclusions that happen
            during bursts into burst exclusions"""
            shorties = np.zeros((len(pe_reason), ), dtype=int)
            for index in np.arange(1, len(pe_reason)):
                if pe_reason[index] == 1 and pe_reason[index-1] == 2:
                    pe_reason[index] = 2
                    shorties[index] = 1

            bursty_indices = np.flatnonzero(pe_reason == 2)
            bursty_gaps = np.diff(bursty_indices)
            bursty_breaks = np.flatnonzero(bursty_gaps > 1)

            if bursty_indices.any():
                if len(bursts) < burst_count + len(bursty_breaks) + 1:
                    # resize the bursts table
                    bursts = self._resize_array(bursts, burst_count, self.subsampleto)

                bursts['zmw'][burst_count] = read.holeNumber
                if dset_type == 'subreads':
                    bursts['seqType'] = 'H'
                elif dset_type == 'scraps':
                    bursts['seqType'] = read.scrapType
                else:
                    raise IOError('dset type must be either subreads or scraps')
                bursts['qStart'][burst_count] = read.qStart
                bursts['qEnd'][burst_count] = read.qEnd
                start_frames = read.peer.get_tag('sf')
                p2b = fm.pls2base(fm.s2npl(read.peer.get_tag('pc')))
                bursts['burstStart'][burst_count] = bursty_indices[0]
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
                        bursts['numShorties'][burst_count] = np.sum(shorties[bs:be])
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
                bursts['numShorties'][burst_count] = np.sum(shorties[bs:be])
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
            
        bursts = bursts[bursts['zmw'] != 0]
        reads = reads[reads['zmw'] != 0]
        return bursts, reads
