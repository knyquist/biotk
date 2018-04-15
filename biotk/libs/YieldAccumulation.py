from pbcore.io import (CmpH5Reader,
                       AlignmentSet)
import numpy as np
import os

class SequencingYield:
    """
    Class for characterizing the yield of a 
    sequencing run
    """

    def __init__(self, aset_path):
        (self.aset,
         self.is_cmph5) = self._openAset(aset_path)

    def calculate_yield_by_time(self):
        """
        Return yield vs. time. Calculation will be different
        depending on whether data are from cmp.h5 vs bam
        """
        if self.is_cmph5:
            time, base_yield = self.ybt_cmph5()
        else:
            time, base_yield = self.ybt_bam()

    def ybt_cmph5(self):
        """
        Return yield vs. time for cmph5 datasets
        """
        min_time, max_time, time_interval = 0, 1800, 5
        time = np.arange(min_time, max_time, time_interval)
        yield_per_time = np.zeros(time.shape, dtype=int)
        for alignment in self.aset:
            advance_time = alignment.IPD() + alignment.PulseWidth()
            advance_time[advance_time == 65534] = np.round(
                np.mean(
                    advance_time[
                        advance_time != 65534]))
            advance_time = np.divide(advance_time,
                                     self.aset.readGroupTable['FrameRate'] * 60,
                                     dtype=float)
            start_frame = np.cumsum(advance_time)
            counts, bin_edges = np.histogram(start_frame, time)
            yield_per_time[0:-1] += counts
        cumulative_yield = np.cumsum(yield_per_time)
        time = time + 0.5 * time_interval
        max_index = np.argmax(cumulative_yield)
        return time[0:max_index+1], cumulative_yield[0:max_index+1]

    def _openAset(self, aset_path):
        ext = os.path.splitext(aset_path)[-1]
        if ext == '.h5':
            return self._openCmpH5(aset_path)
        elif ext == '.xml':
            return self._openAlignmentSet(aset_path)
        else:
            raise IOError('Did not recognize filename extension') 

    def _openCmpH5(self, aset_path):
        print aset_path
        return CmpH5Reader(aset_path), True

    def _openAlignmentSet(self, aset_path):
        return AlignmentSet(aset_path), False
