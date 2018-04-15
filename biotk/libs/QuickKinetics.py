from pbcore.io import SubreadSet
import pandas as pd
import numpy as np
import operator

class kinetics:
    """
    Quickly extract representative IPDs and PWs
    from a subreadset
    """
    def __init__(self, sset_path,
                       nreads=None,
                       samples_per_read=None,
                       unique_zmws=False):
        self.sset_path = sset_path
        self.sset = SubreadSet(self.sset_path)
        self.nreads = nreads
        self.samples_per_read = samples_per_read
        self.unique_zmws = unique_zmws

    def _getUniqueSubreadIndices(self, index):
        fn = lambda df: np.random.choice(df.index, size=1, replace=False)
        indices = index.groupby(by=['holeNumber'], as_index=False).apply(fn)
        indices = np.concatenate(indices).ravel()
        return index.loc[indices]

    def _getSubreadIndices(self, index):
        if self.unique_zmws == True:
            # randomly choose subread from each ZMW to get kinetics
            index = self._getUniqueSubreadIndices(index)

        if self.nreads != None:
            # randomly subsample reads to get kinetics from
            if len(index) > self.nreads:
                index = index.sample(self.nreads)

        return index.index.values


    def summarizeKinetics(self):
        # get indices of subreads to scrape kinetics from
        index = pd.DataFrame.from_records(self.sset.index)

        # We now have an index struct that contains .pbi info
        # If the unique_zmws flag was triggered, it only
        # contains one randomly selected read per ZMW. If
        # the nreads flag was triggered, we randomly subsampled.
        # Now we iterate through the read indices, and extract
        # kinetics.
        sampled_read_ixs = self._getSubreadIndices(index)
        if self.samples_per_read != None:
            storage_per_read = self.samples_per_read
        else:
            storage_per_read = 2000
        kinetics_table = np.zeros(len(sampled_read_ixs)*storage_per_read,
                                  dtype=[('base', 'S1'),
                                         ('IPD', int),
                                         ('PW', int)])
        cntr = 0
        for j, read_ix in enumerate(sampled_read_ixs):
            read = self.sset[read_ix]
            kin = np.array(zip(list(read.read(aligned=False)),
                           read.IPD(aligned=False),
                           read.PulseWidth(aligned=False)),
                           dtype=[('base', 'S1'),
                                  ('IPD', int),
                                  ('PW', int)])
            if self.samples_per_read != None:
                kin = np.random.choice(kin, self.samples_per_read)
            l = len(kin)
            if l + cntr > kinetics_table.shape[0]:
                kinetics_table = self._resize_array(kinetics_table,
                                                    (len(sampled_read_ixs)-j)*self.samples_per_read)
            kinetics_table[cntr:cntr+l] = kin
            cntr += kin.shape[0]

        kinetics_summary = pd.DataFrame.from_records(kinetics_table)
        gb = kinetics_summary.groupby(by='base')
        mean = gb.mean()
        median = gb.median()
        std = gb.std()
        std.rename({'IPD': 'IPD_std',
                    'PW': 'PW_std'},
                   axis='columns',
                   inplace=True)
        summary = mean.merge(median, left_index=True, right_index=True, suffixes=('_mean', '_median'), )
        summary = summary.merge(std, left_index=True, right_index=True)

        return summary

    def _resize_array(self, arr, increase_by):
        """
        Resize NumPy array if necessary
        """
        new_size = tuple(map(operator.add, arr.shape, (increase_by, )))
        arr = np.resize(arr, new_size)
        return arr