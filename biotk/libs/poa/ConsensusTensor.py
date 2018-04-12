import POA as poa
import numpy as np
import logging

logging.basicConfig()
log = logging.getLogger(__name__)

class PoaConsensusTensorList(poa.PoaWithFeatures):
    """
    Class constructs list of ConsensusTensor objects from list of
    subreads
    """
    def __init__(self, subreads,
                       ref=None,
                       context_width=0,
                       collection_mode='standard',
                       subsample_count=None):
        """
        Initialize PoaConsensusTensorList object. Each argument has default values.
        Defaults to no contexts, standard collection, and no subsampling.

        :param context_width: see definition in ConsensusTensor class docstring
        :param collection_mode: two options, 'standard' and 'equal-state'. The 'equal-state'
                                mode returns a roughly equal proportion of {A, T, G, C, -}
                                examples. Standard mode makes no attempt to do so.
        :param subsample_count: returns a particular number of randomly-selected
                                ConsensusTensor objects. If more samples are requested than
                                exist, returns max possible number.
        """
        poa.PoaWithFeatures.__init__(self, subreads,
                                           ref)
        self.context_width = context_width
        self.collection_mode = collection_mode
        self._check_collection_mode()
        self.subsample_count = self._subsample_count(subsample_count)
        self.consensus_tensor_list = self.makeConsensusTensors()

    def _subsample_count(self, subsample_count):
        """
        Populate class value of subsample count. If None, make max length
        :return:
        """
        if subsample_count is None:
            return len(self.PluralityConsensus[1])
        else:
            return subsample_count

    def _check_collection_mode(self):
        """
        Make sure collection is either 'standard' or 'equal-state'. No other options
        are supported.
        :return:
        """
        if (self.collection_mode != 'standard') and (self.collection_mode != 'equal-state'):
            raise ValueError("Collection mode must be either 'standard' or 'equal-state'. "
                             "No other options are currently supported.")

    def _get_tensor_loci(self):
        """
        Use subsample_count and collection_mode to return loci indices
        where tensor info should be populated
        :return:
        """
        # acceptable loci must be within one context_width of the
        # edge of the MSA. Otherwise, the tensor can't be fully
        # populated
        start = self.context_width  # inclusive
        end = len(self.PluralityConsensus[1]) - self.context_width  # inclusive
        seq = np.array(list(self.PluralityConsensus[1][start:end]))

        # equal-state mode does its best to balance the proportion
        # of each class. May be useful for training classifiers
        if self.collection_mode == 'equal-state':
            if self.subsample_count < (end - start):
                n = self.subsample_count / 5  # number of tensors to grab of each type
            else:
                n = (end - start) / 5
            tensor_types = ['A', 'C', 'G', 'T', '-']
            loci = np.zeros((n*5, ), dtype=int)
            loci_ix = 0
            for type in tensor_types:
                ixs = np.flatnonzero(seq == type)
                if ixs.size == 0:
                    continue
                if n > len(ixs):
                    to_replace = True
                else:
                    to_replace = False
                ixs = np.random.choice(ixs, n, replace=to_replace)
                loci[loci_ix:loci_ix+len(ixs)] = ixs
                loci_ix = loci_ix + len(ixs)

        # standard mode just returns the classes as they
        # appear. No balancing.
        elif self.collection_mode == 'standard':
            if self.subsample_count < (end - start):
                n = self.subsample_count
                loci = np.random.choice(np.arange(start, end, 1), n, replace=False)
            else:
                loci = np.arange(start, end, 1)

        return loci


    def makeConsensusTensors(self):
        """
        Given a partial-order alignment with IPD and PW features, construct
        list of ConsensusTensors

        subsample_count and collection_mode are invoked first, to return the indices
        of the bases which should have their tensors recorded.
        :return:
        """
        loci = self._get_tensor_loci()
        if self.refMSA is not None:
            labels = np.array(list(self.refMSA[1]))[loci]
        else:
            labels = np.array(list(self.PluralityConsensus[1]))[loci]
        tensors = np.empty((len(loci), ), dtype=object)
        for index, locus in enumerate(loci):
            ixs = np.arange(locus - self.context_width,
                            locus + self.context_width + 1)
            data = self.feature_vector[1][:, ixs, :]
            label = labels[index]
            tensors[index] = ConsensusTensor(data=data, label=label)

        return tensors

class ConsensusTensor:
    """
    Class defining tensor summary data
    objects for doing consensus calling

                /--------------------|
               / IPD duration (cum.)/|
              /--------------------/ |
             / PW duration (cum.) /| |
            /--------------------/ | |
           / fraction calls     /| | |
    ---   /--------------------/ | | |
     |    | -                 |  | | |
     |    |-------------------|  | | |
     |    | A                 |  | | /  --
     |    |-------------------|  | |/   /
     y    | T                 |  | /   /
     |    |-------------------|  |/   z
     |    | G                 |  /   /
     |    |------------------ | /   /
     |    | C                 |/   /
    ---   |-------------------/   --

          |--------- x --------|

    dimension of x set by context width (0='no flanking bases',
                                         1='1 flanking base on each side',
                                         2='2 flanking bases on each side',
                                         etc.)
    dimension of y set by number of bases (always 5, including a deleted base)
    dimension of z set by number of features (i.e. fraction, IPD, and PW)

    The tensor is structured (row, col, depth) -> (y, x, z), as a normal
    numpy array

    The cumulative durations (IPD and PW) are z-scored according to by-base
    kinetic distributions
    """
    def __init__(self, data,
                       label):
        """
        Initialize ConsensusTensor object
        """
        self.tensor = self._populateTensor(data)  # np.zeros((5, 1 + 2 * context_width, 3))
        self.label = label

    def _populateTensor(self, data):
        """
        Use the data to populate and return the consensus tensor
        :return:
        """
        nrows = 5  # number of states {'-', 'A', 'T', 'G', 'C'}
        ncols = data.shape[1]  # 1 + 2 * context_width
        nlayers = data.shape[2]  # number of features (fraction, IPD, PW)
        tensor = np.zeros((nrows, ncols, nlayers), dtype=float)
        for col_index, col in enumerate(data[:, :, 0].T):
            bases, indices, counts = np.unique(col, return_inverse=True, return_counts=True)
            tensor[bases, col_index, 0] = np.divide(counts, np.sum(counts), dtype=float)
            for base in bases:
                row_indices = np.flatnonzero(data[:, col_index, 0] == base)
                tensor[base, col_index, 1] = np.sum(data[row_indices, col_index, 1])
                tensor[base, col_index, 2] = np.sum(data[row_indices, col_index, 2])

        return tensor