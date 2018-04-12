import sys
sys.path.append('/pbi/dept/enzymology/Kristofor/repos/poapy') # fix this
import numpy as np
import poagraph
import seqgraphalignment
import logging
logging.basicConfig()
log = logging.getLogger(__name__)

try: # try to import ubuntu-compiled ld, then mac, then pure python one (slow)
    import ubuntu.levenshtein_distance as ld
    log.info('Running ubuntu-compiled levenshtein code...')
except ImportError: # ubuntu errored
    try:
        import mac.levenshtein_distance as ld
        log.info('Running mac-compiled levenshtein code...')
    except ImportError: # mac one errored
        import levenshtein_distance as ld
        log.info('Running pure python levenshtein code. Slow!')


class POA:
    """
    Generate partial-order alignment from collection of subreads
    """
    def __init__(self, subreads,
                       ref=None):
        self.subreads = subreads
        # if ref is populated, perform MSA against the
        # reference sequence. Otherwise align against
        # the last subread
        self.reference = ref

    def generatePoaGraph(self):
        """
        Given list of subreads, generate MSA using POA graphs
        """
        subreads = list(self.subreads)
        if self.reference:
            # the seeded sequence is the reference
            root_subread = self.reference
            self.root_subread = root_subread
            root_seq = str(root_subread.sequence)
            root_label = root_subread.header
        else:
            # no reference provided, align against subread
            root_subread = subreads.pop()
            self.root_subread = root_subread
            root_seq = root_subread.read(aligned=False)
            root_label = root_subread.qName
        graph = poagraph.POAGraph(root_seq, label=root_label)
        for subread in subreads:
            # uses levenshtein distance to determine if sequence should
            # be reverse-complemented before being added to the POA MSA
            subread_seq = self._check_direction(subread.read(aligned=False), root_seq)

            subread_label = subread.qName
            alignment = seqgraphalignment.SeqGraphAlignment(subread_seq,
                                                            graph,
                                                            fastMethod=True,
                                                            globalAlign=True)
            graph.incorporateSeqAlignment(alignment, subread_seq, label=subread_label)

        return graph

    def _check_direction(self, seq, root_seq):
        """
        Use Levenshtein distance to check whether seq should be 
        reverse-complemented before aligning to root_seq
        """
        og_score = ld.levenshtein_distance(seq, root_seq)
        rc = self._reverse_complement(seq)
        rc_score = ld.levenshtein_distance(rc, root_seq)
        if og_score < rc_score:
            return seq
        else:
            return rc

    def _reverse_complement(self, seq):
        """
        return reverse-complement of sequence
        """
        complement = {'A': 'T',
                      'G': 'C',
                      'T': 'A',
                      'C': 'G'}
        rev = seq[::-1]
        rc = ''.join([complement[base] for base in rev])
        return rc


class PoaWithFeatures(POA):
    """
    Takes as input adapter-flanked subreads from
    a particular ZMW, produces a multiple-sequence
    alignment using partial-order graphs, and ties
    in IPD and PW info to each raw subread.
    """
    def __init__(self, subreads,
                       ref=None):
        POA.__init__(self, subreads,
                           ref)
        self.subread_names = np.array([s.readName for s in self.subreads])
        self.PoaGraph = self.generatePoaGraph()  # perform POA MSA
        self.PoaStrings = self.PoaGraph.generateAlignmentStrings()  # convert graph to strings
        # only allow one consensus sequence
        while self.PoaStrings[-1][0] != 'Consensus0':
            self.PoaStrings.pop()

        if ref is not None:
            self.MSAs = self.PoaStrings[1:-1]
            self.refMSA = self.PoaStrings[0]
        else:
            self.MSAs = self.PoaStrings[0:-1]
            self.refMSA = None
        self.PluralityConsensus = self.PoaStrings[-1]
        self.feature_vector = self.foldInFeatures()
        self.baseOrder = ['-ATGC']

    def _baseOrder(self, base):
        """
        Construct dictionary encoding base string to integer
        :return:
        """
        base_order = {'-': 0, 'A': 1, 'T': 2, 'G': 3, 'C': 4}
        return base_order[base]

    def foldInFeatures(self):
        """
        For each alignment, connect the by-base features.
        Each deleted base will have 0 stored for each feature.
        By-base features include PW and IPD.

        The first layer in the feature_vector encodes the base-identity
        The second layer in the feature_vector encodes the PW (in frames)
        The third layer in the feature_vector encodes the IPD (in frames)
        """
        if len(self.MSAs) > 0:
            feature_vector = np.zeros((len(self.MSAs), len(self.MSAs[0][1]), 3), dtype=int)
            read_labels = []
        for row_index, msa in enumerate(self.MSAs):
            read_name = msa[0]
            read_labels.append(read_name)
            sequence = msa[1]
            base_ixs = np.flatnonzero(np.array(list(sequence)) != '-')
            raw_subread_ix = np.flatnonzero(self.subread_names == read_name)[0]
            subread = self.subreads[raw_subread_ix]
            feature_vector[row_index, :, 0] = map(self._baseOrder, sequence)
            feature_vector[row_index, base_ixs, 1] = subread.PulseWidth(aligned=False)
            feature_vector[row_index, base_ixs, 2] = subread.IPD(aligned=False)

        return read_labels, feature_vector