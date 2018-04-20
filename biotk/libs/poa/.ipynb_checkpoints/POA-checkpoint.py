import sys
sys.path.append('/pbi/dept/enzymology/Kristofor/repos/poapy')
import random
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
    def __init__(self, subreads):
        self.subreads = subreads

    def generatePoaGraph(self):
        """
        Given list of subreads, generate MSA using POA graphs
        """
        subreads = list(self.subreads)
        root_subread = subreads.pop()
        self.root_subread = root_subread
        root_seq = root_subread.read(aligned=False)
        root_label = root_subread.qName
        graph = poagraph.POAGraph(root_seq, label=root_label)
        for subread in subreads:
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