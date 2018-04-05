from biotk.libs.poa import POA as poa
from pbcore.io import (SubreadSet,
                       ReferenceSet)
import pandas as pd
import numpy as np
from nose.tools import with_setup


class TestPOA:

    def __init__(self, sset_path, ref=None):
        """
        Constructor for testing the POA classes
        :param sset_path:
        """
        self.sset_path = sset_path
        self.sset = SubreadSet(self.sset_path)
        self.zmw_subreads = self._retrieve_single_zmw_subreads()
        self.poa = poa.POA(self.zmw_subreads, ref)

    def _retrieve_single_zmw_subreads(self):
        """
        Given a subreadset, return list of
        subreads for particular zmw
        :return:
        """
        df = pd.DataFrame.from_records(self.sset.index)
        gb = df.groupby(by=['holeNumber'])
        group = list(gb)[-1]  # last group has 4 subreads
        zmw = group[0]
        subreads_ixs = list(group[1].index)
        subreads = self.sset[subreads_ixs]
        ccs_subreads = []  # use only adapter-flanked subreads
        for subread in subreads:
            if subread.contextFlag == 3:
                ccs_subreads.append(subread)
        return ccs_subreads

def setup_func(ref=None):
    """
    set up POA object
    :return:
    """
    test_poa = TestPOA('data/tiny_set_internal.subreadset.xml',
                       ref)
    return test_poa

@with_setup(setup_func)
def test_generatePoaGraph():
    """
    Demonstrate that a poa-generated MSA
    can be generated from a list of subreads
    from a particular ZMW (using pbcore)
    :return:
    """
    # no reference provided
    test_poa = setup_func()
    test_poa.poa.generatePoaGraph()

    # reference provided
    rset = ReferenceSet('data/references/All4mers_InsertOnly.ReferenceSet.xml')
    ref = rset[np.flatnonzero(rset.index['id'] == 'All4mer.V2.105')[0]]
    test_poa = setup_func(ref=ref)
    test_poa.poa.generatePoaGraph()

@with_setup(setup_func)
def test_check_direction():
    """
    Test that levenshtein distance can be used
    to test whether a sequence should be reverse
    complemented before folded into the POA
    :return:
    """
    test_poa = setup_func()
    seq1 = 'ACTGCGG'
    seq2 = 'ACTGCGG'
    seq2_rc = 'CCGCAGT'
    # case 1: fwd
    assert seq1 == test_poa.poa._check_direction(seq1, seq2)
    # case 2: rev
    assert seq1 == test_poa.poa._check_direction(seq2_rc, seq2)

@with_setup(setup_func)
def test_reverse_complement():
    """
    Test that sequence can be reverse complemented
    :return:
    """
    test_poa = setup_func()
    seq = 'ACTG'
    rc = 'CAGT'
    assert rc == test_poa.poa._reverse_complement(seq)
