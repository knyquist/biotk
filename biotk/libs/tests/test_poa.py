from biotk.libs.poa import POA as poa
from pbcore.io import SubreadSet
import pandas as pd
import os

class _test_poa:
    def __init__(self, sset_path):
        """
        Constructor for testing the POA classes
        :param sset_path:
        """
        self.sset_path = sset_path
        self.sset = SubreadSet(self.sset_path)
        self.zmw_subreads = self._retrieve_single_zmw_subreads()
        self.poa = poa.POA(self.zmw_subreads)

    def _retrieve_single_zmw_subreads(self):
        """
        Given a subreadset, return a list of
        subreads for a particular ZMW
        """
        df = pd.DataFrame.from_records(self.sset.index)
        gb = df.groupby(by=['holeNumber'])

        group = list(gb)[-1] # the last group happens to have 4 subreads
        zmw = group[0]
        subread_ixs = list(group[1].index)
        zmw_subreads = self.sset[subread_ixs]
        subreads = []
        for subread in zmw_subreads:
            if subread.contextFlag == 3:
                subreads.append(subread)

        return subreads

def test_poa():
    """
    Given a set of subreads, demonstrate that a multiple-sequence
    alignment can be generated
    """
    test_poa = _test_poa('data/tiny_set_internal.subreadset.xml')
    assert len(test_poa.zmw_subreads) == 4
    graph = poa.POA(test_poa.zmw_subreads).generatePoaGraph()

def test_feature_fold_in():
    """
    Test that features, such as PW and IPD, can be folded into
    the MSA
    """
    test_poa = _test_poa('data/tiny_set_internal.subreadset.xml')
    assert len(test_poa.zmw_subreads) == 4
    poa.PoaWithFeatures(test_poa.zmw_subreads)