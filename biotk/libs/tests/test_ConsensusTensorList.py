from biotk.libs.poa.ConsensusTensor import ConsensusTensorList
from biotk.libs.tests.test_PoaWithFeatures import TestPoaWithFeatures
from pbcore.io import (ReferenceSet, SubreadSet)
import numpy as np
from nose.tools import with_setup


class TestConsensusTensorList(TestPoaWithFeatures):
    def __init__(self, sset_path,
                       context_width=0,
                       collection_mode='standard',
                       subsample_count=None):
        TestPoaWithFeatures.__init__(self, sset_path)
        self.context_width = context_width,
        self.collection_mode = collection_mode,
        self.subsample_count = subsample_count

def setup_func(ref=None):
    """
    set up ConsensusTensorList test object
    :return:
    """
    test_poa_consensus_tensor_list = TestConsensusTensorList('data/'
                                                             'tiny_set_internal.subreadset.xml',
                                                             ref)
    return test_poa_consensus_tensor_list

@with_setup(setup_func)
def test_check_collection_mode():
    """
    Test that collection mode only allows 'standard' or 'full-mode'
    :return:
    """
    mode1 = 'standard'
    mode2 = 'equal-state'
    mode3 = 'something else'
    tpctl = setup_func()
    poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads,
                                          collection_mode=mode1)
    assert poa_tensor_list.collection_mode == mode1
    poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads,
                                          collection_mode=mode2)
    assert poa_tensor_list.collection_mode == mode2
    try:
        poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads,
                                              collection_mode=mode3)
    except ValueError:
        print 'collection mode: ' + "'" + mode3 + "'" + ' is not supported'
        assert poa_tensor_list.collection_mode != mode3

@with_setup(setup_func)
def test_makeConsensusTensors():
    tpctl = setup_func()
    poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads)

@with_setup(setup_func)
def test_makeConsensusTensorsRef_EqualStates():
    """
    Test that ConsensusTensorList can be populated in equal-states
    mode
    :return:
    """
    rset = ReferenceSet('data/references/All4mers_InsertOnly.ReferenceSet.xml')
    ref = rset[np.flatnonzero(rset.index['id'] == 'All4mer.V2.105')[0]]
    tpctl = setup_func(ref)
    poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads,
                                          ref=ref,
                                          collection_mode='equal-state',
                                          subsample_count=30)

@with_setup(setup_func)
def test_makeConsensusTensorsRef_FullMode():
    """
    Test that ConsensusTensorList can be populated in standard mode
    :return:
    """
    rset = ReferenceSet('data/references/All4mers_InsertOnly.ReferenceSet.xml')
    ref = rset[np.flatnonzero(rset.index['id'] == 'All4mer.V2.105')[0]]
    tpctl = setup_func(ref)
    poa_tensor_list = ConsensusTensorList(tpctl.poa.subreads,
                                          ref=ref,
                                          context_width=1,
                                          collection_mode='standard',
                                          subsample_count=15)