from biotk.libs.poa.ConsensusTensor import PoaConsensusTensorList
from biotk.libs.tests.test_PoaWithFeatures import TestPoaWithFeatures
from pbcore.io import (ReferenceSet, SubreadSet)
import numpy as np
from nose.tools import with_setup
import pandas as pd
import random


class TestPoaConsensusTensorList(TestPoaWithFeatures):
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
    test_poa_consensus_tensor_list = TestPoaConsensusTensorList('data/'
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
    poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads,
                                             collection_mode=mode1)
    assert poa_tensor_list.collection_mode == mode1
    poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads,
                                             collection_mode=mode2)
    assert poa_tensor_list.collection_mode == mode2
    try:
        poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads,
                                                 collection_mode=mode3)
    except ValueError:
        assert poa_tensor_list.collection_mode != mode3

@with_setup(setup_func)
def test_makeConsensusTensors():
    tpctl = setup_func()
    poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads)

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
    poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads,
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
    poa_tensor_list = PoaConsensusTensorList(tpctl.poa.subreads,
                                             ref=ref,
                                             context_width=1,
                                             collection_mode='standard',
                                             subsample_count=15)

def test_tmp():
    sset_path = '/pbi/dept/itg/zia/projects/3914/conditions/A/subreads/input.subreadset.xml'
    sset = SubreadSet(sset_path)
    rset_path = '/pbi/dept/enzymology/Kristofor/SoftwareProjects/published_tools/biochemistry-toolkit/biotk/libs/tests/data/references/All4mers_InsertOnly.ReferenceSet.xml'
    rset = ReferenceSet(rset_path)
    df = pd.DataFrame.from_records(sset.index)
    gb = df.groupby(by=['holeNumber'])
    df = gb.get_group(31195616)
    df[df['contextFlag'] == 3]
    subreads = [sset[i] for i in df.index]
    tmp = PoaConsensusTensorList(subreads, context_width=1, collection_mode='equal-state', subsample_count=10)