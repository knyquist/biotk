from biotk.libs.poa import POA as poa
from biotk.libs.tests.test_POA import TestPOA
from nose.tools import with_setup

class TestPoaWithFeatures(TestPOA):
    def connectIpdAndPw(self):
        poa_with_features = poa.PoaWithFeatures(self.zmw_subreads)
        return poa_with_features

def setup_func():
    """
    set up POA object
    :return:
    """
    test_poa = TestPoaWithFeatures('data/tiny_set_internal.subreadset.xml')
    return test_poa

@with_setup(setup_func)
def test_foldInFeatures():
    """
    Test that POA can be generated and IPD
    and PW can be attached to each base in
    their respective alignments.
    """
    tpwf = setup_func()
    tpwf.connectIpdAndPw().foldInFeatures()

