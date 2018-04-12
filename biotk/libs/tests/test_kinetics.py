from biotk.libs.QuickKinetics import kinetics
import pandas as pd

class TestKinetics:

    def __init__(self, sset_path,
                       nreads=None,
                       samples_per_read=None,
                       unique_zmws=False):
        self.sset_path = sset_path
        self.kinetics = kinetics(self.sset_path,
                                 nreads,
                                 samples_per_read,
                                 unique_zmws)

def setup_func(nreads=None,
               samples_per_read=None,
               unique_zmws=False):
    test_kin = TestKinetics('data/tiny_set_internal.subreadset.xml',
                            nreads,
                            samples_per_read,
                            unique_zmws)
    return test_kin

def test__getSubreadIndices():
    test_kin = setup_func(nreads=10,
                          samples_per_read=2,
                          unique_zmws=True)
    index = pd.DataFrame.from_records(test_kin.kinetics.sset.index)
    ixs = test_kin.kinetics._getSubreadIndices(index)

def test_getKinetics():
    test_kin = setup_func(nreads=10,
                          unique_zmws=False)
    test_kin.kinetics.getKinetics()
