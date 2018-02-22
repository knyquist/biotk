import biotk.libs.BurstMetrics as bm
import os

class _test_burst_metrics:
    def __init__(self, sset_path, subsampleto=None):
        self.burst_metrics = bm.PpaBurstMetrics(sset_path,
                                                subsampleto=subsampleto)

def test_burst_metrics_initialization():
    """
    Given a subreadset, test that a BurstMetrics
    object can be initialized

    Returns a BurstMetrics object
    """
    path = os.path.dirname(os.path.abspath(__file__))
    sset_name = 'data/tiny_set_internal.subreadset.xml'
    sset_path = os.path.join(path, sset_name)
    bm = _test_burst_metrics(sset_path, subsampleto=50)


