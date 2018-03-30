import biotk.libs.YieldAccumulation as ya
import os

class _test_yield_accumulation:
    def __init__(self, aset_path):
        self.SequencingYield = ya.SequencingYield(aset_path)

def test_cmph5_alignments():
    path = os.path.dirname(os.path.abspath(__file__))
    aset = 'data/unrolled.aligned_reads.cmp.h5'
    aset_path = os.path.join(path, aset)
    test_seq_yield = _test_yield_accumulation(aset_path)
    assert test_seq_yield.SequencingYield.is_cmph5 is True
    test_seq_yield.SequencingYield.calculate_yield_by_time()
