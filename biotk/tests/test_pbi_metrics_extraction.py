import biotk.scripts.extractIndexMetrics as im
from pbcore.io import AlignmentSet
import os

class _test_index_metrics:
    def __init__(self, aset_path, output_path, mode):
        self.aset_path = aset_path
        self.aset = AlignmentSet(self.aset_path)
        self.output_path = output_path
        self.mode = mode # operation mode, summarize or full-mode
        self.generatePbiCsv = im.generatePbiCsv(self.aset_path,
                                                self.output_path,
                                                self.mode)

def test_index_metrics_full():
    """
    Given an alignmentset, test that the PBI index metrics
    can be saved to a CSV file
    """
    path = os.path.dirname(os.path.abspath(__file__))
    aset_name = 'data/tinyset.alignmentset.xml'
    aset_path = os.path.join(path, aset_name)
    output_path = '/tmp/pbi.csv'
    tim = _test_index_metrics(aset_path, output_path, 'full-mode')

def test_index_metrics_summarize():
    """
    Given an alignmentset, test that the PBI index metrics
    can be saved to a CSV file
    """
    path = os.path.dirname(os.path.abspath(__file__))
    aset_name = 'data/tinyset.alignmentset.xml'
    aset_path = os.path.join(path, aset_name)
    output_path = '/tmp/pbi.csv'
    tim = _test_index_metrics(aset_path, output_path, 'summarize')
