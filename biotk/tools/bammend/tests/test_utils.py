import tools.bammend.utils as bu
import tempfile
import os

class _test_utils_module: # underscore because class is used by test defs below
    """Class for testing the functions in the utils module"""

    def __init__(self, path, input_bam, input_csv, output_bam):
        self.path = path
        self.tmp_path = tempfile.mkdtemp()
        self.input_bam = input_bam
        self.input_csv = input_csv
        self.output_bam = self.tmp_path + output_bam

    def delete_tmp_dir(self):
        """Remove files generated during test"""
        os.rmdir(self.tmp_path)

    def open_input_bam(self):
        """Tests opening input bam with pysam"""
        if bu.open_input_bam(self.input_bam):
            pass

    def prepare_output_bam(self):
        """Tests using input bam to template the output bam with pysam
           by confirming that headers are the same
        """
        i = bu.open_input_bam(self.input_bam)
        o = bu.prepare_output_bam(self.output_bam, i)
        assert i.header == o.header
        i.close()
        o.close()

    def generate_subreadset(self):
        """Tests that pbcore can open, index, and generate a SubreadSet
           from a pysam-generated Bamfile
        """
        i = bu.open_input_bam(self.input_bam)
        o = bu.prepare_output_bam(self.output_bam, i)
        i.close()
        o.close()
        if bu.generate_subreadset(self.output_bam):
            pass

def test_utils_internal_mode():
    """Test utils module with internal-mode data"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_bam = os.path.sep.join([path,
                                  'data',
                                  'internal_mode',
                                  'tiny_set_internal.subreads.bam'])
    input_csv = os.path.sep.join([path,
                                  'data',
                                  'internal_mode',
                                  'tiny_set_testfilter.csv'])
    output_bam = 'tiny_set_internal_output.subreads.bam'
    test = _test_utils_module(path, input_bam, input_csv, output_bam)
    test.open_input_bam()
    test.prepare_output_bam()
    test.generate_subreadset()
    test.delete_tmp_dir()

def test_utils_customer_mode():
    """Test utils module with customer-mode data"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_bam = os.path.sep.join([path,
                                  'data',
                                  'customer_mode',
                                  'tiny_set_customer.subreads.bam'])
    input_csv = os.path.sep.join([path,
                                  'data',
                                  'customer_mode',
                                  'tiny_set_testfilter.csv'])
    output_bam = 'tiny_set_customer_output.subreads.bam'
    test = _test_utils_module(path, input_bam, input_csv, output_bam)
    test.open_input_bam()
    test.prepare_output_bam()
    test.generate_subreadset()
    test.delete_tmp_dir()


