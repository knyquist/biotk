import tools.bammend.bammend as bm
import tempfile
import os

class _test_bammend_module: # underscore because class is used by test defs below
    """Class for testing the functions in the bammend module"""

    def __init__(self, path, input_bam, input_csv, output_bam):
        self.path = path
        self.tmp_path = tempfile.mkdtemp()
        self.input_bam = input_bam
        self.input_csv = input_csv
        self.output_bam = self.tmp_path + output_bam

    def delete_tmp_dir(self):
        """Remove files generated during test"""
        os.rmdir(self.tmp_path)

    def is_internal_mode(self, supposed_to_be_internal):
        """Tests whether detected internal-mode status is correct"""
        assert bm.is_internal_mode(self.input_bam) == supposed_to_be_internal

    def open_annotation_csv(self):
        """Tests whether basecall rejection CSV input can be opened
           and validated. Includes implicit test of bm.validate_bammend_csv
        """
        bm.open_annotation_csv(self.input_csv)

    def reject_basecalls(self):
        """Tests reject_basecalls routine. Includes implicit test of
           bm.filter_subread as well as bm.apply_subreadset_filters
        """
        if bm.reject_basecalls(self.input_bam, self.input_csv, self.output_bam):
            pass

def test_bammend_internal_mode():
    """Test bammend module with internal-mode data"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_bam = (path + '/data/internal_mode/tiny_set_' +
                        'internal.subreads.bam')
    input_csv = (path + '/data/internal_mode/tiny_set_' +
                        'testfilter.csv')
    output_bam = 'tiny_set_internal_output.subreads.bam'
    test = _test_bammend_module(path, input_bam, input_csv, output_bam)
    test.is_internal_mode(True)
    test.open_annotation_csv()
    test.reject_basecalls()
    test.delete_tmp_dir()

def test_bammend_customer_mode():
    """Test bammend module with customer-mode data"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_bam = (path + '/data/customer_mode/tiny_set_' +
                        'customer.subreads.bam')
    input_csv = (path + '/data/customer_mode/tiny_set_' +
                        'testfilter.csv')
    output_bam = 'tiny_set_customer_output.subreads.bam'
    test = _test_bammend_module(path, input_bam, input_csv, output_bam)
    test.is_internal_mode(False)
    test.open_annotation_csv()
    test.reject_basecalls()
    test.delete_tmp_dir()
