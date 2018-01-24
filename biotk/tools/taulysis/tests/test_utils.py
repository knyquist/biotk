import biotk.tools.taulysis.utils as tu
import tempfile
import os

class _test_utils_module: # underscore because class is used by test defs below
    """Class for testing the functions in the utils module"""

    def __init__(self, path,
                       alignment_set_path,
                       basefile_path=None,
                       is_legacy=False): # default to Sequel data
        self.path = path
        self.tmp_path = tempfile.mkdtemp()
        self.alignment_set_path = alignment_set_path
        self.basefile_path = basefile_path
        self.is_legacy = is_legacy

    def open_alignment_set(self):
        """Test opening alignment set with pbcore"""
        if tu.open_alignment_set(self.alignment_set_path, self.is_legacy):
            pass

    def get_smrtbellsize(self):
        """Test retrieving smrtbellsize from reference
           information in alignment set
        """
        aset = tu.open_alignment_set(self.alignment_set_path, self.is_legacy)
        reference_info = tu.get_smrtbellsize(aset)
        return reference_info

    def override_movie_limited_flag(self, should_override):
        """Test movie-limited override function"""
        aset = tu.open_alignment_set(self.alignment_set_path, self.is_legacy)
        assert (tu.override_movie_limited_flag(aset, self.is_legacy) ==
                should_override)

    def is_internal_mode(self, is_internal):
        """Test internal-mode check"""
        aset = tu.open_alignment_set(self.alignment_set_path, self.is_legacy)
        assert tu.is_internal_mode(aset, self.is_legacy) == is_internal

    def open_base_file(self):
        """Test opening basefile, if provided"""
        if self.basefile_path is not None:
            tu.open_base_file(self.basefile_path)
        else:
            pass # basefile was not provided

def test_utils_legacy_RS():
    """Test utils module with legacy RS data"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_cmph5 = os.path.sep.join([path,
                                    'data',
                                    'legacy_RS',
                                    'unrolled.aligned_reads.cmp.h5'])
    input_bash5 = os.path.sep.join([path,
                                    'data',
                                    'legacy_RS',
                                    'unrolled.reads.bas.h5'])
    is_legacy = True
    test = _test_utils_module(path, input_cmph5, input_bash5, is_legacy)
    test.open_alignment_set()
    assert test.get_smrtbellsize()[0]['SMRTBellSize'] == 698
    test.override_movie_limited_flag(should_override=False)
    test.is_internal_mode(is_internal=True)
    test.open_base_file()

def test_utils_internal_bam():
    """Test utils module with internal-mode bam"""
    path = os.path.dirname(os.path.abspath(__file__))
    input_bam = os.path.sep.join([path,
                                  'data',
                                  'internal_mode_bam',
                                  ('tiny_set_unrolled_internal_mode.'
                                   'alignmentset.bam')])
    test = _test_utils_module(path, input_bam)
    test.open_alignment_set()
    assert test.get_smrtbellsize()[0]['SMRTBellSize'] == 698
    test.override_movie_limited_flag(should_override=False)
    test.is_internal_mode(is_internal=True)
    test.open_base_file()

def test_utils_customer_bam():
    """Test utils module with customer bam"""
    # need a good test customer bam
    pass

