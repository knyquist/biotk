import numpy as np
from numpy.lib.recfunctions import append_fields as append_to_recarray
from pbcore.io import (AlignmentSet,
                       BasH5Reader)
import logging
logging.basicConfig(level=logging.INFO)

def open_base_file(basefile_path):
    """Open basefile, if using cmp.h5 legacy format
       and the the basefile path was provided.
    """
    bas_reader = BasH5Reader(basefile_path)
    return bas_reader

def is_internal_mode(alignment_set, is_legacy):
    """See if alignmentset was collected in internal-mode."""
    if is_legacy:
        is_internal = True
    else:
        is_internal = 'PulseCall' in alignment_set.pulseFeaturesAvailable()
    return is_internal

def override_movie_limited_flag(alignment_set, is_legacy):
    """Checks to see if alignment set was collected with internal-mode.
       If not, analysis that includes the movie limited flag will not
       be performed, even if the movie-length was user inputted.
    """
    is_internal = is_internal_mode(alignment_set, is_legacy)
    if is_internal:
        logging.info('Alignmentset contains internal-mode features')
        return False
    else:
        logging.info('Alignmentset was collected without internal-mode. '
                     'It is not possible to determine which alignments '
                     'were movie-limited. Analysis will be done without '
                     'tracking movie-limited alignments, even if movie '
                     'length was inputted.')
        return True

def open_alignment_set(alignments_path, is_legacy):
    """Open AlignmentSet with pbcore"""
    aset = AlignmentSet(alignments_path)
    if 'RefGroupID' in aset.index.dtype.names:
        is_legacy_rpt = True
        logging.info('Alignmentset found to be cmp.h5 format.')
    else:
        is_legacy_rpt = False
        logging.info('Alignmentset found to be bam format.')

    is_file_type_correct = is_legacy == is_legacy_rpt
    if not is_file_type_correct:
        raise AssertionError('The aligmentset does not match expectation. '
                             'Is it a bam file and the --legacy, i.e. '
                             'cmp.h5, flag was enabled (or vice-versa)?')
    return aset

def get_smrtbellsize(aset):
    """Retrieve SMRTBell size from reference info in AlignmentSet"""
    reference_info_table = aset.referenceInfoTable
    reference_information = reference_info_table[['ID', 'FullName']]
    logging.info('If FutureWarning popped up on previous line, it is bullshit.')
    sbell_sizes = []
    for reference_id, full_reference_name in reference_information:
        rsplit = full_reference_name.split('_unrolled_circular_')
        if len(rsplit) == 1: # try other way
            # I don't like doing this, but the reference list is contaminated
            rsplit = full_reference_name.split('_circular_unrolled_')
        if len(rsplit) == 1: # try final way, again I hate doing this
            rsplit = full_reference_name.split('_circular_')
        if len(rsplit) == 2:
            ref_name, ref_size = rsplit
        else:
            raise ValueError('Reference names do not follow unrolled style. '
                             'Splitting reference name by _unrolled_circular_, ' 
                             '_circular_unrolled_, or _circular_, should '
                             'produce two substrings, one with the '
                             'name and the other with length info.')
        x_unrolled = int(ref_size.split('x')[0])
        total_len = int(ref_size.split('l')[1])
        size = np.divide(total_len, x_unrolled, dtype=int)
        sbell_sizes.append(size)
    sbell_sizes = np.array(sbell_sizes, dtype=int)
    reference_information = append_to_recarray(reference_information.copy(),
                                               'SMRTBellSize',
                                               data=sbell_sizes,
                                               asrecarray=True)
    return reference_information
