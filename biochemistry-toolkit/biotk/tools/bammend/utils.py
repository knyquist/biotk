import os
import pysam
from pbcore.io import SubreadSet

def open_input_bam(input_bam_path):
    """Open bamfile with pysam."""
    bam_reader = pysam.Samfile(input_bam_path, 'r', check_sq=False)
    return bam_reader

def prepare_output_bam(output_bam_path, template_bam_reader):
    """Prepare output bamfile with pysam, using original bam as template."""
    output_bam = pysam.Samfile(output_bam_path,
                               'wb',
                               template=template_bam_reader)
    return output_bam

def generate_subreadset(output_bam_path):
    """Save output bamfile and subreadset with pbcore."""
    sset = SubreadSet(output_bam_path, generateIndices=True)
    sset_output_name = ('.'.join(output_bam_path.split('.')[0:-1]) +
                        '.bammended.subreadset.xml')
    sset.newUuid()
    sset.write(sset_output_name)
    return True
