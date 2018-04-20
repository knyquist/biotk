import utils as bu
from pbcore.io import SubreadSet
import csv
import numpy as np
import array
import itertools
import pricompare.FastMetrics as fm
import logging
logging.basicConfig(level=logging.INFO)

def is_internal_mode(subreadset_path):
    """See if subreadset was collected in internal-mode."""
    sset = SubreadSet(subreadset_path)
    is_internal = 'PulseCall' in sset.pulseFeaturesAvailable()
    return is_internal

def validate_annotation_csv(annotation_csv_reader):
    """Validate that bammend csv has correct information"""
    if sorted(list(annotation_csv_reader.next().keys())) == sorted(['holeNumber',
                                                                    'annotationStartIndex',
                                                                    'annotationEndIndex']):
        return True
    else:
        raise KeyError(('Basecall-rejection CSV does not have '
                        'correct column names'))

def open_annotation_csv(annotation_csv_path):
    """Open CSV file marking indices to reject basecalls."""
    annotation_csv = {}
    with open(annotation_csv_path, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
        validate_annotation_csv(reader)
        for row in reader:
            zmw = int(row['holeNumber'])
            start = int(row['annotationStartIndex'])
            end = int(row['annotationEndIndex'])
            if zmw not in annotation_csv:
                annotation_csv[zmw] = {'annotationStartIndex': [],
                                       'annotationEndIndex': []}
            annotation_csv[zmw]['annotationStartIndex'].append(start)
            annotation_csv[zmw]['annotationEndIndex'].append(end)
    
    for key in annotation_csv:
        annotation_csv[key]['annotationStartIndex'] = np.array(annotation_csv[key]['annotationStartIndex'])
        annotation_csv[key]['annotationEndIndex'] = np.array(annotation_csv[key]['annotationEndIndex'])

    return annotation_csv


def filter_subread(subread,
                   rejected_query_indices,
                   is_internal):
    """Filter subread by query indices"""
    if is_internal:
        by_base_tags = ['dq', 'dt', 'iq', 'mq', 'sq', 'st', 'ip', 'pw']
    else:
        by_base_tags = ['ip', 'pw']

    read_query_indices = fm.query2bases(subread.peer.seq,
                                        subread.qStart)
    read_base_indices = read_query_indices - subread.qStart
    rejected_base_indices = np.flatnonzero(
                                np.in1d(read_query_indices,
                                        rejected_query_indices))
    accepted_base_indices = np.setdiff1d(read_base_indices,
                                         rejected_base_indices)
    for tag in by_base_tags:
        data = subread.peer.get_tag(tag)
        if type(data[0]) == str:
            accepted_data = [data[i] for i in accepted_base_indices]
            accepted_data = ''.join(accepted_data)
            subread.peer.set_tag(tag, accepted_data)
        else:
            data_typecode = data.typecode
            data = np.array(data)
            accepted_data = data[accepted_base_indices]
            accepted_data = array.array(data_typecode, accepted_data)
            subread.peer.set_tag(tag, accepted_data)
    if is_internal:
        # edit pulsecalls to make rejected pulses lowercase
        pulse_calls = fm.s2npl(subread.peer.get_tag('pc'))
        pulse_call_base_indices = fm.pls2base(pulse_calls)
        revised_pulse_calls = ''
        pulse_index = 0
        rejected_pulse_count = 0
        for base_index in pulse_call_base_indices:
            if base_index == fm.NA: # pulsecall was already rejected
                if not pulse_calls[pulse_index].islower():
                    raise RuntimeError('This pulse was already supposed to be '
                                       'rejected, but it is uppercase')
                revised_pulse_calls += pulse_calls[pulse_index].lower()
            else:
                if base_index in rejected_base_indices:
                    revised_pulse_calls += pulse_calls[pulse_index].lower()
                    rejected_pulse_count += 1
                else:
                    revised_pulse_calls += pulse_calls[pulse_index]
            pulse_index += 1
        subread.peer.set_tag('pc', revised_pulse_calls)
        # edit read sequence and quality score to remove rejected pulses
        subread.peer.seq = ''.join([base for base in subread.peer.get_tag('pc')
                                    if base.isupper()])
        subread.peer.qual = ''.join(['=' for base in subread.peer.get_tag('pc')
                                     if base.isupper()])
        deficit = rejected_pulse_count
    elif not is_internal:
        # pulsecalls info was not stored, edit seq and qual though
        seq = subread.peer.seq
        seq_dict = {} # dict structure for fast edits
        for i, base in enumerate(subread.peer.seq):
            seq_dict[i] = base
        subread.peer.seq = ''.join([seq_dict[index]
                                    for index in accepted_base_indices])
        subread.peer.qual = ''.join(['=' for index in accepted_base_indices])
        deficit = len(seq) - len(subread.peer.seq)
    # correct the query end and name
    revised_query_end = subread.qEnd - deficit
    subread.peer.set_tag('qe', revised_query_end)
    query_name = subread.peer.query_name.split('_')
    revised_query_name = query_name
    revised_query_name[-1] = str(revised_query_end)
    revised_query_name = '_'.join(revised_query_name)
    subread.peer.query_name = revised_query_name
    return subread

def reject_basecalls(read_bam_path, annotation_csv_path, out_bam_path):
    """Use index CSV to reject basecalls of subreads."""
    is_internal = is_internal_mode(read_bam_path) # check collection mode
    reject_lexicon = open_annotation_csv(annotation_csv_path)

    input_bam_reader = bu.open_input_bam(read_bam_path)
    output_bam_writer = bu.prepare_output_bam(out_bam_path,
                                              input_bam_reader)
    # iterate through subreadset and make basecall rejections
    # write to new bam on-the-fly
    subreadset = SubreadSet(read_bam_path)
    for s_index, subread in enumerate(subreadset):
        if s_index % 10 == 0: # log progress
            logging.info('Finished processing ' + str(s_index) + ' ZMWs')
        to_filter = False
        if subread.holeNumber in reject_lexicon: # reject_lexicon is dict w/ zmw as keys
            to_filter = True
            bases_to_reject = reject_lexicon[subread.holeNumber]
            # find all rows in annotation csv that reference this subread
            rows = ((subread.qStart <= bases_to_reject['annotationStartIndex']) &
                    (subread.qEnd >= bases_to_reject['annotationEndIndex']))
            reject_query_starts = bases_to_reject['annotationStartIndex'][rows]
            reject_query_ends = bases_to_reject['annotationEndIndex'][rows]
            if (not reject_query_starts.any() and
                not reject_query_ends.any()):
                    to_filter = False # subread does not need filtering

        if to_filter:
            n_reject_pulses = np.sum(reject_query_ends - reject_query_starts + 1)
            rejected_query_indices = np.zeros((n_reject_pulses,),
                                              dtype=int)
            q_ix = 0
            for start, end in itertools.izip(reject_query_starts,
                                             reject_query_ends):
                rejected_query_indices[q_ix:(q_ix + end - start + 1)] = (
                    np.arange(start, end + 1, 1))
                q_ix += end - start + 1
            rejected_query_indices.sort
            filtered_subread = filter_subread(subread,
                                              rejected_query_indices,
                                              is_internal)
            output_bam_writer.write(filtered_subread.peer)
            if filtered_subread.peer.seq is not None:
                if len(filtered_subread.peer.seq) != (filtered_subread.qEnd -
                                                      filtered_subread.qStart):
                    raise RuntimeError('Sequence must have same length '
                                       'as query region')
        if not to_filter:
            output_bam_writer.write(subread.peer)

    input_bam_reader.close()
    output_bam_writer.close()

    return bu.generate_subreadset(out_bam_path)
