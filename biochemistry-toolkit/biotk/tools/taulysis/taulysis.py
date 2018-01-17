import utils as tu
import legacy_utils as lu
import plots as pl
import numpy as np
import os
import csv
import pickle
import logging
import random
import pricompare.FastMetrics as fm
logging.basicConfig(level=logging.INFO)


def get_bam_alignment_polrates(alignment,
                               framerate,
                               coarse_grain_binsize,
                               tstart_justified,
                               pol_rates):
    """
    Retrieve polymerase rates by template-position, using
    left-justified template-positions
    """
    basecall_times = fm.aln2sec(alignment.aStart,
                                fm.s2npl(alignment.read()),
                                fm.s2npl(alignment.peer.get_tag('pc')),
                                np.array(alignment.peer.get_tag('sf'),
                                         dtype=int),
                                framerate=framerate)
    read_indices = np.flatnonzero(basecall_times > 0)
    times = basecall_times[read_indices]
    template_positions = alignment.referencePositions()[read_indices]
    template_positions = template_positions - (alignment.tStart - tstart_justified)
    tend_justified = template_positions[-1]
    # return digitized polrates by coarse_grain_bin
    first_bin = tstart_justified - (tstart_justified % coarse_grain_binsize)
    last_bin = tend_justified + (coarse_grain_binsize - tend_justified % coarse_grain_binsize)
    bins = np.arange(first_bin,
                     last_bin + coarse_grain_binsize, # include final bin
                     coarse_grain_binsize)
    bin_centers = np.arange(first_bin + coarse_grain_binsize/2,
                            last_bin + coarse_grain_binsize + coarse_grain_binsize/2, # include final bin
                            coarse_grain_binsize)
    indices = np.digitize(template_positions, bins)
    binned_polrates = {}
    for index, bin in enumerate(bin_centers):
        t = times[indices == index]
        p = template_positions[indices == index]
        if len(t) > 1:
            polrate = np.divide(p[-1] - p[0],
                                t[-1] - t[0],
                                dtype=float)
            binned_polrates[bin] = polrate
    return binned_polrates

def get_bam_alignment_timespan(alignment, framerate):
    """Retrieve start and end time of alignment for internal-mode bam"""
    # use pricompare from primary-toolkit to get alignment timespan
    basecall_times = fm.aln2sec(alignment.aStart,
                                fm.s2npl(alignment.read()),
                                fm.s2npl(alignment.peer.get_tag('pc')),
                                np.array(alignment.peer.get_tag('sf'),
                                         dtype=int),
                                framerate=framerate)
    start_time = np.divide(basecall_times[0], 60, dtype=float) # to minutes
    end_time = np.divide(basecall_times[-1], 60, dtype=float) 
    return start_time, end_time

def get_cmp_alignment_timespan(alignment, bas_reader):
    """Retrieve start and end time of alignment for legacy cmp.h5 alignment
       Requires cross-referencing against associated bas file. I will pull
       the start/end time of the HQ region and assume that's close enough.
    """
    unaligned_read = bas_reader[alignment.holeNumber]
    hqr_start_time = unaligned_read.zmwMetric('HQRegionStartTime')
    hqr_end_time = unaligned_read.zmwMetric('HQRegionEndTime')
    # convert to minutes. stored units are seconds
    hqr_start_time = np.divide(hqr_start_time, 60, dtype=float)
    hqr_end_time = np.divide(hqr_end_time, 60, dtype=float)
    return hqr_start_time, hqr_end_time

def left_justify_template(template_start, template_end, smrtbellsize):
    """Use SMRTBell size to left justify the template start and end"""
    smrtbellsize = smrtbellsize
    shift = template_start / smrtbellsize
    justified_template_start = template_start - smrtbellsize * shift
    justified_template_end = template_end - smrtbellsize * shift
    return (justified_template_start,
            justified_template_end)

def log_movie_limit_determination_rule(threshold, reference_name):
    """If movie-limiting determination is happening, log details"""
    logging.info(('Time info is available. Movie-limited flagging active'))
    logging.info(('For alignments mapping to reference ' +
                  reference_name + ', all alignments that ' +
                  'terminate within ' + str(threshold) +
                  ' minute(s) of the end of the movie are flagged as ' +
                  'movie-limited.'))

def get_template_positions(alignment_set,
                           alignment_ixs,
                           bas_reader,
                           reference_info,
                           movie_length,
                           movie_limited_threshold,
                           template_min_start,
                           template_max_start,
                           max_start_time,
                           framerate,
                           coarse_grain_binsize,
                           is_legacy):
    """For a list of alignments, retrieve tStarts, tEnds, and whether
       the alignment was movie-limited (if flagged). If not flagged
       no alignments are presumed to be movie-limited. This code is a
       mess of if/else statements because I have to support the legacy
       cmp.h5 file format. Without that, the code would be quite simple.
    """
    
    # subsample alignments to speed up measurement
    subsample_to = 10000
    if len(alignment_ixs) > subsample_to:
        logging.info(('Found ' + str(len(alignment_ixs)) +
                      ' alignments on reference ' +
                      str(reference_info['FullName'])))
        logging.info('Since more than ' + str(subsample_to) + ' alignments found, subsampling.')
        alignment_ixs = np.array(random.sample(alignment_ixs, subsample_to),
                                 dtype=int)

    smrtbellsize = reference_info['SMRTBellSize']
    summary_data = np.recarray((len(alignment_ixs), ),
                               dtype=[('holeNumber', int),
                                      ('tStart', int),
                                      ('tEnd', int),
                                      ('tStartJustified', int),
                                      ('tEndJustified', int),
                                      ('startTime', float),
                                      ('endTime', float),
                                      ('isMovieLimited', bool)])
    pol_rates = {}

    movie_limit_logged = False
    for index, alignment_ix in enumerate(alignment_ixs): # iterate for time information
        # assume alignment is not movie-limited until otherwise proven
        alignment = alignment_set[alignment_ix]
        is_movie_limited = False
        # get time information bam format
        if not is_legacy and movie_length is not None:
            start_time, end_time = get_bam_alignment_timespan(alignment,
                                                              framerate)
            if not movie_limit_logged:
                log_movie_limit_determination_rule(movie_limited_threshold,
                                                   str(reference_info['FullName']))
                movie_limit_logged = True
            if float(movie_length) - end_time < movie_limited_threshold:
                is_movie_limited = True
        elif not is_legacy:
            start_time, end_time = (np.nan, np.nan)
            is_movie_limited = False
        # get time information cmph5 legacy
        elif is_legacy and bas_reader and movie_length:
            start_time, end_time = get_cmp_alignment_timespan(alignment,
                                                              bas_reader)
            if not movie_limit_logged:
                log_movie_limit_determination_rule(movie_limited_threshold,
                                                   str(reference_info['FullName']))
                movie_limit_logged = True
            if float(movie_length) - end_time < movie_limited_threshold:
                is_movie_limited = True
        else:
            start_time, end_time = (np.nan, np.nan)
            is_movie_limited = False

        # left justify template positions
        (justified_template_start,
         justified_template_end) = left_justify_template(alignment.tStart,
                                                         alignment.tEnd,
                                                         smrtbellsize)
        summary_data['holeNumber'][index] = alignment.holeNumber
        summary_data['tStart'][index] = alignment.tStart
        summary_data['tEnd'][index] = alignment.tEnd
        summary_data['tStartJustified'][index] = justified_template_start
        summary_data['tEndJustified'][index] = justified_template_end
        summary_data['startTime'][index] = start_time
        summary_data['endTime'][index] = end_time
        summary_data['isMovieLimited'][index] = is_movie_limited

        # pipe in the polymerase rate logic
        if not is_legacy and movie_length is not None:
            # data are internal-mode and sequel-format
            # calculate polymerase rates
            if (justified_template_start > int(template_min_start)) & (justified_template_start < int(template_max_start)) & (start_time < int(max_start_time)):
                aln_pol_rate = get_bam_alignment_polrates(alignment,
                                                          framerate,
                                                          coarse_grain_binsize,
                                                          justified_template_start,
                                                          pol_rates)
                for tpos in aln_pol_rate:
                    if tpos not in pol_rates:
                        pol_rates[tpos] = [aln_pol_rate[tpos], 1] # store [running mean, count]
                    else:
                        N = pol_rates[tpos][1] + 1
                        prev_mean = pol_rates[tpos][0]
                        new_value = aln_pol_rate[tpos]
                        pol_rates[tpos][0] = prev_mean + np.divide(new_value - prev_mean,
                                                                   N,
                                                                   dtype=float)
                        pol_rates[tpos][1] = N
                    
    return summary_data, pol_rates

def clean_template_positions(termination_info,
                             reference_info,
                             template_min_start,
                             template_max_start,
                             maximum_start_time): # minutes
    """Use template start min/max and start time max to filter for suitable
       alignments for survival-by-template-position analysis.
    """
    # must start within x bases of beginning
    reference_name = reference_info['FullName']
    logging.info(('For alignments mapping to reference ' + str(reference_name) +
                  ', only alignments that start after base position ' +
                  str(template_min_start) + ' and before base position ' +
                  str(template_max_start) + ' and within ' +
                  str(maximum_start_time) + ' minutes of acquistion start ' +
                  'are included for analysis. ' +
                  'Time-based filtering ignored when time info not available.'))
    start_times = termination_info['startTime']
    if np.isnan(start_times).all(): # time info not stored
        start_times[:] = 0.
    pass_filter = ((termination_info['tStartJustified'] >
                    int(template_min_start)) &
                   (termination_info['tStartJustified'] <
                    int(template_max_start)) &
                   (start_times <
                    maximum_start_time))
    terminations_filtered = termination_info[['tStart',
                                              'tEnd',
                                              'tStartJustified',
                                              'tEndJustified',
                                              'startTime',
                                              'isMovieLimited']][pass_filter]

    return terminations_filtered

def calculate_survival(termination_info):
    """Calculate survival curve by template position"""
    not_movie_limited = termination_info['isMovieLimited'] == False
    template_starts = termination_info['tStartJustified'][not_movie_limited]
    template_ends = termination_info['tEndJustified'][not_movie_limited]
    template_starts, template_ends = zip(*np.sort(zip(template_starts,
                                                      template_ends),
                                                  axis=0))
    number_alignments = len(termination_info) # all alignments for total count
    survivor_count = np.array([number_alignments for row in template_ends],
                              dtype=int) # initialize array to tally survivors
    for index in np.arange(len(template_ends)):
        survivor_count[index::] -= 1
    # normalize
    survivor_count = np.divide(survivor_count, number_alignments, dtype=float)
    return template_starts, template_ends, survivor_count

def coarse_grain_survival(template_position, survival, coarse_grain_binsize):
    """Coarse grain the survival curve given a window size"""
    template_position = np.array(template_position)
    survival = np.array(survival)
    start = template_position[0]
    end = start + coarse_grain_binsize
    nbins = np.divide(template_position[-1] - template_position[0],
                      coarse_grain_binsize,
                      dtype=int)
    cg_template_position = np.zeros((nbins, ), dtype=int)
    cg_survival = np.zeros((nbins, ), dtype=float)
    for index in np.arange(len(cg_template_position)):
        in_slice = ((template_position >= start) &
                    (template_position < end))
        if in_slice.any():
            cg_template_position[index] = np.mean(template_position[in_slice])
            cg_survival[index] = np.mean(survival[in_slice])
        else: # no data found, interpolate template position and use previous survival value
            cg_template_position[index] = (cg_template_position[index-1] +
                                           coarse_grain_binsize / 2)
            cg_survival[index] = cg_survival[index-1]
        start = end
        end = end + coarse_grain_binsize
    return cg_template_position, cg_survival

def calculate_termination_rate(template_position,
                               survival,
                               coarse_grain_binsize):
    """Calculate termination rate from survival curve"""
    # coarse-grain survival curve for calculating derivative
    (cg_template_position,
     cg_survival) = coarse_grain_survival(template_position,
                                          survival,
                                          coarse_grain_binsize)
    # estimate derivative of survival curve using coarse-grained data
    dx = np.diff(cg_template_position)
    dy = np.diff(cg_survival)
    np.seterr(divide='warn', invalid='warn') # allow, but log troubled division
    dydx = np.divide(dy, dx, dtype=float)
    dydx_survival_scaled = np.divide(dydx,
                                     cg_survival[1:],
                                     dtype=float)
    termination_rate = -dydx_survival_scaled
    return cg_template_position[1:], termination_rate 

def fit_taus(template_position,
             termination_rate,
             first_adapter_hit,
             second_adapter_hit,
             movie_length):
    """Fit termination rate to get Tau1, Tau2, and TauRC"""
    logging.info(('First adapter hit expected to be at template position ' +
                  str(first_adapter_hit) + '. Second adapter hit expected '
                  'to be at template position ' + str(second_adapter_hit) + '.'
                  ' If these values are incorrect, rerun with flags'
                  '--first-adapter-hit and --second-adapter-hit enabled.'))
    adapter_buffer = 333 # data exclusion buffer around adapters
    if movie_length:
        # fit whole rolling-circle region if movie-limited alignments flagged
        rolling_circle_max = np.max(template_position)
    else:
        # fit one-insert length to avoid movie-limited region
        logging.info(('Movie-limited traces not flagged, fitting one insert\n'
                      'length past second adapter to estimate rolling-circle\n'
                      'tau.'))
        rolling_circle_max = (second_adapter_hit +
                              second_adapter_hit - first_adapter_hit)

    first_pass_filter = template_position < first_adapter_hit - adapter_buffer
    second_pass_filter = ((template_position >
                           first_adapter_hit + adapter_buffer) &
                          (template_position <
                           second_adapter_hit - adapter_buffer))
    rolling_circle_filter = ((template_position >
                              second_adapter_hit + adapter_buffer) &
                             (template_position < rolling_circle_max))
    first_pass_terminations = termination_rate[first_pass_filter]
    second_pass_terminations = termination_rate[second_pass_filter]
    rolling_circle_terminations = termination_rate[rolling_circle_filter]

    first_pass_terminations = first_pass_terminations[first_pass_terminations != 0]
    second_pass_terminations = second_pass_terminations[second_pass_terminations != 0]
    rolling_circle_terminations = rolling_circle_terminations[rolling_circle_terminations != 0]

    first_pass_taus = np.divide(1., first_pass_terminations, dtype=float)
    second_pass_taus = np.divide(1., second_pass_terminations, dtype=float)
    rolling_circle_taus = np.divide(1., rolling_circle_terminations, dtype=float)



    tau_1 = np.nanmedian(first_pass_taus)
    tau_2 = np.nanmedian(second_pass_taus)
    tau_rc = np.nanmedian(rolling_circle_taus)
    # catch divide-by-zero edge case
    # occurence will be logged by numpy
    if tau_1 == np.infty or np.isnan(tau_1):
        tau_1 = 1
    if tau_2 == np.infty or np.isnan(tau_2):
        tau_2 = 1
    if tau_rc == np.infty or np.isnan(tau_rc):
        tau_rc = 1

    return (int(tau_1), # use mean to protect from divide-by-zero overestimation
            int(tau_2), # use median to insulate from large skew
            int(tau_rc))

def save_taus(alignment_set,
              alignment_indices,
              cohort_information,
              reference_info,
              taus,
              results_file,
              fieldnames):
    """
    Save results.
    """
    results = {
        'referenceName': str(reference_info['FullName']),
        'nTotalAlignments': len(alignment_indices),
        'nCohortAlignments': len(cohort_information),
        'percentMovieLimited': np.divide(
            len(cohort_information[cohort_information['isMovieLimited'] == True]),
            len(cohort_information),
            dtype=float),
        'tau1': taus[0],
        'tau2': taus[1],
        'tauRC': taus[2]
    }
    with open(results_file, 'a') as file:
        appender = csv.DictWriter(file, fieldnames=fieldnames)
        appender.writerow(results)

def get_taus(alignment_set,
             bas_reader,
             reference_information,
             output_directory,
             movie_length,
             movie_limited_threshold,
             first_adapter_template_position,
             second_adapter_template_position,
             template_min_start,
             template_max_start,
             coarse_grain_binsize,
             is_legacy):
    """Perform tau analysis, by reference
       Generate the following plots:
            1) template start histogram
            2) start time histogram
            3) template-start position by start-time
            4) survival by template-position
            5) termination rate by template-position
            6) termination rates w/ 1st, 2nd, and RC-circle tau fits
       Save a summary CSV storing tau values, by reference
    """
    framerate = alignment_set.resourceReaders()[0].readGroupTable.FrameRate[0]

    # set up results file by writing header
    results_file = os.path.join(output_directory, 'results_summary.csv')
    with open(results_file, 'wb') as file:
        fieldnames = ['referenceName',
                      'nTotalAlignments',
                      'nCohortAlignments',
                      'percentMovieLimited',
                      'tau1',
                      'tau2',
                      'tauRC']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        
    for ref in reference_information:
        logging.info(('Processing alignments from reference ' + 
                      str(ref['FullName'])))
        # row in recarray w/ 'ID', 'FullName', and 'SMRTBellSize'
        if is_legacy:
            aln_to_ref_indices = np.flatnonzero(alignment_set.index['RefGroupID'] ==
                                                int(ref['ID']))
        else:
            aln_to_ref_indices = np.flatnonzero(alignment_set.index['tId'] ==
                                                int(ref['ID']))

        if len(aln_to_ref_indices) > 0: # check that there are alignments
            max_start_time = 5
            (template_position_info,
             pol_rates) = get_template_positions(alignment_set,
                                                 aln_to_ref_indices,
                                                 bas_reader,
                                                 ref,
                                                 movie_length,
                                                 float(movie_limited_threshold),
                                                 template_min_start,
                                                 template_max_start,
                                                 max_start_time,
                                                 framerate,
                                                 coarse_grain_binsize,
                                                 is_legacy)
            pickle.dump(template_position_info, open((output_directory +
                                                      'info_' +
                                                      str(ref['FullName']) +
                                                      '.pkl'), 'wb'))
            pickle.dump(pol_rates, open((output_directory +
                                         'rates_' +
                                         str(ref['FullName']) +
                                         '.pkl'), 'wb'))
            # remove alignments that started late
            clean_termination_info = clean_template_positions(template_position_info,
                                                              ref,
                                                              template_min_start,
                                                              template_max_start,
                                                              max_start_time)
            # generate template start histogram
            plot_descriptors = ['templateStartHistogram']
            pl.plot_histogram(data=template_position_info['tStart'],
                              output_directory=output_directory,
                              title=ref['FullName'],
                              xlabel='Template start (bp)',
                              ylabel='Counts',
                              bins=[0,
                                    np.max(template_position_info['tStart']),
                                    10], # bases binwidth
                              descriptors=plot_descriptors)
            plot_descriptors = ['templateStartJustifiedHistogram']
            pl.plot_histogram(data=template_position_info['tStartJustified'],
                              output_directory=output_directory,
                              title=ref['FullName'],
                              xlabel='Template start (bp)',
                              ylabel='Counts',
                              bins=[0,
                                    np.max(template_position_info['tStartJustified']),
                                    10], # bases binwidth
                              descriptors=plot_descriptors)

            if movie_length:
                # generate start time histogram
                plot_descriptors = ['startTimeHistogram']
                pl.plot_histogram(data=template_position_info['startTime'],
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  xlabel='Start time (min)',
                                  ylabel='Counts',
                                  bins=[0,
                                        movie_length,
                                        1], # minute binwidth
                                  descriptors=plot_descriptors)
                # generate tStart vs. start time plots, if time info exists
                plot_descriptors = ['tStart_vs_startTime_unjustified']
                pl.plot_line_plot(x=template_position_info['startTime'],
                                  y=template_position_info['tStart'],
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  datapoint_mode='markers',
                                  xlabel='Start time (min)',
                                  ylabel='Template start position (bp)',
                                  yaxis_range=None,
                                  descriptors=plot_descriptors)
                plot_descriptors = ['tStart_vs_startTime_justified']
                pl.plot_line_plot(x=template_position_info['startTime'],
                                  y=template_position_info['tStartJustified'],
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  datapoint_mode='markers',
                                  xlabel='Start time (min)',
                                  ylabel='Template start position (bp)',
                                  yaxis_range=None,
                                  descriptors=plot_descriptors)
                
                if not is_legacy and pol_rates:
                    # plot pol rates
                    xyz = []
                    for tpos in pol_rates:
                        rate, count = pol_rates[tpos]
                        xyz.append((tpos, rate, count))
                    xyz.sort(key=lambda tup: tup[0]) # sort for line plot
                    tpos, rate, count = zip(*xyz)
                    plot_descriptors = ['polrate_by_template_position',
                                        'tStartRange',
                                        str(template_min_start),
                                        str(template_max_start)]
                    pl.plot_line_plot(x=tpos,
                                      y=rate,
                                      output_directory=output_directory,
                                      title=ref['FullName'],
                                      datapoint_mode='lines',
                                      xlabel='Template position (bp)',
                                      ylabel='Mean pol rate (bases/sec)',
                                      yaxis_range=[0, 4],
                                      descriptors=plot_descriptors)
                    plot_descriptors = ['numalns_by_template_position',
                                        'tStartRange',
                                        str(template_min_start),
                                        str(template_max_start)]
                    pl.plot_line_plot(x=tpos,
                                      y=count,
                                      output_directory=output_directory,
                                      title=ref['FullName'],
                                      datapoint_mode='lines',
                                      xlabel='Template position (bp)',
                                      ylabel='Number of alignments',
                                      yaxis_range=[0, len(aln_to_ref_indices)],
                                      descriptors=plot_descriptors)

            minimum_alignment_number = 500
            if len(clean_termination_info) > minimum_alignment_number:
                # produce survival curve
                (template_start_position,
                 template_end_position,
                 survival) = calculate_survival(clean_termination_info)
                # plot survival curve
                plot_descriptors = ['survival_by_template_position',
                                    'tStartRange',
                                    str(template_min_start),
                                    str(template_max_start)]
                pl.plot_line_plot(x=template_end_position,
                                  y=survival,
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  datapoint_mode='lines',
                                  xlabel='Template position (bp)',
                                  ylabel='Fraction survivors left',
                                  yaxis_range=[0, 1],
                                  descriptors=plot_descriptors)
                # produce termination rates
                (template_position,
                 termination_rate) = calculate_termination_rate(template_end_position,
                                                                survival,
                                                                coarse_grain_binsize)
                # plot termination rates
                plot_descriptors = ['termination_rates',
                                    'tStartRange',
                                    str(template_min_start),
                                    str(template_max_start)]
                pl.plot_line_plot(x=template_position,
                                  y=termination_rate,
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  datapoint_mode='lines',
                                  xlabel='Template position (bp)',
                                  ylabel='Termination rate (per bp)',
                                  yaxis_range=None,
                                  descriptors=plot_descriptors)
                # plot taus
                plot_descriptors = ['termination_taus',
                                    'tStartRange',
                                    str(template_min_start),
                                    str(template_max_start)]
                pl.plot_line_plot(x=template_position[termination_rate != 0],
                                  y=np.divide(1,
                                              termination_rate[termination_rate != 0],
                                              dtype=float),
                                  output_directory=output_directory,
                                  title=ref['FullName'],
                                  datapoint_mode='lines',
                                  xlabel='Template position (bp)',
                                  ylabel='Tau (bp)',
                                  yaxis_range=None,
                                  descriptors=plot_descriptors)
                # fit termination rates to get taus
                if first_adapter_template_position is None:
                    first_adapter_template_position = np.divide(ref['SMRTBellSize'],
                                                                2,
                                                                dtype=int)
                if second_adapter_template_position is None:
                    second_adapter_template_position = ref['SMRTBellSize']
                tau_1, tau_2, tau_rc = fit_taus(template_position,
                                                termination_rate,
                                                int(first_adapter_template_position),
                                                int(second_adapter_template_position),
                                                movie_length)
                logging.info(('\n\nTau1 is ' + str(tau_1) + ' bases\n'
                              'Tau2 is ' + str(tau_2) + ' bases\n'
                              'TauRC is ' + str(tau_rc) + ' bases\n'))
                # plot termination rate with fits
                plot_descriptors = ['fitted_termination_rates']
                pl.plot_fitted_line_plot(x=template_position,
                                         y=termination_rate,
                                         t_1=[np.min(template_position),
                                              first_adapter_template_position,
                                              np.divide(1., tau_1, dtype=float)],
                                         t_2=[first_adapter_template_position,
                                              second_adapter_template_position,
                                              np.divide(1., tau_2, dtype=float)],
                                         t_rc=[second_adapter_template_position,
                                               np.max(template_position),
                                               np.divide(1., tau_rc, dtype=float)],
                                         output_directory=output_directory,
                                         title=ref['FullName'],
                                         xlabel='Template position (bp)',
                                         ylabel='Termination rate (per bp)',
                                         descriptors=plot_descriptors)
                # plot termination taus with fits
                plot_descriptors = ['fitted_termination_taus']
                pl.plot_fitted_line_plot(x=template_position[termination_rate != 0],
                                         y=np.divide(1.,
                                                     termination_rate[
                                                        termination_rate != 0],
                                                     dtype=float),
                                         t_1=[np.min(template_position),
                                              first_adapter_template_position,
                                              tau_1],
                                         t_2=[first_adapter_template_position,
                                              second_adapter_template_position,
                                              tau_2],
                                         t_rc=[second_adapter_template_position,
                                               np.max(template_position),
                                               tau_rc],
                                         output_directory=output_directory,
                                         title=ref['FullName'],
                                         xlabel='Template position (bp)',
                                         ylabel='Termination tau (bp)',
                                         descriptors=plot_descriptors)
                taus = (tau_1, tau_2, tau_rc)
                save_taus(alignment_set,
                          aln_to_ref_indices,
                          clean_termination_info,
                          ref,
                          taus,
                          results_file,
                          fieldnames)
            else:
                logging.info(('Less than ' + str(minimum_alignment_number) +
                              ' alignments exist against reference ' +
                              str(ref['FullName']) + '. Analysis skipped.'))
        else:
            logging.info(('No alignments produced on reference ' +
                          str(ref['FullName']) + '.'))

    return None
