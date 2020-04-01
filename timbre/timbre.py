
import numpy as np
from urllib.request import urlopen
from urllib.error import URLError
import json
from os.path import expanduser
from scipy import interpolate
from hashlib import md5

from Chandra.Time import DateTime
import xija

home = expanduser("~")

results_dtype = [('msid', '|U20'),
                 ('date', '|U8'),
                 ('datesecs', np.float64),
                 ('obsid1', np.int32),
                 ('sequence1', np.int32),
                 ('duration1_fraction', np.float64),
                 ('obsid2', np.int32),
                 ('sequence2', np.int32),
                 ('limit', np.float64),
                 ('pitch1', np.float64),
                 ('pitch2', np.float64),
                 ('roll1', np.float64),
                 ('roll2', np.float64),
                 ('ccd_count1', np.int8),
                 ('ccd_count2', np.int8),
                 ('fep_count1', np.int8),
                 ('fep_count2', np.int8),
                 ('clocking1', np.bool),
                 ('clocking2', np.bool),
                 ('vid_board1', np.bool),
                 ('vid_board2', np.bool),
                 ('sim_z1', np.int32),
                 ('sim_z2', np.int32),
                 ('t_dwell1', np.float64),
                 ('t_dwell2', np.float64),
                 ('min_temp', np.float64),
                 ('mean_temp', np.float64),
                 ('max_temp', np.float64),
                 ('min_pseudo', np.float64),
                 ('mean_pseudo', np.float64),
                 ('max_pseudo', np.float64),
                 ('converged', np.bool),
                 ('unconverged_hot', np.bool),
                 ('unconverged_cold', np.bool),
                 ('hotter_state', np.int8),
                 ('colder_state', np.int8)]


pseudo_names = dict(
    zip(['aacccdpt', 'pftank2t', '1dpamzt', '4rt700t', '1deamzt'], ['aca0', 'pf0tank2t', 'dpa0', 'oba0', None]))


def load_model_specs():
    """ Load Xija model parameters for all available models.

    Returns:
        dictionary: A dictionary containing the model specifications for all available Xija models

    Note:
        This will need to be updated as new models are approved or existing models are renamed.

    """

    def get_model(branch, internet):
        """ Load parameters for a single Xija model.

        Args:
            branch (str): Relative location of model file, starting from the chandra_models/chandra_models/xija/
                directory
            internet (bool): Availability of an internet connection, for accessing github.com

        Returns:
            dictionary: JSON file stored as a dictionary, containing Xija model parameters

        """

        url = 'https://raw.githubusercontent.com/sot/chandra_models/master/chandra_models/xija/'
        local_dir = '/AXAFLIB/chandra_models/chandra_models/xija/'

        if internet:
            model_spec_url = url + branch  # aca/aca_spec.json'
            with urlopen(model_spec_url) as url:
                response = url.read()
                f = response.decode('utf-8')
        else:
            with open(home + local_dir + branch) as fid:  # 'aca/aca_spec.json', 'rb') as fid:
                f = fid.read()

        md5_hash = md5(f.encode('utf-8')).hexdigest()

        return json.loads(f), md5_hash


    model_specs = {}

    internet = True
    try:
        _ = urlopen('https://github.com')
    except URLError:
        internet = False

    model_specs['aacccdpt'], model_specs['aacccdpt_hash'] = get_model('aca/aca_spec.json', internet)
    model_specs['1deamzt'], model_specs['1deamzt_hash'] = get_model('dea/dea_spec.json', internet)
    model_specs['1dpamzt'], model_specs['1dpamzt_hash'] = get_model('dpa/dpa_spec.json', internet)
    model_specs['fptemp'], model_specs['fptemp_hash'] = get_model('acisfp/acisfp_spec.json', internet)
    model_specs['1pdeaat'], model_specs['1pdeaat_hash'] = get_model('psmc/psmc_spec.json', internet)
    model_specs['pftank2t'], model_specs['pftank2t_hash'] = get_model('pftank2t/pftank2t_spec.json', internet)
    model_specs['tcylaft6'], model_specs['tcylaft6_hash'] = get_model('tcylaft6/tcylaft6_spec.json', internet)
    model_specs['4rt700t'], model_specs['4rt700t_hash'] = get_model('fwdblkhd/4rt700t_spec.json', internet)
    model_specs['pline03t'], model_specs['pline03t_hash'] = get_model('pline/pline03t_model_spec.json', internet)
    model_specs['pline04t'], model_specs['pline04t_hash'] = get_model('pline/pline04t_model_spec.json', internet)
    model_specs['pm1thv2t'], model_specs['pm1thv2t_hash'] = get_model('mups_valve/pm1thv2t_spec.json', internet)
    model_specs['pm2thv1t'], model_specs['pm2thv1t_hash'] = get_model('mups_valve/pm2thv1t_spec.json', internet)

    return model_specs


def c_to_f(temp):
    """ Convert Celsius to Fahrenheit

    Args:
        temp (int, float, numpy.ndarray, list, tuple): Temperature in Celsius

    Returns:
        (int, float, numpy.ndarray, list, tuple): Temperature in Fahrenheit

    """
    if type(temp) is list or type(temp) is tuple:
        return [c * 1.8 + 32 for c in temp]
    else:
        return temp * 1.8 + 32.0


def f_to_c(temp):
    """ Convert Fahrenheit to Celsius

    Args:
        temp (int, float, numpy.ndarray, list, tuple): Temperature in Fahrenheit

    Returns:
        (int, float, numpy.ndarray, list, tuple): Temperature in Celsius

    """
    if type(temp) is list or type(temp) is tuple:
        return [(c - 32) / 1.8 for c in temp]
    else:
        return (temp - 32.0) / 1.8


def setup_model(msid, t0, t1, model_spec, init):
    """ Create Xija model object

    This function creates a Xija model object with initial parameters, if any. This function is intended to create a
    streamlined method to creating Xija models that can take both single value data and time defined data
    (e.g. [pitch1, pitch2, pitch3], [time1, time2, time3]), defined in the `init` dictionary.

    Args:
        msid (str): Primary MSID for model; in this case it can be anything as it is only being used to name the model,
            however keeping the convention to name the model after the primary MSID being predicted reduces confusion
        t0 (str, float, int): Start time for model prediction; this can be any format that Chandra.Time.DateTime accepts
        t1 (str, float, int): End time for model prediction; this can be any format that Chandra.Time.DateTime accepts
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters, can be empty

    Returns:
        (xija.model.XijaModel): Xija model object

    Example:

        model_specs = load_model_specs()

        init = {'1dpamzt': 35., 'dpa0': 35., 'eclipse': False, 'roll': 0, 'vid_board': True, 'pitch':155,
                          'clocking': True, 'fep_count': 5, 'ccd_count': 5, 'sim_z': 100000}

        model = setup_model('1dpamzt', '2019:001:00:00:00', '2019:010:00:00:00', model_specs['1dpamzt'], init)

    Notes:
        This does not run the model, only sets up the model to be run.

        Any parameters not specified in `init` will either need to be pulled from telemetry or explicitly defined
        outside of this function before running the model.

    """

    model = xija.ThermalModel(msid, start=t0, stop=t1, model_spec=model_spec)
    for key, value in init.items():
        if isinstance(value, dict):
            model.comp[key].set_data(value['data'], value['times'])
        else:
            model.comp[key].set_data(value)

    return model


def run_profile(times, schedule, msid, model_spec, init, pseudo=None):
    """ Run a Xija model for a given time and state profile.

    Args:
        times (numpy.ndarray): Array of time values, in seconds from '1997:365:23:58:56.816' (Chandra.Time.DateTime
            epoch)
        schedule (dict): Dictionary of pitch, roll, etc. values that match the time values specified above in `times`
        msid (str): Primary MSID for model being run
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters, can be empty but not recommended
        pseudo (:obj:`str`, optional): Name of one or more pseudo MSIDs used in the model, if any, only necessary if one
            wishes to retrieve model results for this pseudo node, if it exists

    Returns:
        dict: Dictionary of results, keys are node names (e.g. 'aacccdpt', 'aca0'), values are Xija model component
            objects

    Example:

        times = np.array(DateTime(['2019:001:00:00:00', '2019:001:12:00:00', '2019:002:00:00:00',
                                   '2019:003:00:00:00']).secs)

        pitch = np.array([150, 90, 156, 156])

        schedule = {'pitch': pitch}

        model_specs = load_model_specs()

        init = {'1dpamzt': 20., 'dpa0': 20., 'eclipse': False, 'roll': 0, 'vid_board': True,
                          'clocking': True, 'fep_count': 5, 'ccd_count': 5, 'sim_z': 100000}

        results = run_profile(times, pitch, '1dpamzt', model_specs['1dpamzt'], init, pseudo='dpa0')

    Notes:
        Any parameters specified in `init` will be overwritten by those specified in the body of this function, if they
        happen to be defined in both places.

    """

    model = setup_model(msid, times[0], times[-1], model_spec, init)

    for key, value in schedule.items():
        model.comp[key].set_data(value, times=times)

    model.make()
    model.calc()
    tmsid = model.get_comp(msid)
    results = {msid: tmsid}

    if pseudo is not None:
        results[pseudo] = model.get_comp(pseudo_names[msid])

    return results


def calc_binary_schedule(datesecs, state1, state2, t_dwell1, t_dwell2, msid, model_spec, init, duration=2592000.,
                         t_backoff=1725000., pseudo=None):
    """ Simulate a schedule that switches between two states

    This runs the model over a "binary" schedule. This function is intended to be used to optimize the `t_dwell2`
    parameter so that the predicted temperature during the last `t_backoff` number of seconds peaks within a tolerance
    of a limit (limit used and specified in a different function).

    Args:
        datesecs (float, int): Date for start of simulation, in seconds from '1997:365:23:58:56.816'
            (Chandra.Time.DateTime epoch)
        state1 (dict): States for fixed dwell (pitch, roll, ccds, etc.)
        state2 (dict): States for variable dwell (pitch, roll, ccds, etc.)
        t_dwell1 (float, int): Fixed dwell duration in seconds, this is in the SCALED format
        t_dwell2 (list, tuple): Variable dwell duration in seconds (this is the parameter that is optimized), the
            optimization routine returns a list, though in this case we are only interested in the first value (only
            value?), this is in the SCALED format
        msid (str): Primary MSID for model being run
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters
        duration (:obj:`float`, optional): Duration for entire simulated schedule, defaults to 30 days (in seconds)
        t_backoff (:obj:`float`, optional): Duration for tail end of simulated schedule used to determine convergence,
            defaults to 10 days (in seconds)
        pseudo (:obj:`str`, optional): Name of one or more pseudo MSIDs used in the model, if any, only necessary if one
            wishes to retrieve model results for this pseudo node, if it exists

    Returns:
        dict: Dictionary of results, keys are node names (e.g. 'aacccdpt', 'aca0'), values are Xija model component
            objects, this is the same object returned by `run_profile`
        ndarray: Numpy array of time values input into Xija (may not exactly match Xija output)
        ndarray: Numpy array of state markers, matching the time array output (may not exactly match Xija output)

    Notes:
        Keys in state1 must match keys in state2.
        Keys in state1 must match Xija component names (e.g. 'pitch', 'ccd_count', 'sim_z')

    """

    num = np.int(duration / (t_dwell1 + t_dwell2))
    reltimes = np.cumsum([1, t_dwell1 - 1, 1, t_dwell2 - 1] * num)
    times = np.array(reltimes) - reltimes[0] + datesecs - t_backoff

    schedule = dict(zip(state1.keys(), []))
    for key, value in state1.items():
        layout = [state1[key], state1[key], state2[key], state2[key]] * num
        schedule[key] = np.array(layout)

    statekey = [1, 1, 2, 2] * num
    statekey = np.array(statekey)

    model_results = run_profile(times, schedule, msid, model_spec, init, pseudo=pseudo)

    return model_results, times, statekey


def create_opt_fun(datesecs, dwell1_state, dwell2_state, t_dwell1, msid, model_spec, init, t_backoff,
                   duration):
    """ Generate a Xija model function with preset values, for use with an optimization routine.

    Args:
        datesecs (float, int): Date for start of simulation, in seconds from '1997:365:23:58:56.816'
            (Chandra.Time.DateTime epoch)
        dwell1_state (dict): States for fixed dwell (pitch, roll, ccds, etc.)
        dwell2_state (dict): States for variable dwell (pitch, roll, ccds, etc.)
        t_dwell1 (float, int): Fixed dwell duration in seconds, this is in the SCALED format
        msid: msid (str): Primary MSID for model being run
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters, can be empty
        t_backoff (float): Duration for tail end of simulated schedule, used to determine convergence
        duration (float): Duration for entire simulated schedule, defaults to 30 days (in seconds)

    Returns:
        function: Function generated from specified parameters, to be passed to optimization routine

    Notes:

        Keys in dwell1_state must match keys in dwell2_state.

        Keys in dwell1_state must match Xija component names (e.g. 'pitch', 'ccd_count', 'sim_z')

    """
    def opt_binary_schedule(t):
        model_results, _, _ = calc_binary_schedule(datesecs, dwell1_state, dwell2_state, t_dwell1, t, msid,
                                                   model_spec, init, duration=duration, t_backoff=t_backoff)

        model_temps = model_results[msid].mvals
        model_times = model_results[msid].times
        ind = model_times > (model_times[-1] - t_backoff)
        dmax = np.max(model_temps[ind])
        dmin = np.min(model_temps[ind])
        dmean = np.mean(model_temps[ind])

        return t, dmax, dmean, dmin

    return opt_binary_schedule


def find_second_dwell(date, dwell1_state, dwell2_state, t_dwell1, msid, limit, model_spec, init,
                      duration=2592000, t_backoff=1725000, n_dwells=10., max_dwell=None, pseudo=None):
    """ Determine the required dwell time at pitch2 to balance a given fixed dwell time at pitch1, if any exists.

    Args:
        date (float, int, str): Date for start of simulation, in seconds from '1997:365:23:58:56.816', or any other
            format readable by Chandra.Time.DateTime
        dwell1_state (dict): States for fixed dwell (pitch, roll, ccds, etc.)
        dwell2_state (dict): States for variable dwell (pitch, roll, ccds, etc.)
        t_dwell1 (float, int): Fixed dwell duration in seconds
        msid: msid (str): Primary MSID for model being run
        limit (float): Temperature limit for primary MSID in model for this simulation
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters, can be empty
        duration (float): Duration for entire simulated schedule, defaults to 30 days (in seconds)
        t_backoff (float): Duration for tail end of simulated schedule used to determine convergence, defaults to 10
            days (in seconds)
        n_dwells (int): Number of second dwell possibilities to run (more dwells = finer resolution)
        max_dwell (float): Maximum duration for second dwell, can be tuned to provide better results
        pseudo (:obj:`str`, optional): Name of one or more pseudo MSIDs used in the model, if any, only necessary if one
            wishes to retrieve model results for this pseudo node, if it exists - To be implemented at a later date

    Returns:
        dict: Dictionary of results information
        ndarray: Numpy array of maximum, mean, and minimum temperatures for each simulation generated, within the last
            `t_backoff` duration (e.g. the last two thirds of `duration`).

    """

    # def opt_dwell2():
    #     t_vals = (1.0e-6, max_dwell / 2, max_dwell)
    #
    #     output = np.array([opt_fun(t) for t in t_vals], dtype=[('duration2', np.float64), ('max', np.float64),
    #                                                            ('mean', np.float64), ('min', np.float64)])
    #     output_sorted = np.sort(output, order='max')  # low to high
    #     ind = np.searchsorted(output_sorted['max'], limit)
    #
    #     n = 0
    #     while ((np.max(t_vals) - np.min(t_vals)) > 10) and (n < 20):
    #         n = n + 1
    #         if ind == 0:
    #             return output_sorted[0]
    #
    #         elif ind == 1:
    #             t_vals = [output_sorted['duration2'][0], np.mean(output_sorted['duration2'][:2]),
    #                       output_sorted['duration2'][1]]
    #             new = np.array([opt_fun(t_vals[1]), ], dtype=[('duration2', np.float64), ('max', np.float64),
    #                                                           ('mean', np.float64), ('min', np.float64)])
    #             output = np.hstack((output_sorted[:2], new))
    #             output_sorted = np.sort(output, order='max')  # low to high
    #             ind = np.searchsorted(output_sorted['max'], limit)
    #
    #         elif ind == 2:
    #             t_vals = [output_sorted['duration2'][1], np.mean(output_sorted['duration2'][1:]),
    #                       output_sorted['duration2'][2]]
    #             new = np.array([opt_fun(t_vals[1]), ], dtype=[('duration2', np.float64), ('max', np.float64),
    #                                                           ('mean', np.float64), ('min', np.float64)])
    #             output = np.hstack((output_sorted[1:], new))
    #             output_sorted = np.sort(output, order='max')  # low to high
    #             ind = np.searchsorted(output_sorted['max'], limit)
    #
    #         else:
    #             return None  # this should never happen, maybe throw an error here
    #
    #     f_dwell_2_time = interpolate.interp1d(output_sorted['max'], output_sorted['duration2'], assume_sorted=False)
    #     f_min_temp = interpolate.interp1d(output_sorted['max'], output_sorted['min'], assume_sorted=False)
    #     f_mean_temp = interpolate.interp1d(output_sorted['max'], output_sorted['mean'], assume_sorted=False)
    #
    #     results['max_temp'] = limit
    #     results['dwell_2_time'] = f_dwell_2_time(limit)
    #     results['min_temp'] = f_min_temp(limit)
    #     results['mean_temp'] = f_mean_temp(limit)
    #     results['converged'] = True
    #
    #     # print('Number of Iterations: {}'.format(n))
    #
    #     return None

    datesecs = DateTime(date).secs

    if max_dwell is None:
        # This ensures three "cycles" of the two dwell states, within the portion of the schedule used for evaluation
        # (t_backoff).
        # Subtract 1000 sec for extra padding.
        max_dwell = (t_backoff - t_dwell1) / 3 - 1000

    results = {'converged': False, 'unconverged_hot': False, 'unconverged_cold': False,
               'min_temp': np.nan, 'mean_temp': np.nan, 'max_temp': np.nan, 'temperature_limit': limit,
               'dwell_2_time': np.nan, 'dwell_2_time_limit': max_dwell, 'min_pseudo': np.nan, 'mean_pseudo': np.nan,
               'max_pseudo': np.nan, 'hotter_state': np.nan, 'colder_state': np.nan}

    # Ensure t_dwell1 is a float, may not be necessary anymore
    t_dwell1 = np.float(t_dwell1)

    opt_fun = create_opt_fun(datesecs, dwell1_state, dwell2_state, t_dwell1, msid, model_spec, init,
                             t_backoff, duration)

    # dwell2_range = np.logspace(1.0e-6, 1, n_dwells, endpoint=True) / n_dwells
    # dwell2_range = max_dwell * (dwell2_range - dwell2_range[0]) / (dwell2_range[-1] - dwell2_range[0])
    # output = np.array([opt_fun(t) for t in dwell2_range],
    #                   dtype=[('duration2', np.float64), ('max', np.float64), ('mean', np.float64), ('min', np.float64)])

    output = np.array([opt_fun(t) for t in [1.0e-6, max_dwell]],
                      dtype=[('duration2', np.float64), ('max', np.float64), ('mean', np.float64), ('min', np.float64)])

    if np.all(output['max'] < limit):
        results['dwell_2_time'] = np.nan
        ind = np.argmax(output['max'])
        results['max_temp'] = output['max'][ind]
        results['min_temp'] = output['min'][ind]
        results['mean_temp'] = output['mean'][ind]
        results['converged'] = False
        results['unconverged_cold'] = True

    elif np.all(output['max'] > limit):
        results['dwell_2_time'] = np.nan
        ind = np.argmin(output['max'])
        results['max_temp'] = output['max'][ind]
        results['min_temp'] = output['min'][ind]
        results['mean_temp'] = output['mean'][ind]
        results['converged'] = False
        results['unconverged_hot'] = True

    else:

        # --------------------------------------------------------------------------------------------------------------
        n_dwells_1 = n_dwells
        dwell2_range = np.logspace(1.0e-6, 1, n_dwells_1, endpoint=True) / n_dwells_1
        dwell2_range = max_dwell * (dwell2_range - dwell2_range[0]) / (dwell2_range[-1] - dwell2_range[0])
        output = np.array([opt_fun(t) for t in dwell2_range],
                          dtype=[('duration2', np.float64), ('max', np.float64), ('mean', np.float64),
                                 ('min', np.float64)])

        output_sorted = np.sort(output, order='max') # low to high
        ind = np.searchsorted(output_sorted['max'], limit)

        if ind == 0:
            # np.searchsorted finds the first suitable location by default, so if ind == 0, then the duration must fall
            # at the bounded value. This is not true if ind == -1 (the last value).
            results['max_temp'] = limit
            results['dwell_2_time'] = output['duration2'][ind]
            results['min_temp'] = output['min'][ind]
            results['mean_temp'] = output['mean'][ind]
            results['converged'] = True

        else:
            n_dwells_2 = n_dwells
            t_bound = (output_sorted['duration2'][ind - 1], output_sorted['duration2'][ind])
            dwell2_range = np.linspace(np.min(t_bound), np.max(t_bound), n_dwells_2, endpoint=True)
            output = np.array([opt_fun(t) for t in dwell2_range],
                              dtype=[('duration2', np.float64), ('max', np.float64), ('mean', np.float64),
                                     ('min', np.float64)])

            f_dwell_2_time = interpolate.interp1d(output['max'], output['duration2'], kind='quadratic', assume_sorted=False)
            f_min_temp = interpolate.interp1d(output['max'], output['min'], kind='quadratic', assume_sorted=False)
            f_mean_temp = interpolate.interp1d(output['max'], output['mean'], kind='quadratic', assume_sorted=False)

            results['max_temp'] = limit
            results['dwell_2_time'] = f_dwell_2_time(limit)
            results['min_temp'] = f_min_temp(limit)
            results['mean_temp'] = f_mean_temp(limit)
            results['converged'] = True

        # --------------------------------------------------------------------------------------------------------------

        # opt_dwell2()

    if output['max'][0] > output['max'][-1]:
        results['hotter_state'] = 1
        results['colder_state'] = 2
    else:
        results['hotter_state'] = 2
        results['colder_state'] = 1

    return results, output





def run_state_pairs(msid, model_spec, init, limit, date, state_pairs, max_dwell=None, n_dwells=10, pseudo=None,
                    shared_data=None):
    """ Determine dwell balance times for a set of cases.

    Args:
        msid: msid (str): Primary MSID for model being run
        model_spec (dict, string): Dictionary of model parameters or file location where parameters can be imported
        init (dict): Dictionary of Xija model initialization parameters, can be empty
        limit (float): Temperature limit for primary MSID in model for this simulation
        date (float, int, str): Date for start of simulation, in seconds from '1997:365:23:58:56.816', or any other
            format readable by Chandra.Time.DateTime
        state_pairs: Iterable of dictionary pairs, where each pair of dictionaries contain dwell1 and dwell2 states, see
            state_pair section below for further details
        max_dwell (float): Maximum duration for second dwell, can be tuned to provide better results
        n_dwells (int): Number of second dwell possibilities to run (more dwells = finer resolution)
        pseudo (:obj:`str`, optional): Name of one or more pseudo MSIDs used in the model, if any, only necessary if one
            wishes to retrieve model results for this pseudo node, if it exists - To be implemented at a later date,
            not currently used
        shared_data (list): Shared list of results, used when running multiple `run_state_pairs` threads in parallel via
            the multiprocessing package

    Returns:
        Structured numpy array of results


    State Pairs Data Structure:
        The first dictionary in the pair has the following structure, the minimum fields required are shown:
        sequence1: Numerical value to identify the first dwell segment, can be user defined, does not need to match
            LTS value
        obsid1: Obsid for first dwell segment
        duration1_fraction: The fraction of total time for this Obsid in the current week represented by 'duration1'
        duration1: Fixed dwell time for first dwell segment
        pitch: Pitch for first dwell segment

        The second dictionary in the pair has the following structure, the minimum fields required are shown:
        sequence2: Numerical value to identify the second dwell segment, can be user defined, does not need to match
            LTS value
        obsid2: Obsid for second dwell segment
        pitch: Pitch for second dwell segment

        State information that does not change from dwell1 to dwell2 can be specified in the model initialization
        dictionary. State information that does change from dwell1 to dwell2 should be specified in the state pairs
        dictionary described above. Dictionary names for states should match those expected by Xija (e.g. fep_count,
        roll, sim_z).


    Example:

        model_init = {'aacccdpt': {'aacccdpt': -10., 'aca0': -10., 'eclipse': False}, }

        model_specs = load_model_specs()
        date = '2019:001:00:00:00'
        t_dwell1 = 20000.
        msid = 'aacccdpt'
        limit = -9.5

        state_pairs = (({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                         'pitch': 144.2}, {'sequence2': 2000, 'obsid2': 22222, 'pitch': 154.95}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 90.2}, {'sequence2': 3000, 'obsid2': 33333,'pitch': 148.95}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 50}, {'sequence2': 4000, 'obsid2': 44444,'pitch': 140}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 90}, {'sequence2': 5000, 'obsid2': 55555,'pitch': 100}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 75}, {'sequence2': 6000, 'obsid2': 66666,'pitch': 130}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 170}, {'sequence2': 7000, 'obsid2': 77777,'pitch': 90}),
                       ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1,
                          'pitch': 90}, {'sequence2': 8000, 'obsid2': 88888,'pitch': 170}))

        results = run_state_pairs(msid, model_specs[msid], model_init[msid], limit, date, state_pairs)

    """

    duration = 30 * 24 * 3600.
    t_backoff = 2 * duration / 3
    datestr = DateTime(date).date[:8]
    datesecs = DateTime(date).secs

    results = []

    num = np.float(len(state_pairs))
    for n, pair in enumerate(state_pairs):

        if np.mod(n, 1000) == 0:
            print("Running simulations for state pair #: {} out of {}".format(n + 1, num))

        dwell1_state = pair[0]
        dwell2_state = pair[1]

        # These need to be pulled out of the dwell state data structures, otherwise they will cause errors later when
        # constructing the state "schedule".
        duration1 = dwell1_state.pop('duration1')
        duration1_fraction = dwell1_state.pop('duration1_fraction')
        sequence1 = dwell1_state.pop('sequence1')
        obsid1 = dwell1_state.pop('obsid1')

        sequence2 = dwell2_state.pop('sequence2')
        obsid2 = dwell2_state.pop('obsid2')

        dwell_results, output = find_second_dwell(date, dwell1_state, dwell2_state, duration1, msid, limit, model_spec,
                                                  init, duration=duration, t_backoff=t_backoff, n_dwells=n_dwells,
                                                  max_dwell=max_dwell, pseudo=None)

        row = (msid,
               datestr,
               datesecs,
               obsid1,
               sequence1,
               duration1_fraction,
               obsid2,
               sequence2,
               limit,
               dwell1_state['pitch'],
               dwell2_state['pitch'],
               dwell1_state['roll'] if 'roll' in dwell1_state else 0,
               dwell2_state['roll'] if 'roll' in dwell1_state else 0,
               dwell1_state['ccd_count'] if 'ccd_count' in dwell1_state else 0,
               dwell2_state['ccd_count'] if 'ccd_count' in dwell1_state else 0,
               dwell1_state['fep_count'] if 'fep_count' in dwell1_state else 0,
               dwell2_state['fep_count'] if 'fep_count' in dwell1_state else 0,
               dwell1_state['clocking'] if 'clocking' in dwell1_state else 0,
               dwell2_state['clocking'] if 'clocking' in dwell1_state else 0,
               dwell1_state['vid_board'] if 'vid_board' in dwell1_state else 0,
               dwell2_state['vid_board'] if 'vid_board' in dwell1_state else 0,
               dwell1_state['sim_z'] if 'sim_z' in dwell1_state else 0,
               dwell2_state['sim_z'] if 'sim_z' in dwell1_state else 0,
               duration1,
               dwell_results['dwell_2_time'],
               dwell_results['min_temp'],
               dwell_results['mean_temp'],
               dwell_results['max_temp'],
               dwell_results['min_pseudo'],
               dwell_results['mean_pseudo'],
               dwell_results['max_pseudo'],
               dwell_results['converged'],
               dwell_results['unconverged_hot'],
               dwell_results['unconverged_cold'],
               dwell_results['hotter_state'],
               dwell_results['colder_state'])

        results.append(row)

    results = np.array(results, dtype=results_dtype)

    if shared_data is not None:
        shared_data.append(results)
    else:
        return results


if __name__ == '__main__':

    t1 = DateTime().secs

    model_init = {'aacccdpt': {'aacccdpt': -10., 'aca0': -10., 'eclipse': False},
                  'pftank2t': {'pftank2t': f_to_c(95.), 'pf0tank2t': f_to_c(95.), 'eclipse': False},
                  'tcylaft6': {'tcylaft6': f_to_c(120.), 'cc0': f_to_c(120.), 'eclipse': False},
                  '4rt700t': {'4rt700t': f_to_c(95.), 'oba0': f_to_c(95.), 'eclipse': False},
                  '1dpamzt': {'1dpamzt': 35., 'dpa0': 35., 'eclipse': False, 'vid_board': True, 'clocking': True,
                              'dpa_power': 0.0, 'sim_z': 100000},
                  '1deamzt': {'1deamzt': 35., 'eclipse': False, 'vid_board': True, 'clocking': True, 'dpa_power': 0.0,
                              'sim_z': 100000}}

    model_specs = load_model_specs()

    date = '2019:001:00:00:00'
    scale_factor = 1000.
    t_dwell1 = 20000.

    msid = 'aacccdpt'
    limit = -9.5

    state_pairs = (({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 144.2}, {'sequence2': 2000, 'obsid2': 22222, 'pitch': 154.95}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 90.2}, {'sequence2': 3000, 'obsid2': 33333,'pitch': 148.95}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 50}, {'sequence2': 4000, 'obsid2': 44444,'pitch': 140}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 90}, {'sequence2': 5000, 'obsid2': 55555,'pitch': 100}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 75}, {'sequence2': 6000, 'obsid2': 66666,'pitch': 130}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 170}, {'sequence2': 7000, 'obsid2': 77777,'pitch': 90}),
                   ({'sequence1': 1000, 'obsid1': 99999, 'duration1_fraction': 1.0, 'duration1': t_dwell1, 'pitch': 90}, {'sequence2': 8000, 'obsid2': 88888,'pitch': 170}))
    results = run_state_pairs(msid, model_specs[msid], model_init[msid], limit, date, state_pairs, n_dwells=20)

    print(results)
    print('MD5 sum for ACA model: {}'.format(model_specs['aacccdpt_hash']))

    t2 = DateTime().secs

    print('\nRunning {} state pairs tooks {} seconds'.format(len(state_pairs), t2 - t1))