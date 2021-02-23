import numpy as np
from timbre import run_state_pairs, f_to_c
from pathlib import Path
import astropy
import matplotlib.pyplot as plt
import json
from os import listdir

# shortpath to my model spec folder
XIJA_PATH = "C:/Users/JKRIST~1/DOCUME~1/MATLAB/FOT_TO~2/THERMA~1/CHANDR~1.33/CHANDR~1/xija"
# map of model parameter to dtype
PARAM_DTYPES = {'pitch': np.float64, 'roll': np.float64, 'fep_count': np.int8,
                'ccd_count': np.int8, 'clocking': np.bool, 'vid_board': np.bool, 'sim_z': np.int32}
# Default Properties for any models in this analysis
DEFAULT_PARAMS = {'roll': np.float64(0.0), 'fep_count': np.int8(3), 'ccd_count': np.int16(3),
                  'clocking': np.bool(True), 'vid_board': np.bool(True), 'sim_z': np.int32(100000)}
# data type specifications for reading in *.csv files
CONVERTERS = {k: [v] for k, v in zip(
    ['time', 'model', 'pitch', 'duration', 'chips', 'sim_z', 'group id', 'side id', 'is_cold'],
    [astropy.io.ascii.convert_numpy(np.float), astropy.io.ascii.convert_numpy(np.str),
     astropy.io.ascii.convert_numpy(np.float), astropy.io.ascii.convert_numpy(np.float),
     astropy.io.ascii.convert_numpy(np.uint8), astropy.io.ascii.convert_numpy(np.int32),
     astropy.io.ascii.convert_numpy(np.str), astropy.io.ascii.convert_numpy(np.bool),
     astropy.io.ascii.convert_numpy(np.bool)])}
# initialize model_params with the defaults
model_specific_params = DEFAULT_PARAMS.copy()
# initialize balance group id to 0
balance_group_id = 0


# define the data unique to each model
MODEL_INFO = {
    'aacccdpt': {
        'limit': -6.5, 'spec_file': Path(XIJA_PATH + '/aca/aca_spec.json'),
        'model_init': {'aacccdpt': -6.5, 'aca0': -6.5, 'eclipse': False},
        'params': []},
    'pftank2t': {
        'limit': f_to_c(110), 'spec_file': Path(XIJA_PATH + '/pftank2t/pftank2t_spec.json'),
        'model_init': {'pftank2t': f_to_c(110), 'pf0tank2t': f_to_c(110), 'eclipse': False},
        'params': ['roll']},
    '1deamzt': {
        'limit': 37.5, 'spec_file': Path(XIJA_PATH + '/dea/dea_spec.json'),
        'model_init': {'1deamzt': 37.5, 'dea0': 37.5, 'eclipse': False, 'dpa_power': 0.0},
        'params': ['roll', 'fep_count', 'ccd_count', 'clocking', 'vid_board', 'sim_z']},
    '1dpamzt': {
        'limit': 37.5, 'spec_file': Path(XIJA_PATH + '/dpa/dpa_spec.json'),
        'model_init': {'1dpamzt': 37.5, 'dpa0': 37.5, 'eclipse': False, 'dpa_power': 0.0},
        'params': ['roll', 'fep_count', 'ccd_count', 'clocking', 'vid_board', 'sim_z']}}


def build_state_pairs(starting_pitch, model_params):
    """
    builds a list of state pairs.  These pairs vary in the destination pitch.
    """
    return [({'pitch': starting_pitch, **model_params}, {'pitch': p, **model_params}) for p in range(40, 181, 5)]


def find_balance_group(pitch, dwell_time, date, models=None):
    """
    determines the dwells that balance the given dwell/pitch combination as
    well as other dwells/pitches that are equivalent to the given inputs
    """
    if models is None:
        model_info = MODEL_INFO
    else:
        model_info = {k: v for k, v in zip(MODEL_INFO.keys(), MODEL_INFO.values()) if k in models}
    results = []
    for model in model_info:
        # Here is the list of all the parameters needed for this model
        model_param_keys = np.append(model_info[model]['params'], 'pitch')
        # build up the state pairs based on the constant parameters and the given pitch
        model_params = {k: v for k, v in zip(model_specific_params.keys(), model_specific_params.values())
                        if k in model_param_keys}
        state_pairs = build_state_pairs(pitch, model_params)
        # build a dtype array for this model
        state_pair_dtype = {k: v for k, v in zip(PARAM_DTYPES.keys(), PARAM_DTYPES.values()) if k in model_param_keys}
        # determine the limit type
        limit_type = model_info[model].get('limit_type', 'max')
        # run the state pairs
        print(f'starting {model} anaysis\n')
        tmp_results = run_state_pairs(model,
                                      model_info[model]['spec_file'],
                                      model_info[model]['model_init'],
                                      model_info[model]['limit'],
                                      date, dwell_time,
                                      state_pairs, state_pair_dtype,
                                      limit_type=limit_type)
        # capture results
        results.append(tmp_results)
    return results


def get_converged_table(results_data):
    table_data = astropy.table.Table(results_data)
    idx = table_data['converged']
    return table_data[idx]


def set_model_specific_params(or_data, or_idx, chip_cnt):
    """
    Set the model_specific_params for the current option
    """
    model_specific_params['fep_count'] = chip_cnt
    model_specific_params['ccd_count'] = chip_cnt
    model_specific_params['sim_z'] = or_data['or_simpos'][or_idx]


def get_groups(pitch, duration, time_str):
    global balance_group_id
    PLOT = False
    msids = [k for k in MODEL_INFO.keys()]
    if PLOT:
        fig, axs = plt.subplots(len(msids))
    column_names = ['time', 'model', 'pitch', 'duration', 'chips', 'sim_z', 'group id', 'side id', 'is_cold']
    column_types = [np.float, str, np.float, np.float, np.uint8, np.int32, str, np.bool, np.bool]
    results = astropy.table.Table(names=column_names, dtype=column_types)
    # the groups that balances the starting pitch/dwell for all models
    group_b = find_balance_group(pitch, duration, time_str)
    for r_idx in range(0, len(msids)):
        msid = msids[r_idx]
        # the converged group that balances the starting pitch/dwell for the current model
        g2 = get_converged_table(group_b[r_idx])
        tmp = g2[['datesecs', 'msid', 'pitch2', 't_dwell2']]
        tmp.rename_columns(['datesecs', 'msid', 'pitch2', 't_dwell2'], ['time', 'model', 'pitch', 'duration'])
        tmp.add_column(g2['colder_state'] == 2, name='is_cold')
        tmp.add_column(model_specific_params['fep_count'], name='chips')
        tmp.add_column(model_specific_params['sim_z'], name='sim_z')
        tmp.add_column(f'{balance_group_id:#0{6}x}', name='group id')
        tmp.add_column(True, name='side id')
        results = astropy.table.vstack([results, tmp])
        # add in the state data for the actual OR
        tmp_row = list(tmp[['time', 'model', 'chips', 'sim_z', 'group id', 'side id', 'is_cold']][0])
        tmp_row.insert(2, pitch)
        tmp_row.insert(3, duration)
        # flip these because we are starting with data for the group that balances the OR
        tmp_row[-2:] = np.logical_not(tmp_row[-2:])
        results.add_row(tmp_row)
        # the group that includes the starting pitch/dwell
        group_a = find_balance_group(g2['pitch2'][0], g2['t_dwell2'][0], time_str, models=msid)
        g1 = get_converged_table(group_a[0])
        tmp = g1[['datesecs', 'msid', 'pitch2', 't_dwell2']]
        tmp.rename_columns(['datesecs', 'msid', 'pitch2', 't_dwell2'], ['time', 'model', 'pitch', 'duration'])
        tmp.add_column(g1['colder_state'] == 2, name='is_cold')
        tmp.add_column(model_specific_params['fep_count'], name='chips')
        tmp.add_column(model_specific_params['sim_z'], name='sim_z')
        tmp.add_column(f'{balance_group_id:#0{6}x}', name='group id')
        balance_group_id += 1
        tmp.add_column(False, name='side id')
        results = astropy.table.vstack([results, tmp])
        results.sort(['time', 'model', 'pitch', 'group id', 'side id'])
        if PLOT:
            axs[r_idx].plot(g2['pitch2'], g2['t_dwell2'], 'g.', g1['pitch2'], g1['t_dwell2'], 'y.')
            axs[r_idx].set_title(msid)
            axs[r_idx].set_xlabel('pitch')
            axs[r_idx].set_ylabel('duration (s)')
            axs[r_idx].set_xlim([40, 180])
            axs[r_idx].set_ylim([-10000, 100000])
    if PLOT:
        plt.show()
    return results


input_path = Path('C:/Users/jkristoff/Documents/MATLAB/FOT_Tools/timbre_inputs.json')
with open(input_path, 'r') as f:
    or_data = json.loads(f.read())

num_ors = len(or_data['or_simpos'])
num_times = len(or_data['sched_t'])
results = astropy.table.Table()
# loop through all the observations
for or_idx in range(0, num_ors):
    pitch = or_data['pitch_deg'][or_idx]
    duration = or_data['dwell_t'][or_idx]
    # determine if there are any optional chips
    num_opts = or_data['opt_chips'][or_idx] + 1
    # for each chip configuration
    for opt_idx in range(0, num_opts):
        chip_cnt = or_data['ccd_cnt'][or_idx] - opt_idx
        # loop through all the available times
        for time_idx in range(0, num_times):
            time_str = or_data['sched_t'][time_idx]
            # setup the model specific parameters and get the balance groups
            set_model_specific_params(or_data, or_idx, chip_cnt)
            tmp = get_groups(pitch[time_idx], duration, time_str)
            # write out the data in case of failure
            tmp.write(f'results_{tmp["group id"][0]}.csv')
            # combine the results
            if len(results) == 0:
                results = tmp
            else:
                results = astropy.table.vstack(results, tmp)
results.write('final_results.csv')


def read_csv_data(file_name):
    return astropy.table.Table.read(Path(file_name), format='csv', converters=CONVERTERS)


def merge_csv_data(directory_name):
    all_files = np.array(listdir(directory_name))
    csv_idx = np.array([file.find('.csv') > 0 for file in all_files], dtype=np.bool)
    csv_files = all_files[csv_idx]
    results = []
    for file_name in csv_files:
        tmp_table = read_csv_data(file_name)
        if len(results) == 0:
            results = tmp_table
        else:
            results = astropy.table.vstack([results, tmp_table])
    return results


q = merge_csv_data('.')

# if the pitches in a given group (g1 or g2) are included in all models
# then the pitches are homogenous (either all cooling or all heating).
# otherwise, they are heterogenous, and heat some parts while cooling others.

# for Homogenous pitches, an energy direction can be assigned (e.g. Heating or Cooling)
# and the equivalent duration for the pitch can be calculated as:
# eq_dur = min(all_durations_at_pitch) if energy_direction is Heating
# eq_dur = max(all_durations_at_pitch) if energy_direction is Cooling
# this biases the equalibrium towards cooling more (although this might need to be
# reconsidered for models that have a 'min' limit_type)

# Maybe each pitch should have a duration mean and standard deviation, instead of a
# specific value, so I can use it as a continuous variable in an optimization, e.g.
# given obs_0 has been scheduled, and a list of pitches and durations, which observation
# fits most closely to the mean.  Or, should I use the pre-computed balance data to see
# which observation has obs_0 closest to it's balancing pitch/duration curves.

# QUESTION1: Minimum duration, maximum duration, mean duration, or none of the above?
# QUESTION2: Can multiple models be combined?  Is yes, how to do it?

# treat each model independently and prioritize by distance from reference?
# mission planners want to "stay as close to the limits as possible" but I think
# this actually means that they want to use as much time on an observation as possible.
#
# Concept: a model far from its limit has a larger "buffer" and can tolerate a larger error
#          a model close to its limit has a smaller "buffer" and cannot tolerate much error
#
# at each scheduling step, we may need to:
# 1. calculate the difference between each model's current value and its limit
# 2. if any difference indicates a limit has been crossed, stop and mark the last observation as
#  invalid (will want to calculate the amount of the observation that could be scheduled also)
# 3. calculate a tolerance based on the difference to the reference for each model.
# 4. use the calculated tolerance in a loose matching to determine if the next observation
#  matches the balancing data for the last observation.
# 4.a. calculate the mse of the next observation and store it for prioritization.
# 5. do step 4 for all the remaining observations.
# 6. organize the list of observations by:
#   1. the observations that loosly matched the balance group data (prioritied by lowest mse on top)
#   2.
