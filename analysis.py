import numpy as np
from timbre import run_state_pairs, f_to_c
from pathlib import Path
import astropy

# shortpath to my model spec folder
XIJA_PATH = "C:/Users/JKRIST~1/DOCUME~1/MATLAB/FOT_TO~2/THERMA~1/CHANDR~1.33/CHANDR~1/xija"
# map of model parameter to dtype
PARAM_DTYPES = {'pitch': np.float64, 'roll': np.float64, 'fep_count': np.int8,
                'ccd_count': np.int8, 'clocking': np.bool, 'vid_board': np.bool, 'sim_z': np.int32}
# Properties no changing for any models in this analysis
CONST_PARAMS = {'roll': np.float64(0.0), 'fep_count': np.int8(4), 'ccd_count': np.int8(4),
                'clocking': np.bool(True), 'vid_board': np.bool(True), 'sim_z': np.int32(100000)}


def build_state_pairs(starting_pitch, model_params):
    """
    builds a list of state pairs.  These pairs vary in the destination pitch.
    """
    return [({'pitch': starting_pitch, **model_params}, {'pitch': p, **model_params}) for p in range(45, 180, 5)]


def find_balance_group(pitch, dwell_time, date):
    """
    determines the dwells that balance the given dwell/pitch combination as
    well as other dwells/pitches that are equivalent to the given inputs
    """
    # define the data unique to each model
    model_info = {
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

    results = []
    for model in model_info:
        # Here is the list of all the parameters needed for this model
        model_param_keys = np.append(model_info[model]['params'], 'pitch')
        # build up the state pairs based on the constant parameters and the given pitch
        model_params = {k: v for k, v in zip(CONST_PARAMS.keys(), CONST_PARAMS.values()) if k in model_param_keys}
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


results1 = find_balance_group(127, 10000, '2021:050:12:30:00')
r1 = astropy.table.Table(results1[0])
idx = r1['converged']
r1 = r1[idx]

results2 = find_balance_group(r1[0]['pitch2'], r1[0]['t_dwell2'], '2021:050:12:30:00')
r2 = astropy.table.Table(results2[0])
idx = r2['converged']
r2 = r2[idx]
