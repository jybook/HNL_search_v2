import numpy as np
import h5py
import argparse

from icecube import icetray, dataclasses, dataio, photonics_service, phys_services, recclasses
from I3Tray import I3Units, I3Tray

from collections import OrderedDict
import itertools
import random

from icecube.hdfwriter import I3HDFWriter

import time
import glob

## Create ability to change settings from terminal ##
parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", 
    "--input",
    type=str, 
    dest="input_file", 
    help="path and name of the input file"
)
parser.add_argument(
    "-n", 
    "--name", 
    type=str, 
    dest="output_name", 
    help="name for output file (no path)"
)
parser.add_argument(
    "-o", 
    "--outdir", 
    type=str, 
    dest="output_dir",
    help="path of ouput file"
)
parser.add_argument(
    "--model_dir",
    type=str,
    dest="model_dir",
    help="path for where to pull the model from",
)
parser.add_argument(
    "--modelname_list",
    nargs="+",
    default=[],
    dest="modelname_list",
    help="name of the models to use",
)
parser.add_argument(
    "--variable_list",
    nargs="+",
    default=[],
    dest="variable_list",
    help="names of variables to predict",
)
parser.add_argument(
    "-e",
    "--epochs_list",
    nargs="+",
    default=[],
    dest="epochs_list",
    help="epoch numbers for models to use",
)
parser.add_argument(
    "-f",
    "--factor_list",
    nargs="+",
    default=[100.0, 1, 1, 1, 1, 1, 1, 1],
    dest="factor_list",
    help="transformation factors to adjust output by",
)
parser.add_argument(
    "--cleaned",
    type=str,
    default="True",
    dest="cleaned",
    help="True if wanted to use SRTTWOfflinePulsesDC",
)
parser.add_argument(
    "--charge_min",
    type=float,
    default=0.25,
    dest="charge_min",
    help="minimum charge pulse to keep, remove < this",
)
parser.add_argument(
    "--hdf5",
    action="store_false",
    dest="write_hdf5",
    default=True,
    help="Also store hdf5 file"
)
parser.add_argument(
    "--real_sim",
    type=int,
    dest="real_sim",
    default=1,
    help="Is this the real (model dependent) simulation?"
)
parser.add_argument(
    "--gcd",
    dest="gcd",
    default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
    help="path and filename of gcd",
)

args = parser.parse_args()

input_file = args.input_file
output_dir = args.output_dir
output_name = args.output_name
gcdfile = args.gcd
write_hdf5  = args.write_hdf5
real_sim    = bool(args.real_sim)

if args.cleaned == "True" or args.cleaned == "true":
    use_cleaned_pulses = True
else:
    use_cleaned_pulses = False
charge_min = args.charge_min

variable_list = args.variable_list
scale_factor = args.factor_list
scale_factor_list = np.array(scale_factor, dtype=float)
modelname_list = args.modelname_list
epoch_list = args.epochs_list
epoch_list = np.array(epoch_list, dtype=float)
print(len(epoch_list))
model_path = args.model_dir

accepted_names = ["energy", "zenith", "class", "vertex", "muon", "muonL4", "nDOM", "ending", "muonV3"]
for var in variable_list:
    assert (
        var in accepted_names
    ), "Variable must be one of the accepted names, check parse arg help for variable for more info"

model_name_list = []
num_variables = len(variable_list)
for variable_index in range(num_variables):
    if variable_list[variable_index] == "nDOM":
        model_name = "Not_CNN"
    else:
        if len(epoch_list) == 0:
            model_name = args.model_dir + "/" + modelname_list[variable_index] + ".hdf5"
        else:
            if epoch_list[variable_index] is None:
                model_name = (
                    args.model_dir + "/" + modelname_list[variable_index] + ".hdf5"
                )
            else:
                model_name = "%s/%s/%s_%sepochs_model.hdf5" % (
                    args.model_dir,
                    modelname_list[variable_index],
                    modelname_list[variable_index],
                    epoch_list[variable_index],
                )
    model_name_list.append(model_name)
    print(
        "Predicting: %s,\nOutput transformation scale factor: %.2f.,\nUsing model: %s"
        % (variable_list[variable_index], scale_factor_list[variable_index], model_name)
    )


########## Functions needed to store as hdf5 ##########

HDF5_KEYS = [
    # EventProperties
    "decay_channel",
    "distance",
    "distanceMax",
    "distanceMin",
    "finalStateX",
    "finalStateY",
    "final_state_particle0",
    "final_state_particle1",
    "primary_type",
    "lifetime",
    "mHNL",
    "outgoing_neutrino_energy",
    "totalEnergy",
    "physical",
    "total_column_depth",
    # I3MCTree
    # true HNL variables
    "HNL_true_x",
    "HNL_true_y",
    "HNL_true_z",
    "HNL_true_energy",
    "HNL_true_zenith",
    "HNL_true_azimuth",
    "HNL_true_time",
    # true primary variables
    "true_x",
    "true_y",
    "true_z",
    "true_energy",
    "true_zenith",
    "true_azimuth",
    "true_time",
    # true first (DIS) cascade variables
    "casc0_true_x",
    "casc0_true_y",
    "casc0_true_z",
    "casc0_true_energy",
    "casc0_true_zenith",
    "casc0_true_azimuth",
    "casc0_true_time",
    # true second (HNL decay) cascade variables
    "casc1_true_x",
    "casc1_true_y",
    "casc1_true_z",
    "casc1_true_energy",
    "casc1_true_zenith",
    "casc1_true_azimuth",
    "casc1_true_time",
    "nan_decay_energy",
    # weights
    "LeptonInjectorWeight",
    "LifetimeWeight_1e-03",
    "OneWeight",
    "ReferenceWeight_1e-03",
    # total true quantities
    'true_decay_length',
    'true_total_energy',
    # flercnn quantities
    'FLERCNN_BDT_ProbNu',
    'FLERCNN_coszen',
    'FLERCNN_energy',
    'FLERCNN_nDOM',
    'FLERCNN_prob_muon_v3',
    'FLERCNN_prob_track',
    'FLERCNN_vertex_rho36',
    'FLERCNN_vertex_x',
    'FLERCNN_vertex_y',
    'FLERCNN_vertex_z',
    'FLERCNN_ending_x',
    'FLERCNN_ending_y',
    'FLERCNN_ending_z',
    'FLERCNN_zenith',
    # additional length seed options
    'L5_DirectHitsC_DirTrackLength',
    'FiniteRecoFit_length',
    # PID usable variables
    'L7_PIDClassifier_FullSky_ProbTrack',
    'L7_PIDClassifier_Upgoing_ProbTrack',
    'retro_crs_prefit__zero_dllh_mean',
    'retro_crs_prefit__max_llh',
]

# remove duplicates from the HDF5 keys
HDF5_KEYS = list(
    set(HDF5_KEYS)  # note that this does NOT preserve order, but doesn't matter
)

# true variables - model independent simulation case
def write_true_data_to_keys(frame):
    """
    Function to write true cascade information, energy and decay length to frame
    """
    cascades = frame["I3MCTree"].get_daughters(frame["I3MCTree"].primaries[0])

    # add true cascades, energy and decay length to the frame
    frame["true_decay_length"] = dataclasses.I3Double(
        phys_services.I3Calculator.distance(cascades[0], cascades[1])
    )
    frame["true_total_energy"] = dataclasses.I3Double(
        cascades[0].energy + cascades[1].energy
    )
    for count, cascade in enumerate(cascades):
        frame["casc{}_true_x".format(count)] = dataclasses.I3Double(cascade.pos.x)
        frame["casc{}_true_y".format(count)] = dataclasses.I3Double(cascade.pos.y)
        frame["casc{}_true_z".format(count)] = dataclasses.I3Double(cascade.pos.z)
        frame["casc{}_true_energy".format(count)] = dataclasses.I3Double(cascade.energy)
        frame["casc{}_true_zenith".format(count)] = dataclasses.I3Double(cascade.dir.zenith)
        frame["casc{}_true_azimuth".format(count)] = dataclasses.I3Double(cascade.dir.azimuth)
        frame["casc{}_true_time".format(count)] = dataclasses.I3Double(cascade.time)

    return True

# true variables - model specific HNL simulation case
def store_mc_true_variables(frame):
    if not frame.Has("I3MCTree"):
        return False
    mctree = frame["I3MCTree"]
    p_true = mctree.primaries[0]

    p_daughters = mctree.get_daughters(p_true)
    assert len(p_daughters) == 2

    for p_daughter in p_daughters:
        if p_daughter.type == dataclasses.I3Particle.Hadrons:
            casc_0_true = p_daughter
        else:
            hnl_true = p_daughter

    hnl_daughters = mctree.get_daughters(hnl_true)
    assert len(hnl_daughters) > 0

    nan_energy = False

    for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):
        if not count_hnl_daughters:
            casc_1_true = hnl_daughter
        else:
            assert casc_1_true.pos == hnl_daughter.pos
            if np.isnan(hnl_daughter.energy):
                nan_energy = True
            casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy

    frame["true_x"] = dataclasses.I3Double(p_true.pos.x)
    frame["true_y"] = dataclasses.I3Double(p_true.pos.y)
    frame["true_z"] = dataclasses.I3Double(p_true.pos.z)
    frame["true_energy"] = dataclasses.I3Double(p_true.energy)
    frame["true_zenith"] = dataclasses.I3Double(p_true.dir.zenith)
    frame["true_azimuth"] = dataclasses.I3Double(p_true.dir.azimuth)
    frame["true_time"] = dataclasses.I3Double(p_true.time)

    frame["HNL_true_x"] = dataclasses.I3Double(hnl_true.pos.x)
    frame["HNL_true_y"] = dataclasses.I3Double(hnl_true.pos.y)
    frame["HNL_true_z"] = dataclasses.I3Double(hnl_true.pos.z)
    frame["HNL_true_energy"] = dataclasses.I3Double(hnl_true.energy)
    frame["HNL_true_zenith"] = dataclasses.I3Double(hnl_true.dir.zenith)
    frame["HNL_true_azimuth"] = dataclasses.I3Double(hnl_true.dir.azimuth)
    frame["HNL_true_time"] = dataclasses.I3Double(hnl_true.time)

    frame["casc0_true_x"] = dataclasses.I3Double(casc_0_true.pos.x)
    frame["casc0_true_y"] = dataclasses.I3Double(casc_0_true.pos.y)
    frame["casc0_true_z"] = dataclasses.I3Double(casc_0_true.pos.z)
    frame["casc0_true_energy"] = dataclasses.I3Double(casc_0_true.energy)
    frame["casc0_true_zenith"] = dataclasses.I3Double(casc_0_true.dir.zenith)
    frame["casc0_true_azimuth"] = dataclasses.I3Double(casc_0_true.dir.azimuth)
    frame["casc0_true_time"] = dataclasses.I3Double(casc_0_true.time)

    frame["casc1_true_x"] = dataclasses.I3Double(casc_1_true.pos.x)
    frame["casc1_true_y"] = dataclasses.I3Double(casc_1_true.pos.y)
    frame["casc1_true_z"] = dataclasses.I3Double(casc_1_true.pos.z)
    frame["casc1_true_energy"] = dataclasses.I3Double(casc_1_true.energy)
    frame["casc1_true_zenith"] = dataclasses.I3Double(casc_1_true.dir.zenith)
    frame["casc1_true_azimuth"] = dataclasses.I3Double(casc_1_true.dir.azimuth)
    frame["casc1_true_time"] = dataclasses.I3Double(casc_1_true.time)

    frame["nan_decay_energy"] = icetray.I3Bool(nan_energy)

    return True

# function to get EventProperties (LI) - this won't run with the standard icetray software
def store_LI_event_properties(frame):
    if not frame.Has("EventProperties"):
        return False
    event_properties = frame["EventProperties"]

    frame["mHNL"] = dataclasses.I3Double(event_properties.mHNL)
    frame["distance"] = dataclasses.I3Double(event_properties.distance)
    frame["distanceMin"] = dataclasses.I3Double(event_properties.distanceMin)
    frame["distanceMax"] = dataclasses.I3Double(event_properties.distanceMax)
    frame["lifetime"] = dataclasses.I3Double(event_properties.lifetime)
    frame["totalEnergy"] = dataclasses.I3Double(event_properties.totalEnergy)
    frame["finalStateX"] = dataclasses.I3Double(event_properties.finalStateX)
    frame["finalStateY"] = dataclasses.I3Double(event_properties.finalStateY)
    frame['primary_type'] =  dataclasses.I3Double(event_properties.initialType)
    frame['final_state_particle0'] =  dataclasses.I3Double(event_properties.finalType1)
    frame['final_state_particle1'] =  dataclasses.I3Double(event_properties.finalType2)
    frame['total_column_depth'] =  dataclasses.I3Double(event_properties.totalColumnDepth)
    frame["decay_channel"] = dataclasses.I3Double(event_properties.decay_channel)
    frame["outgoing_neutrino_energy"] = dataclasses.I3Double(event_properties.outgoingNeutrinoEnergy)
    frame["physical"] = dataclasses.I3Double(event_properties.physical)

    return True

# function to get weights
def store_weights(frame):
    if not frame.Has("I3MCWeightDict"):
        return False
    weight_dict = frame["I3MCWeightDict"]

    frame["LeptonInjectorWeight"] = dataclasses.I3Double(weight_dict['LeptonInjectorWeight'])
    frame["LifetimeWeight_1e-03"] = dataclasses.I3Double(weight_dict['LifetimeWeight_1e-03'])
    frame["OneWeight"] = dataclasses.I3Double(weight_dict['OneWeight'])
    frame["ReferenceWeight_1e-03"] = dataclasses.I3Double(weight_dict['ReferenceWeight_1e-03'])

    return True

########## End


def get_observable_features(
    frame,
    low_window=-500,
    high_window=4000,
    use_cleaned_pulses=True,
    charge_min=charge_min,
):
    """
    Load observable features from IceCube files
    Receives:
        frame = IceCube object type from files
    Returns:
        observable_features: Observables dictionary
    """
    if use_cleaned_pulses:
        ice_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
            frame, "SRTTWOfflinePulsesDC"
        )
    else:
        ice_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
            frame, "SplitInIcePulses"
        )

    # First cut: Check if there are 8 cleaned pulses > 0.25 PE
    cleaned_ice_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
        frame, "SRTTWOfflinePulsesDC"
    )
    count_cleaned_pulses = 0
    clean_pulses_8_or_more = False
    for omkey, pulselist in cleaned_ice_pulses:
        if clean_pulses_8_or_more == True:
            break
        for pulse in pulselist:
            a_charge = pulse.charge

            # Cut any pulses < 0.25 PE
            if a_charge < charge_min:
                continue

            # Count number pulses > 0.25 PE in event
            count_cleaned_pulses += 1
            if count_cleaned_pulses >= 8:
                clean_pulses_8_or_more = True
                break
    # Look inside ice pulses and get stats on charges and time
    # DC = deep core which is certain strings/DOMs in IceCube
    store_string = []
    IC_near_DC_strings = [
        17,
        18,
        19,
        25,
        26,
        27,
        28,
        34,
        35,
        36,
        37,
        38,
        44,
        45,
        46,
        47,
        54,
        55,
        56,
    ]
    # IC_near_DC_strings = [26, 27, 35, 36, 37, 45, 46]
    DC_strings = [79, 80, 81, 82, 83, 84, 85, 86]

    # Five summary variables: sum charges, time first pulse, Time of last pulse, Charge weighted mean time of pulses, Charge weighted standard deviation of pulse times
    array_DC = np.zeros(
        [len(DC_strings), 60, 5]
    )  # [string, dom_index, charge & time summary]
    array_IC_near_DC = np.zeros(
        [len(IC_near_DC_strings), 60, 5]
    )  # [string, dom_index, charge & time summary]
    num_pulses_per_dom = np.zeros([len(DC_strings), 60, 1])
    count_outside = 0
    charge_outside = 0
    count_inside = 0
    charge_inside = 0

    # Config 1011 is SMT3
    # dataclasses.TriggerKey(source, ttype, config_id)
    triggers = frame["I3TriggerHierarchy"]
    trigger_time = None
    num_extra_DC_triggers = 0
    for trig in triggers:
        key_str = str(trig.key)
        s = key_str.strip("[").strip("]").split(":")
        if len(s) > 2:
            config_id = int(s[2])
            if config_id == 1011:
                if trigger_time:
                    num_extra_DC_triggers += 1
                trigger_time = trig.time

    if trigger_time == None:
        shift_time_by = 0
    else:
        shift_time_by = trigger_time

    # Start by making all times negative shift time (to distinguish null from 0)
    array_DC[..., 1:] = -20000
    array_IC_near_DC[..., 1:] = -20000

    # Only go through pulse series if we're keeping it
    # if clean_pulses_8_or_more == True:
    for omkey, pulselist in ice_pulses:
        dom_index = omkey.om - 1
        string_val = omkey.string
        timelist = []
        chargelist = []

        DC_flag = False
        IC_near_DC_flag = False

        for pulse in pulselist:
            charge = pulse.charge

            # Cut any pulses < 0.25 PE
            if charge < charge_min:
                continue

            # Quantize pulse chargest to make all seasons appear the same
            quanta = 0.05
            charge = (np.float64(charge) // quanta) * quanta + quanta / 2.0

            if string_val not in store_string:
                store_string.append(string_val)

            # Check IceCube near DeepCore DOMs
            if (string_val in IC_near_DC_strings) and dom_index < 60:
                string_index = IC_near_DC_strings.index(string_val)
                timelist.append(pulse.time)
                chargelist.append(charge)
                IC_near_DC_flag = True

            # Check DeepCore DOMS
            elif (string_val in DC_strings) and dom_index < 60:  # dom_index >=10
                string_index = DC_strings.index(string_val)
                timelist.append(pulse.time)
                chargelist.append(charge)
                DC_flag = True

            else:
                count_outside += 1
                charge_outside += charge

        if DC_flag == True or IC_near_DC_flag == True:
            charge_array = np.array(chargelist)
            time_array = np.array(timelist)
            time_array = [(t_value - shift_time_by) for t_value in time_array]
            time_shifted = [(t_value - time_array[0]) for t_value in time_array]
            time_shifted = np.array(time_shifted)

            # Remove pulses so only those in certain time window are saved
            original_num_pulses = len(timelist)
            time_array_in_window = list(time_array)
            charge_array_in_window = list(charge_array)
            for time_index in range(0, original_num_pulses):
                time_value = time_array[time_index]
                if time_value < low_window or time_value > high_window:
                    time_array_in_window.remove(time_value)
                    charge_array_in_window.remove(charge_array[time_index])
            charge_array = np.array(charge_array_in_window)
            time_array = np.array(time_array_in_window)
            assert len(charge_array) == len(
                time_array
            ), "Mismatched pulse time and charge"
            if len(charge_array) == 0:
                continue

            # Original Stats
            count_inside += len(chargelist)
            charge_inside += sum(chargelist)

            # Check that pulses are sorted in time
            for i_t, time in enumerate(time_array):
                assert time == sorted(time_array)[i_t], "Pulses are not pre-sorted!"

            # Charge weighted mean and stdev
            weighted_avg_time = np.average(time_array, weights=charge_array)
            weighted_std_time = np.sqrt(
                np.average((time_array - weighted_avg_time) ** 2, weights=charge_array)
            )

        if DC_flag == True:
            array_DC[string_index, dom_index, 0] = sum(chargelist)
            array_DC[string_index, dom_index, 1] = time_array[0]
            array_DC[string_index, dom_index, 2] = time_array[-1]
            array_DC[string_index, dom_index, 3] = weighted_avg_time
            array_DC[string_index, dom_index, 4] = weighted_std_time

            num_pulses_per_dom[string_index, dom_index, 0] = len(chargelist)

        if IC_near_DC_flag == True:
            array_IC_near_DC[string_index, dom_index, 0] = sum(chargelist)
            array_IC_near_DC[string_index, dom_index, 1] = time_array[0]
            array_IC_near_DC[string_index, dom_index, 2] = time_array[-1]
            array_IC_near_DC[string_index, dom_index, 3] = weighted_avg_time
            array_IC_near_DC[string_index, dom_index, 4] = weighted_std_time

    return (
        array_DC,
        array_IC_near_DC,
        trigger_time,
        num_extra_DC_triggers,
        clean_pulses_8_or_more,
    )


def apply_transform(
    features_DC,
    features_IC,
    labels=None,
    energy_factor=100.0,
    track_factor=200.0,
    transform="MaxAbs",
):
    from flercenn_auxiliary.scaler_transformations import TransformData
    from flercenn_auxiliary.scaler_transformations import new_transform

    static_stats = [25.0, 4000.0, 4000.0, 4000.0, 2000.0]
    low_stat_DC = static_stats
    high_stat_DC = static_stats
    low_stat_IC = static_stats
    high_stat_IC = static_stats

    features_DC = new_transform(features_DC)
    features_DC = TransformData(
        features_DC, low_stats=low_stat_DC, high_stats=high_stat_DC, scaler=transform
    )
    features_IC = new_transform(features_IC)
    features_IC = TransformData(
        features_IC, low_stats=low_stat_IC, high_stats=high_stat_IC, scaler=transform
    )

    return features_DC, features_IC


def cnn_test(
    features_DC,
    features_IC,
    load_model_name,
    output_variables=1,
    DC_drop_value=0.2,
    IC_drop_value=0.2,
    connected_drop_value=0.2,
    model_type="energy",
):
    if (
        model_type == "class"
        or model_type == "muon"
        or model_type == "muonL4"
        or model_type == "muonV3"
    ):
        from flercenn_auxiliary.cnn_model_classification import make_network
    else:
        from flercenn_auxiliary.cnn_model import make_network

    model_DC = make_network(
        features_DC,
        features_IC,
        output_variables,
        DC_drop_value,
        IC_drop_value,
        connected_drop_value,
    )
    model_DC.load_weights(load_model_name)

    Y_test_predicted = model_DC.predict([features_DC, features_IC])

    return Y_test_predicted


def read_files(filename, gcd_filename=None, use_cleaned_pulses=True, charge_min=0.25):
    """
    Read list of files, make sure they pass L5 cuts, create truth labels
    Receives:
        filename_list = list of strings, filenames to read data from
    Returns:
        output_features_DC = dict with input observable features from the DC strings
        output_features_IC = dict with input observable features from the IC strings
        output_labels = dict with output labels  (energy, zenith, azimith, time, x, y, z,
                        tracklength, isTrack, flavor ID, isAntiNeutrino, isCC)
        output_reco_labels = dict with PegLeg output labels (energy, zenith, azimith, time, x, y, z)
        output_initial_stats = array with info on number of pulses and sum of charge "inside" the strings used
                                vs. "outside", i.e. the strings not used (pulse count outside, charge outside,
                                pulse count inside, charge inside) for finding statistics
        output_num_pulses_per_dom = array that only holds the number of pulses seen per DOM (finding statistics)
        output_trigger_times = list of trigger times for each event (used to shift raw pulse times)
    """

    if gcd_filename is not None:
        gcd_file = dataio.I3File(gcd_filename)
        print("Using GCD file: %s" % gcd_filename)
        pass2_cal = gcd_file.pop_frame(icetray.I3Frame.Calibration)["I3Calibration"]

    print("reading file: {}".format(filename))
    event_file = dataio.I3File(filename)

    output_features_DC = []
    output_features_IC = []
    output_headers = []
    skipped_triggers = 0
    skipped_8hits = 0
    skip_event = []

    for frame in event_file:
        if frame.Stop == icetray.I3Frame.Physics:
            header = frame["I3EventHeader"]

            if header.sub_event_stream != "InIceSplit":
                continue

            if gcd_filename is not None:
                frame["I3Calibration"] = pass2_cal

            (
                DC_array,
                IC_near_DC_array,
                trig_time,
                extra_triggers,
                clean_pulses_8_or_more,
            ) = get_observable_features(
                frame, use_cleaned_pulses=use_cleaned_pulses, charge_min=charge_min
            )

            # Cut events with...
            # Multiple SMT3 tiggers or no SMT3 trigger
            # Less than 8 hits in the cleaned pulse series
            skip = False
            if extra_triggers > 0 or trig_time == None:
                skipped_triggers += 1
                skip = True
            # if clean_pulses_8_or_more == False:
            #    skipped_8hits +=1
            #    skip = True

            skip_event.append(skip)
            header_numbers = np.array(
                [float(header.run_id), float(header.sub_run_id), float(header.event_id)]
            )
            output_headers.append(header_numbers)
            output_features_DC.append(DC_array)
            output_features_IC.append(IC_near_DC_array)

        # close the input file once we are done
    del event_file

    output_headers = np.asarray(output_headers)
    output_features_DC = np.asarray(output_features_DC)
    output_features_IC = np.asarray(output_features_IC)
    skip_event = np.asarray(skip_event)
    print(
        "Number events with 0 or > 1 SMT3 triggers: %i, Number events with less than 8 hits: %i"
        % (skipped_triggers, skipped_8hits)
    )

    return output_features_DC, output_features_IC, output_headers, skip_event


def test_write(
    filename_list,
    model_name_list,
    output_dir,
    output_name,
    model_factor_list=[100.0, 1.0, 1.0, 1.0],
    model_type_list=["energy", "class", "zenith", "vertex", "muon"],
    gcd_file=None,
    use_cleaned_pulses=True,
    charge_min=0.25,
):
    for a_file in filename_list:
        if output_name is None:
            basename = a_file.split("/")[-1]
            basename = basename[:-7]
            output_name = str(basename) + "_FLERCNN"
        outfile_path = output_dir + output_name + ".i3.zst"
        outfile = dataio.I3File(outfile_path, "w")
        print("Writing to %s" % (outfile_path))

        DC_array, IC_near_DC_array, header_array, skip_event = read_files(
            a_file,
            gcd_filename=gcd_file,
            use_cleaned_pulses=use_cleaned_pulses,
            charge_min=charge_min,
        )
        print(DC_array.shape, IC_near_DC_array.shape)

        if DC_array.shape[0] == 0:
            print("THERE ARE P-FRAME EVENTS IN THIS FILE")
        else:
            DC_array, IC_near_DC_array = apply_transform(DC_array, IC_near_DC_array)

            cnn_predictions = []
            for network in range(num_variables):
                if model_type_list[network] == "nDOM":
                    t0 = time.time()
                    charge_DC = DC_array[:, :, :, 0] > 0
                    charge_IC = IC_near_DC_array[:, :, :, 0] > 0
                    DC_flat = np.reshape(charge_DC, [DC_array.shape[0], 480])
                    IC_flat = np.reshape(charge_IC, [IC_near_DC_array.shape[0], 1140])
                    DOMs_hit_DC = np.sum(DC_flat, axis=-1)
                    DOMs_hit_IC = np.sum(IC_flat, axis=-1)
                    DOMs_hit = DOMs_hit_DC + DOMs_hit_IC
                    t1 = time.time()
                    cnn_predictions.append(DOMs_hit)
                    print(
                        "Time to calculate number DOMs hit on %i events: %f seconds"
                        % (DC_array.shape[0], t1 - t0)
                    )
                else:
                    if (
                        model_type_list[network] == "vertex"
                        or model_type_list[network] == "ending"
                    ):
                        output_var = 3
                    else:
                        output_var = 1
                    t0 = time.time()
                    cnn_predictions.append(
                        cnn_test(
                            DC_array,
                            IC_near_DC_array,
                            model_name_list[network],
                            model_type=model_type_list[network],
                            output_variables=output_var,
                        )
                    )
                    t1 = time.time()
                    print(
                        "Time to run CNN Predict %s on %i events: %f seconds"
                        % (model_type_list[network], DC_array.shape[0], t1 - t0)
                    )

        index = 0
        skipped_write = 0
        event_file = dataio.I3File(a_file)
        for frame in event_file:
            if frame.Stop == icetray.I3Frame.Physics:
                header = frame["I3EventHeader"]

                # Check Header
                if header.sub_event_stream != "InIceSplit":
                    continue
                if float(header.run_id) != header_array[index][0]:
                    print("Run ID is off")
                    continue
                if float(header.sub_run_id) != header_array[index][1]:
                    print("Sub Run ID is off")
                    continue
                if float(header.event_id) != header_array[index][2]:
                    print("Event ID is off")
                    continue

                # Check for multiple triggers or 8 hit flag
                if skip_event[index] == True:
                    skipped_write += 1
                    index += 1
                    continue

                check_overwrite = []
                for network in range(num_variables):
                    factor = model_factor_list[network]
                    model_type = model_type_list[network]
                    prediction = cnn_predictions[network]
                    for check in check_overwrite:
                        assert (
                            check != model_type
                        ), "Rewriting key, need different names"
                    if model_type == "class":
                        key_name = "FLERCNN_prob_track"
                    elif model_type == "muon":
                        key_name = "FLERCNN_prob_muon"
                    elif model_type == "muonL4":
                        key_name = "FLERCNN_prob_muon_v2"
                    elif model_type == "muonV3":
                        key_name = "FLERCNN_prob_muon_v3"
                    else:
                        key_name = "FLERCNN_%s" % model_type

                    if model_type == "vertex" or model_type == "ending":
                        ending = ["_x", "_y", "_z"]
                        for reco_i in range(prediction.shape[-1]):
                            adjusted_prediction = prediction[index][reco_i] * factor
                            key_name_loop = key_name + ending[reco_i]
                            frame[key_name_loop] = dataclasses.I3Double(
                                adjusted_prediction
                            )
                        x = prediction[index][0] * factor
                        y = prediction[index][1] * factor
                        x_origin = 46.290000915527344
                        y_origin = -34.880001068115234
                        r = np.sqrt((x - x_origin) ** 2 + (y - y_origin) ** 2)
                        if model_type == "vertex":
                            frame["FLERCNN_vertex_rho36"] = dataclasses.I3Double(r)
                        if model_type == "ending":
                            frame["FLERCNN_ending_rho36"] = dataclasses.I3Double(r)
                    elif model_type == "nDOM":
                        frame[key_name] = dataclasses.I3Double(prediction[index])
                    else:
                        adjusted_prediction = prediction[index][0] * factor
                        frame[key_name] = dataclasses.I3Double(adjusted_prediction)
                        if model_type == "muon":
                            frame["FLERCNN_prob_nu"] = dataclasses.I3Double(
                                1.0 - adjusted_prediction
                            )
                        if model_type == "muon_v2":
                            frame["FLERCNN_prob_nu_v2"] = dataclasses.I3Double(
                                1.0 - adjusted_prediction
                            )
                        if model_type == "zenith":
                            frame["FLERCNN_coszen"] = dataclasses.I3Double(
                                np.cos(adjusted_prediction)
                            )

                outfile.push(frame)
                index += 1
            else:
                outfile.push(frame)
        print("Removed %i events due to cuts" % skipped_write)
    
    return 0


def write_hdf5_file(
    filename_list,
    output_dir,
    output_name,
    gcd_file=None,
):
    # add hdf5writer to wirte "outfile" to hdf5 file

    for a_file in filename_list:
        infile_path = output_dir + output_name + ".i3.zst"
        outfile_path = output_dir + output_name + ".hdf5"

        # open tray to run through file outfile
        tray = I3Tray()
        tray.context['I3FileStager'] = dataio.get_stagers()
        tray.AddModule('I3Reader','reader',FilenameList=[gcd_file, infile_path],)

        if real_sim:
            # add function to write true data as keys - model dependent case
            tray.AddModule(
                store_mc_true_variables,
                "WriteTrueDataToKeys_Model_Dependent",
                Streams=[icetray.I3Frame.Physics],
            )
        else:
            # add function to write true data as keys - model independent case
            tray.AddModule(
                write_true_data_to_keys,
                "WriteTrueDataToKeys_Model_Independent",
                If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
                )

        # get weights
        tray.AddModule(
            store_weights,
            "store_weights",
        )

        tray.Add(
            I3HDFWriter,
            Output=outfile_path,
            Keys=HDF5_KEYS,
            SubEventStreams=["InIceSplit"],
        )

        tray.Execute()
        tray.Finish()

    return 0


# Construct list of filenames
event_file_names = sorted(glob.glob(input_file))
assert event_file_names, "No files loaded, please check path."
time_start = time.time()
test_write(
    event_file_names,
    model_name_list,
    output_dir,
    output_name,
    model_factor_list=scale_factor_list,
    model_type_list=variable_list,
    gcd_file=gcdfile,
    use_cleaned_pulses=use_cleaned_pulses,
    charge_min=charge_min,
)
if write_hdf5:
    write_hdf5_file(
        event_file_names,
        output_dir,
        output_name,
        gcd_file=gcdfile,
    )
time_end = time.time()
print("Total time: %f" % (time_end - time_start))
