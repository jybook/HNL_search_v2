import h5py
import matplotlib.pyplot as plt
import numpy as np
from modules.plot import *
from modules.utils import *
from modules.DecayWidths import FullWidth
import math
from tqdm import tqdm
import sys, argparse
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, roc_auc_score

def compare_train_test(clf, X_train, y_train, X_test, y_test, w_train, w_test, bins=100):
	decisions = []
	for X,y in ((X_train, y_train), (X_test, y_test)):
		d1 = clf.decision_function(X[y>0.5]).ravel()
		d2 = clf.decision_function(X[y<0.5]).ravel()
		decisions += [d1, d2]
	    
	low = min(np.min(d) for d in decisions)
	high = max(np.max(d) for d in decisions)
	low_high = (low,high)
	
	plt.hist(decisions[0], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='S (train)')
	plt.hist(decisions[1], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='B (train)')
	#plt.hist(decisions[0], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, weights=w_train[y_train>0.5], label='S (train)')
	#plt.hist(decisions[1], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, weights=w_train[y_train<0.5], label='B (train)')
	
	#hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True, weights=w_test[y_test>0.5])
	hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
	scale = len(decisions[2]) / sum(hist)
	err = np.sqrt(hist * scale) / scale
	
	width = (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')
	
	#hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True, weights=w_test[y_test<0.5])
	hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
	scale = len(decisions[2]) / sum(hist)
	err = np.sqrt(hist * scale) / scale
	
	plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')
	
	plt.xlabel("BDT output")
	plt.ylabel("Arbitrary units")
	plt.legend(loc='best')

	plt.savefig('classifier_output.pdf')
	plt.savefig('classifier_output.png')
	plt.show()

def weight_lifetime(tau_min,tau_max,lifetime,lifetime_proper):
    pdf_uniform = 1./(tau_max-tau_min)
    pdf_inverse = (1./(np.log(tau_max)-np.log(tau_min)))*(1./lifetime)
    pdf_exp1 = 1./lifetime_proper
    pdf_exp2 = np.exp(-lifetime/lifetime_proper)
    pdf_exp = pdf_exp1*pdf_exp2
    return pdf_exp/pdf_inverse

def get_weight(inputfile = '', condition = [], U_tau4_sq = 1e-03, mode = 'rate'):
    # Open file
    f = h5py.File(inputfile, 'r')
    nFiles = float(inputfile.split('.')[0].split('_')[-1])
    # Get weight
    c = 299792458.
    hbar = 6.582119569E-25
    mass = f['hnlDict']['Mass'][condition]
    gamma = (f['hnlDict']['Energy'][condition]+mass)/mass
    speed = c*np.sqrt(1-np.power(1./gamma,2))
    tau_min = f['propertiesDict']['distanceMin'][condition]/(gamma*speed)
    tau_max = f['propertiesDict']['distanceMax'][condition]/(gamma*speed)
    lifetime = 1e-09*f['propertiesDict']['lifetime'][condition]
    lifetime_proper = hbar/(FullWidth_vec(mass)*U_tau4_sq)
    # First factor needed to account for the fact that only 20% (1/5) of the generated events were processed
    # Second factor needed to account for the fact that some files (over the 1000 initial, this number changes form data set to data set) were lost on the way
    #weights = 5*(1000.0/nFiles)*U_tau4_sq*weight_lifetime(tau_min,tau_max,lifetime,lifetime_proper)*f['weight']['value'][condition]
    weights = (10000.0/nFiles)*U_tau4_sq*weight_lifetime(tau_min,tau_max,lifetime,lifetime_proper)*f['weight']['value'][condition]
    if mode == 'event': weights = 8*365*24*3600*(1./1000.)*weights
    return weights

# Vectorization
weight_lifetime_vec = np.vectorize(weight_lifetime)
FullWidth_vec = np.vectorize(lambda x: FullWidth(x))

# Define input arguments
parser = argparse.ArgumentParser(description='Plot the selected quantity for the selected particle')
parser.add_argument('-i', '--inputfile', nargs='+', type=str, default='190606_L7_995.h5', help='Name of the file')
args = parser.parse_args()

inputfile = args.inputfile

f = [h5py.File(inputfile[0], 'a'),h5py.File(inputfile[1], 'a'),h5py.File(inputfile[2], 'a'),h5py.File(inputfile[3], 'a'),h5py.File(inputfile[4], 'a')]
nFiles = [float(inputfile[0].split('.')[0].split('_')[-1]),float(inputfile[1].split('.')[0].split('_')[-1]),float(inputfile[2].split('.')[0].split('_')[-1]),float(inputfile[3].split('.')[0].split('_')[-1]),float(inputfile[4].split('.')[0].split('_')[-1])]
level = inputfile[0].split('.')[0].split('_')[-2]

### Conditions ###
# "Reconstructible" double cascade condition
distance = f[0]['propertiesDict']['distance']
radius0 = np.sqrt(np.power(f[0]['cascade0Dict']['X'],2)+np.power(f[0]['cascade0Dict']['Y'],2))
radius1 = np.sqrt(np.power(f[0]['cascade1Dict']['X'],2)+np.power(f[0]['cascade1Dict']['Y'],2))
inDeepCore0 = (radius0 < 150) & (f[0]['cascade0Dict']['Z'] < -150) & (f[0]['cascade0Dict']['Z'] > -500)
inIceCube0 = (radius0 < 600) & (f[0]['cascade0Dict']['Z'] < 500) & (f[0]['cascade0Dict']['Z'] > -500)
inDeepCore1 = (radius1 < 150) & (f[0]['cascade1Dict']['Z'] < -150) & (f[0]['cascade1Dict']['Z'] > -500)
inIceCube1 = (radius1 < 600) & (f[0]['cascade1Dict']['Z'] < 500) & (f[0]['cascade1Dict']['Z'] > -500)
# "Reconstructible" double cascade condition
cdt_DC = (distance > 20) & ((inDeepCore0 & (f[0]['cascade0Dict']['Energy'] > 5) & inDeepCore1 & (f[0]['cascade1Dict']['Energy'] > 5)) | (inDeepCore0 & (f[0]['cascade0Dict']['Energy'] > 5) & inIceCube1 & (f[0]['cascade1Dict']['Energy'] > 20)) | (inDeepCore1 & (f[0]['cascade1Dict']['Energy'] > 5) & inIceCube0 & (f[0]['cascade0Dict']['Energy'] > 20)) | (inIceCube0 & (f[0]['cascade0Dict']['Energy'] > 20) & inIceCube1 & (f[0]['cascade1Dict']['Energy'] > 20)))
# Single cascade condition
cdt_SC = (distance < 20) | ((distance > 20) & ((radius0 > 600) | (radius1 > 600) | (f[0]['cascade0Dict']['Z'] < -500) | (f[0]['cascade1Dict']['Z'] < -500) | (f[0]['cascade0Dict']['Z'] > 500) | (f[0]['cascade1Dict']['Z'] > 500)))
# One cascade in I3 condition
cdt_1C = (radius0 > 600) | (radius1 > 600) | (f[0]['cascade0Dict']['Z'] < -500) | (f[0]['cascade1Dict']['Z'] < -500) | (f[0]['cascade0Dict']['Z'] > 500) | (f[0]['cascade1Dict']['Z'] > 500)
# Track condition
cdt_T = f[0]['retroDict']['PID'] > 0.8
# Condition on the confinement
#conf0 = (f[0]['millipedeDict']['E_cascade0_40m']+f[0]['millipedeDict']['E_cascade1_40m'])/f[0]['millipedeDict']['E_tot']
#conf1 = (f[1]['millipedeDict']['E_cascade0_40m']+f[1]['millipedeDict']['E_cascade1_40m'])/f[1]['millipedeDict']['E_tot']
#conf2 = (f[2]['millipedeDict']['E_cascade0_40m']+f[2]['millipedeDict']['E_cascade1_40m'])/f[2]['millipedeDict']['E_tot']
#conf3 = (f[3]['millipedeDict']['E_cascade0_40m']+f[3]['millipedeDict']['E_cascade1_40m'])/f[3]['millipedeDict']['E_tot']
#conf4 = (f[4]['millipedeDict']['E_cascade0_40m']+f[4]['millipedeDict']['E_cascade1_40m'])/f[4]['millipedeDict']['E_tot']
#cdt_conf = [conf0>0.99,conf1>0.99,conf2>0.99,conf3>0.99,conf4>0.99]
# Condition on the neutrino energy
#cdt_nuE = [f[0]['neutrinoDict']['Energy']>100, f[1]['neutrinoDict']['Energy']>100, f[2]['neutrinoDict']['Energy']>100, f[3]['neutrinoDict']['Energy']>100, f[4]['neutrinoDict']['Energy']>100]
# Condition on the reconstructed cascade angle
#cdt_nuE = [f[0]['neutrinoDict']['Energy']>100, f[1]['neutrinoDict']['Energy']>100, f[2]['neutrinoDict']['Energy']>100, f[3]['neutrinoDict']['Energy']>100, f[4]['neutrinoDict']['Energy']>100]
# "No condition" (always True)
#cdt0 = [np.ones(len(f[0]['neutrinoDict']['Energy']), dtype=bool),
#        np.ones(len(f[1]['neutrinoDict']['Energy']), dtype=bool),
#        np.ones(len(f[2]['neutrinoDict']['Energy']), dtype=bool),
#        np.ones(len(f[3]['neutrinoDict']['Energy']), dtype=bool),
#        np.ones(len(f[4]['leptonDict']['Energy']), dtype=bool)]
# Basic condition
cdt_E = [(f[0]['cascade0TaupedeDict']['Energy'] > 5) & (f[0]['cascade1TaupedeDict']['Energy'] > 5),
         (f[1]['cascade0TaupedeDict']['Energy'] > 5) & (f[1]['cascade1TaupedeDict']['Energy'] > 5),
         (f[2]['cascade0TaupedeDict']['Energy'] > 5) & (f[2]['cascade1TaupedeDict']['Energy'] > 5),
         (f[3]['cascade0TaupedeDict']['Energy'] > 5) & (f[3]['cascade1TaupedeDict']['Energy'] > 5),
         (f[4]['cascade0TaupedeDict']['Energy'] > 5) & (f[4]['cascade1TaupedeDict']['Energy'] > 5)]
# Condition on DeepCore containement
radiusTaupede0 = [np.sqrt(np.power(f[0]['cascade0TaupedeDict']['X'],2)+np.power(f[0]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[1]['cascade0TaupedeDict']['X'],2)+np.power(f[1]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[2]['cascade0TaupedeDict']['X'],2)+np.power(f[2]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[3]['cascade0TaupedeDict']['X'],2)+np.power(f[3]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[4]['cascade0TaupedeDict']['X'],2)+np.power(f[4]['cascade0TaupedeDict']['Y'],2))]
radiusTaupede1 = [np.sqrt(np.power(f[0]['cascade1TaupedeDict']['X'],2)+np.power(f[0]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[1]['cascade1TaupedeDict']['X'],2)+np.power(f[1]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[2]['cascade1TaupedeDict']['X'],2)+np.power(f[2]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[3]['cascade1TaupedeDict']['X'],2)+np.power(f[3]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[4]['cascade1TaupedeDict']['X'],2)+np.power(f[4]['cascade1TaupedeDict']['Y'],2))]
inDeepCore0 = [(radiusTaupede0[0] < 150) & (f[0]['cascade0TaupedeDict']['Z'] < -150) & (f[0]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[1] < 150) & (f[1]['cascade0TaupedeDict']['Z'] < -150) & (f[1]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[2] < 150) & (f[2]['cascade0TaupedeDict']['Z'] < -150) & (f[2]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[3] < 150) & (f[3]['cascade0TaupedeDict']['Z'] < -150) & (f[3]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[4] < 150) & (f[4]['cascade0TaupedeDict']['Z'] < -150) & (f[4]['cascade0TaupedeDict']['Z'] > -500)]
inDeepCore1 = [(radiusTaupede1[0] < 150) & (f[0]['cascade1TaupedeDict']['Z'] < -150) & (f[0]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[1] < 150) & (f[1]['cascade1TaupedeDict']['Z'] < -150) & (f[1]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[2] < 150) & (f[2]['cascade1TaupedeDict']['Z'] < -150) & (f[2]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[3] < 150) & (f[3]['cascade1TaupedeDict']['Z'] < -150) & (f[3]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[4] < 150) & (f[4]['cascade1TaupedeDict']['Z'] < -150) & (f[4]['cascade1TaupedeDict']['Z'] > -500)]
cdt_DC_Reco = [inDeepCore0[0] & inDeepCore1[0], inDeepCore0[1] & inDeepCore1[1], inDeepCore0[2] & inDeepCore1[2], inDeepCore0[3] & inDeepCore1[3], inDeepCore0[4] & inDeepCore1[4]]
# Condition to get meaningful reconstructed variables 
cdt_basic = [(f[0]['cascade0TaupedeDict']['Energy'] > 0) & (f[0]['cascade1TaupedeDict']['Energy'] > 0) & (f[0]['millipedeDict']['E_cascade0_15m'] > 0) & (f[0]['millipedeDict']['E_cascade1_15m'] > 0) & (f[0]['cascade0TaupedeDict']['Energy'] < 1000) & (f[0]['cascade1TaupedeDict']['Energy'] < 5000) & (f[0]['millipedeDict']['E_tot'] < 10000) & (f[0]['taupedeDict']['bestFitLength'] > 0),
             (f[1]['cascade0TaupedeDict']['Energy'] > 0) & (f[1]['cascade1TaupedeDict']['Energy'] > 0) & (f[1]['millipedeDict']['E_cascade0_15m'] > 0) & (f[1]['millipedeDict']['E_cascade1_15m'] > 0) & (f[1]['cascade0TaupedeDict']['Energy'] < 1000) & (f[1]['cascade1TaupedeDict']['Energy'] < 5000) & (f[1]['millipedeDict']['E_tot'] < 10000) & (f[1]['taupedeDict']['bestFitLength'] > 0),
             (f[2]['cascade0TaupedeDict']['Energy'] > 0) & (f[2]['cascade1TaupedeDict']['Energy'] > 0) & (f[2]['millipedeDict']['E_cascade0_15m'] > 0) & (f[2]['millipedeDict']['E_cascade1_15m'] > 0) & (f[2]['cascade0TaupedeDict']['Energy'] < 1000) & (f[2]['cascade1TaupedeDict']['Energy'] < 5000) & (f[2]['millipedeDict']['E_tot'] < 10000) & (f[2]['taupedeDict']['bestFitLength'] > 0),
             (f[3]['cascade0TaupedeDict']['Energy'] > 0) & (f[3]['cascade1TaupedeDict']['Energy'] > 0) & (f[3]['millipedeDict']['E_cascade0_15m'] > 0) & (f[3]['millipedeDict']['E_cascade1_15m'] > 0) & (f[3]['cascade0TaupedeDict']['Energy'] < 1000) & (f[3]['cascade1TaupedeDict']['Energy'] < 5000) & (f[3]['millipedeDict']['E_tot'] < 10000) & (f[3]['taupedeDict']['bestFitLength'] > 0),
             (f[4]['cascade0TaupedeDict']['Energy'] > 0) & (f[4]['cascade1TaupedeDict']['Energy'] > 0) & (f[4]['millipedeDict']['E_cascade0_15m'] > 0) & (f[4]['millipedeDict']['E_cascade1_15m'] > 0) & (f[4]['cascade0TaupedeDict']['Energy'] < 1000) & (f[4]['cascade1TaupedeDict']['Energy'] < 5000) & (f[4]['millipedeDict']['E_tot'] < 10000) & (f[4]['taupedeDict']['bestFitLength'] > 0)]
# Cascade not nan
cdt_nan = np.logical_not(np.isnan(f[0]['cascade0Dict']['Energy'])) & np.logical_not(np.isnan(f[0]['cascade1Dict']['Energy']))
# Condition on CC/NC
cdt_CC = [f[1]['leptonDict']['Type']%2!=0,
          f[2]['leptonDict']['Type']%2!=0,
          f[3]['leptonDict']['Type']%2!=0]
cdt_NC_only = [f[1]['leptonDict']['Type']%2==0,
               f[2]['leptonDict']['Type']%2==0,
               f[3]['leptonDict']['Type']%2==0]
# Final condition (signal and CC)
cdt = [cdt_basic[0] & cdt_DC_Reco[0] & cdt_E[0], cdt_basic[1] & cdt_DC_Reco[1] & cdt_E[1], cdt_basic[2] & cdt_DC_Reco[2] & cdt_E[2], cdt_basic[3] & cdt_DC_Reco[3] & cdt_E[3], cdt_basic[4] & cdt_DC_Reco[4] & cdt_E[4]]
#cdt = [cdt_basic[0] & cdt_DC_Reco[0], cdt_basic[1] & cdt_DC_Reco[1], cdt_basic[2] & cdt_DC_Reco[2], cdt_basic[3] & cdt_DC_Reco[3], cdt_basic[4] & cdt_DC_Reco[4]]
#cdt = [cdt_DC & cdt_conf[0],cdt_basic[1] & cdt_conf[1] & cdt_CC[0],cdt_basic[2] & cdt_conf[2] & cdt_CC[1],cdt_basic[3] & cdt_conf[3] & cdt_CC[2]]
# Final condition (NC)
cdt_NC = [cdt_basic[1] & cdt_DC_Reco[1] & cdt_E[1] & cdt_NC_only[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_E[2] & cdt_NC_only[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_E[3] & cdt_NC_only[2]]
#cdt_NC = [cdt_basic[1] & cdt_DC_Reco[1] & cdt_NC_only[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_NC_only[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_NC_only[2]]
#cdt_NC = [cdt_basic[1] & cdt_conf[1] & cdt_NC_only[0], cdt_basic[2] & cdt_conf[2] & cdt_NC_only[1], cdt_basic[3] & cdt_conf[3] & cdt_NC_only[2]]

# Weighting (choose a value for the mixing parameter - does not apply for U-mass 2D plots)
U_tau4_sq = 1e-03
weights = [get_weight(inputfile[0],cdt[0],U_tau4_sq,''),1000*f[1]['weight']['value'][cdt[1]]/nFiles[1],1000*f[2]['weight']['value'][cdt[2]]/nFiles[2],1000*f[3]['weight']['value'][cdt[3]]/nFiles[3],1000*f[4]['weight']['value'][cdt[4]]/nFiles[4]]

print('Prepare datasets for training and evaluation...')

# Defining discriminating variables to be used
# Signal for training
costheta_S        = np.cos(f[0]['cascade0TaupedeDict']['Zenith'][cdt[0]])
L_S               = f[0]['taupedeDict']['bestFitLength'][cdt[0]]
conf15_S          = (f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]+f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]])/f[0]['millipedeDict']['E_tot'][cdt[0]]
asym15_S          = (f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]-f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]])/(f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]+f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]])
PID_S             = f[0]['retroDict']['PID'][cdt[0]]
rLogL_taupede_S   = f[0]['taupedeDict']['rLogL'][cdt[0]]
rLogL_millipede_S = f[0]['millipedeDict']['rLogL'][cdt[0]]
# Full Signal
F_costheta_S        = np.cos(f[0]['cascade0TaupedeDict']['Zenith'])
F_L_S               = f[0]['taupedeDict']['bestFitLength']
F_conf15_S          = np.divide(f[0]['millipedeDict']['E_cascade0_15m']+f[0]['millipedeDict']['E_cascade1_15m'],f[0]['millipedeDict']['E_tot'],out=np.zeros_like(f[0]['millipedeDict']['E_cascade0_15m']+f[0]['millipedeDict']['E_cascade1_15m']),where=f[0]['millipedeDict']['E_tot']!=0)
F_asym15_S          = np.divide(f[0]['millipedeDict']['E_cascade0_15m']-f[0]['millipedeDict']['E_cascade1_15m'],f[0]['millipedeDict']['E_cascade0_15m']+f[0]['millipedeDict']['E_cascade1_15m'],out=np.zeros_like(f[0]['millipedeDict']['E_cascade0_15m']-f[0]['millipedeDict']['E_cascade1_15m']),where=(f[0]['millipedeDict']['E_cascade0_15m']+f[0]['millipedeDict']['E_cascade1_15m'])!=0)
F_PID_S             = f[0]['retroDict']['PID']
F_rLogL_taupede_S   = f[0]['taupedeDict']['rLogL']
F_rLogL_millipede_S = f[0]['millipedeDict']['rLogL']
# Backgrounds for training
costheta_nue        = np.cos(f[1]['cascade0TaupedeDict']['Zenith'][cdt[1]])
L_nue               = f[1]['taupedeDict']['bestFitLength'][cdt[1]]
conf15_nue          = (f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]+f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])/f[1]['millipedeDict']['E_tot'][cdt[1]]
asym15_nue          = (f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]-f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])/(f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]+f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])
PID_nue             = f[1]['retroDict']['PID_FullSky'][cdt[1]]
rLogL_taupede_nue   = f[1]['taupedeDict']['rLogL'][cdt[1]]
rLogL_millipede_nue = f[1]['millipedeDict']['rLogL'][cdt[1]]
costheta_numu        = np.cos(f[2]['cascade0TaupedeDict']['Zenith'][cdt[2]])
L_numu               = f[2]['taupedeDict']['bestFitLength'][cdt[2]]
conf15_numu          = (f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]+f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])/f[2]['millipedeDict']['E_tot'][cdt[2]]
asym15_numu          = (f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]-f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])/(f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]+f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])
PID_numu             = f[2]['retroDict']['PID_FullSky'][cdt[2]]
rLogL_taupede_numu   = f[2]['taupedeDict']['rLogL'][cdt[2]]
rLogL_millipede_numu = f[2]['millipedeDict']['rLogL'][cdt[2]]
costheta_nutau        = np.cos(f[3]['cascade0TaupedeDict']['Zenith'][cdt[3]])
L_nutau               = f[3]['taupedeDict']['bestFitLength'][cdt[3]]
conf15_nutau          = (f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]+f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])/f[3]['millipedeDict']['E_tot'][cdt[3]]
asym15_nutau          = (f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]-f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])/(f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]+f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])
PID_nutau             = f[3]['retroDict']['PID_FullSky'][cdt[3]]
rLogL_taupede_nutau   = f[3]['taupedeDict']['rLogL'][cdt[3]]
rLogL_millipede_nutau = f[3]['millipedeDict']['rLogL'][cdt[3]]
costheta_mu        = np.cos(f[4]['cascade0TaupedeDict']['Zenith'][cdt[4]])
L_mu               = f[4]['taupedeDict']['bestFitLength'][cdt[4]]
conf15_mu          = (f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]+f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])/f[4]['millipedeDict']['E_tot'][cdt[4]]
asym15_mu          = (f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]-f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])/(f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]+f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])
PID_mu             = f[4]['retroDict']['PID_FullSky'][cdt[4]]
rLogL_taupede_mu   = f[4]['taupedeDict']['rLogL'][cdt[4]]
rLogL_millipede_mu = f[4]['millipedeDict']['rLogL'][cdt[4]]
# Full Backgrounds
F_costheta_nue        = np.cos(f[1]['cascade0TaupedeDict']['Zenith'])
F_L_nue               = f[1]['taupedeDict']['bestFitLength']
F_conf15_nue          = np.divide(f[1]['millipedeDict']['E_cascade0_15m']+f[1]['millipedeDict']['E_cascade1_15m'],f[1]['millipedeDict']['E_tot'],out=np.zeros_like(f[1]['millipedeDict']['E_cascade0_15m']+f[1]['millipedeDict']['E_cascade1_15m']),where=f[1]['millipedeDict']['E_tot']!=0)
F_asym15_nue          = np.divide(f[1]['millipedeDict']['E_cascade0_15m']-f[1]['millipedeDict']['E_cascade1_15m'],f[1]['millipedeDict']['E_cascade0_15m']+f[1]['millipedeDict']['E_cascade1_15m'],out=np.zeros_like(f[1]['millipedeDict']['E_cascade0_15m']-f[1]['millipedeDict']['E_cascade1_15m']),where=(f[1]['millipedeDict']['E_cascade0_15m']+f[1]['millipedeDict']['E_cascade1_15m'])!=0)
F_PID_nue             = f[1]['retroDict']['PID_FullSky']
F_rLogL_taupede_nue   = f[1]['taupedeDict']['rLogL']
F_rLogL_millipede_nue = f[1]['millipedeDict']['rLogL']
F_costheta_numu        = np.cos(f[2]['cascade0TaupedeDict']['Zenith'])
F_L_numu               = f[2]['taupedeDict']['bestFitLength']
F_conf15_numu          = np.divide(f[2]['millipedeDict']['E_cascade0_15m']+f[2]['millipedeDict']['E_cascade1_15m'],f[2]['millipedeDict']['E_tot'],out=np.zeros_like(f[2]['millipedeDict']['E_cascade0_15m']+f[2]['millipedeDict']['E_cascade1_15m']),where=f[2]['millipedeDict']['E_tot']!=0)
F_asym15_numu          = np.divide(f[2]['millipedeDict']['E_cascade0_15m']-f[2]['millipedeDict']['E_cascade1_15m'],f[2]['millipedeDict']['E_cascade0_15m']+f[2]['millipedeDict']['E_cascade1_15m'],out=np.zeros_like(f[2]['millipedeDict']['E_cascade0_15m']-f[2]['millipedeDict']['E_cascade1_15m']),where=(f[2]['millipedeDict']['E_cascade0_15m']+f[2]['millipedeDict']['E_cascade1_15m'])!=0)
F_PID_numu             = f[2]['retroDict']['PID_FullSky']
F_rLogL_taupede_numu   = f[2]['taupedeDict']['rLogL']
F_rLogL_millipede_numu = f[2]['millipedeDict']['rLogL']
F_costheta_nutau        = np.cos(f[3]['cascade0TaupedeDict']['Zenith'])
F_L_nutau               = f[3]['taupedeDict']['bestFitLength']
F_conf15_nutau          = np.divide(f[3]['millipedeDict']['E_cascade0_15m']+f[3]['millipedeDict']['E_cascade1_15m'],f[3]['millipedeDict']['E_tot'],out=np.zeros_like(f[3]['millipedeDict']['E_cascade0_15m']+f[3]['millipedeDict']['E_cascade1_15m']),where=f[3]['millipedeDict']['E_tot']!=0)
F_asym15_nutau          = np.divide(f[3]['millipedeDict']['E_cascade0_15m']-f[3]['millipedeDict']['E_cascade1_15m'],f[3]['millipedeDict']['E_cascade0_15m']+f[3]['millipedeDict']['E_cascade1_15m'],out=np.zeros_like(f[3]['millipedeDict']['E_cascade0_15m']-f[3]['millipedeDict']['E_cascade1_15m']),where=(f[3]['millipedeDict']['E_cascade0_15m']+f[3]['millipedeDict']['E_cascade1_15m'])!=0)
F_PID_nutau             = f[3]['retroDict']['PID_FullSky']
F_rLogL_taupede_nutau   = f[3]['taupedeDict']['rLogL']
F_rLogL_millipede_nutau = f[3]['millipedeDict']['rLogL']
F_costheta_mu        = np.cos(f[4]['cascade0TaupedeDict']['Zenith'])
F_L_mu               = f[4]['taupedeDict']['bestFitLength']
F_conf15_mu          = np.divide(f[4]['millipedeDict']['E_cascade0_15m']+f[4]['millipedeDict']['E_cascade1_15m'],f[4]['millipedeDict']['E_tot'],out=np.zeros_like(f[4]['millipedeDict']['E_cascade0_15m']+f[4]['millipedeDict']['E_cascade1_15m']),where=f[4]['millipedeDict']['E_tot']!=0)
F_asym15_mu          = np.divide(f[4]['millipedeDict']['E_cascade0_15m']-f[4]['millipedeDict']['E_cascade1_15m'],f[4]['millipedeDict']['E_cascade0_15m']+f[4]['millipedeDict']['E_cascade1_15m'],out=np.zeros_like(f[4]['millipedeDict']['E_cascade0_15m']-f[4]['millipedeDict']['E_cascade1_15m']),where=(f[4]['millipedeDict']['E_cascade0_15m']+f[4]['millipedeDict']['E_cascade1_15m'])!=0)
F_PID_mu             = f[4]['retroDict']['PID_FullSky']
F_rLogL_taupede_mu   = f[4]['taupedeDict']['rLogL']
F_rLogL_millipede_mu = f[4]['millipedeDict']['rLogL']

print('Concatenate and shape data sets...')

costheta_B        = np.concatenate((costheta_nue,costheta_numu,costheta_nutau,costheta_mu))
L_B               = np.concatenate((L_nue,L_numu,L_nutau,L_mu))
conf15_B          = np.concatenate((conf15_nue,conf15_numu,conf15_nutau,conf15_mu))
asym15_B          = np.concatenate((asym15_nue,asym15_numu,asym15_nutau,asym15_mu))
PID_B             = np.concatenate((PID_nue,PID_numu,PID_nutau,PID_mu))
rLogL_taupede_B   = np.concatenate((rLogL_taupede_nue,rLogL_taupede_numu,rLogL_taupede_nutau,rLogL_taupede_mu))
rLogL_millipede_B = np.concatenate((rLogL_millipede_nue,rLogL_millipede_numu,rLogL_millipede_nutau,rLogL_millipede_mu))

F_costheta_B        = np.concatenate((F_costheta_nue,F_costheta_numu,F_costheta_nutau,F_costheta_mu))
F_L_B               = np.concatenate((F_L_nue,F_L_numu,F_L_nutau,F_L_mu))
F_conf15_B          = np.concatenate((F_conf15_nue,F_conf15_numu,F_conf15_nutau,F_conf15_mu))
F_asym15_B          = np.concatenate((F_asym15_nue,F_asym15_numu,F_asym15_nutau,F_asym15_mu))
F_PID_B             = np.concatenate((F_PID_nue,F_PID_numu,F_PID_nutau,F_PID_mu))
F_rLogL_taupede_B   = np.concatenate((F_rLogL_taupede_nue,F_rLogL_taupede_numu,F_rLogL_taupede_nutau,F_rLogL_taupede_mu))
F_rLogL_millipede_B = np.concatenate((F_rLogL_millipede_nue,F_rLogL_millipede_numu,F_rLogL_millipede_nutau,F_rLogL_millipede_mu))

signal = np.column_stack((costheta_S,L_S,conf15_S,asym15_S,PID_S,rLogL_taupede_S,rLogL_millipede_S))
bg     = np.column_stack((costheta_B,L_B,conf15_B,asym15_B,PID_B,rLogL_taupede_B,rLogL_millipede_B))
#bg     = np.column_stack((costheta_nue,L_nue,conf15_nue,asym15_nue,PID_nue,rLogL_taupede_nue,rLogL_millipede_nue))
X = np.concatenate((signal,bg))
y = np.concatenate((np.ones(signal.shape[0]),np.zeros(bg.shape[0])))
w = np.concatenate((weights[0],weights[1],weights[2],weights[3],weights[4]))
#w = np.concatenate((weights[0],weights[1]))

print('Split data sets between train and test...')

X_dev,X_eval, y_dev,y_eval, w_dev,w_eval = train_test_split(X, y, w, test_size=0.5, random_state=42)
X_train,X_test, y_train,y_test, w_train,w_test = train_test_split(X_dev, y_dev, w_dev, test_size=0.5, random_state=492)

print('Declare the classifier object...')
clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3), algorithm="SAMME", n_estimators=200, random_state=0)
#clf.fit(X_train, y_train)
print('Fit over the defined data sets...')
clf.fit(X_train, y_train, sample_weight=w_train)

print('Plot the results...')
compare_train_test(clf, X_train, y_train, X_test, y_test, w_train, w_test)

# Add output to signal and BG files
print('Add the classifier to the files...')
signal = np.column_stack((F_costheta_S,F_L_S,F_conf15_S,F_asym15_S,F_PID_S,F_rLogL_taupede_S,F_rLogL_millipede_S))
bg     = np.column_stack((F_costheta_B,F_L_B,F_conf15_B,F_asym15_B,F_PID_B,F_rLogL_taupede_B,F_rLogL_millipede_B))
X      = np.concatenate((signal,bg))
y      = np.concatenate((np.ones(signal.shape[0]),np.full(len(F_costheta_nue),2),np.full(len(F_costheta_numu),3),np.full(len(F_costheta_nutau),4),np.full(len(F_costheta_mu),5)))
bdtOutput_S = clf.decision_function(X[y==1]).ravel()
bdtOutput_nue = clf.decision_function(X[y==2]).ravel()
bdtOutput_numu = clf.decision_function(X[y==3]).ravel()
bdtOutput_nutau = clf.decision_function(X[y==4]).ravel()
bdtOutput_mu = clf.decision_function(X[y==5]).ravel()

if 'bdtOutput' in f[0]: del f[0]['bdtOutput']
if 'bdtOutput' in f[1]: del f[1]['bdtOutput']
if 'bdtOutput' in f[2]: del f[2]['bdtOutput']
if 'bdtOutput' in f[3]: del f[3]['bdtOutput']
if 'bdtOutput' in f[4]: del f[4]['bdtOutput']
f[0].create_dataset('bdtOutput', data=bdtOutput_S)
f[1].create_dataset('bdtOutput', data=bdtOutput_nue)
f[2].create_dataset('bdtOutput', data=bdtOutput_numu)
f[3].create_dataset('bdtOutput', data=bdtOutput_nutau)
f[4].create_dataset('bdtOutput', data=bdtOutput_mu)
f[0].close()
f[1].close()
f[2].close()
f[3].close()
f[4].close()
