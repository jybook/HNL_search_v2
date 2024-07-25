import h5py
import matplotlib.pyplot as plt
import numpy as np
from modules.plot import *
from modules.utils import *
from modules.DecayWidths import FullWidth
import math
from tqdm import tqdm
import sys, argparse

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

def get_weight_fixedM(inputfile = '', condition = [], U_tau4_sq = 1e-03, mass = 1.0, mode = 'rate'):
	# Open file
	f = h5py.File(inputfile, 'r')
	nFiles = float(inputfile.split('.')[0].split('_')[-1])
	# Get weight
	c = 299792458.
	hbar = 6.582119569E-25
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
parser.add_argument('-v', '--variable', type=str, default='energy', help='Name of the variable (see below)')
parser.add_argument('-m', '--mode', type=str, default='rate', help='rate or event')
args = parser.parse_args()

inputfile = args.inputfile
variable = args.variable
mode = args.mode

f = [h5py.File(inputfile[0], 'r'),h5py.File(inputfile[1], 'r'),h5py.File(inputfile[2], 'r'),h5py.File(inputfile[3], 'r'),h5py.File(inputfile[4], 'r')]
nFiles = [float(inputfile[0].split('.')[0].split('_')[-1]),float(inputfile[1].split('.')[0].split('_')[-1]),float(inputfile[2].split('.')[0].split('_')[-1]),float(inputfile[3].split('.')[0].split('_')[-1]),float(inputfile[4].split('.')[0].split('_')[-1])]
level = inputfile[0].split('.')[0].split('_')[-2]

print(nFiles)

labelY = 'Rate (mHz/bin)'
if variable == 'U_mass': labelY = 'Rate (mHz)'
if mode == 'event':
	labelY = 'Events/bin'
	if '_' in variable: labelY = 'Events'
labelX = ''
isLogX = True
isLogY = True

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
conf0 = (f[0]['millipedeDict']['E_cascade0_40m']+f[0]['millipedeDict']['E_cascade1_40m'])/f[0]['millipedeDict']['E_tot']
conf1 = (f[1]['millipedeDict']['E_cascade0_40m']+f[1]['millipedeDict']['E_cascade1_40m'])/f[1]['millipedeDict']['E_tot']
conf2 = (f[2]['millipedeDict']['E_cascade0_40m']+f[2]['millipedeDict']['E_cascade1_40m'])/f[2]['millipedeDict']['E_tot']
conf3 = (f[3]['millipedeDict']['E_cascade0_40m']+f[3]['millipedeDict']['E_cascade1_40m'])/f[3]['millipedeDict']['E_tot']
conf4 = (f[4]['millipedeDict']['E_cascade0_40m']+f[4]['millipedeDict']['E_cascade1_40m'])/f[4]['millipedeDict']['E_tot']
cdt_conf = [conf0>0.99,conf1>0.99,conf2>0.99,conf3>0.99,conf4>0.99]
# Condition on the neutrino energy
#cdt_nuE = [f[0]['neutrinoDict']['Energy']>100, f[1]['neutrinoDict']['Energy']>100, f[2]['neutrinoDict']['Energy']>100, f[3]['neutrinoDict']['Energy']>100, f[4]['neutrinoDict']['Energy']>100]
# Condition on the reconstructed cascade angle
#cdt_nuE = [f[0]['neutrinoDict']['Energy']>100, f[1]['neutrinoDict']['Energy']>100, f[2]['neutrinoDict']['Energy']>100, f[3]['neutrinoDict']['Energy']>100, f[4]['neutrinoDict']['Energy']>100]
# "No condition" (always True)
#cdt_basic = [np.ones(len(f[0]['neutrinoDict']['Energy']), dtype=bool),
#             np.ones(len(f[1]['neutrinoDict']['Energy']), dtype=bool),
#             np.ones(len(f[2]['neutrinoDict']['Energy']), dtype=bool),
#             np.ones(len(f[3]['neutrinoDict']['Energy']), dtype=bool),
#             np.ones(len(f[4]['leptonDict']['Energy']), dtype=bool)]
# Condition on DeepCore containement
radiusTaupede0 = [np.sqrt(np.power(f[0]['cascade0TaupedeDict']['X'],2)+np.power(f[0]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[1]['cascade0TaupedeDict']['X'],2)+np.power(f[1]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[2]['cascade0TaupedeDict']['X'],2)+np.power(f[2]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[3]['cascade0TaupedeDict']['X'],2)+np.power(f[3]['cascade0TaupedeDict']['Y'],2)),np.sqrt(np.power(f[4]['cascade0TaupedeDict']['X'],2)+np.power(f[4]['cascade0TaupedeDict']['Y'],2))]
radiusTaupede1 = [np.sqrt(np.power(f[0]['cascade1TaupedeDict']['X'],2)+np.power(f[0]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[1]['cascade1TaupedeDict']['X'],2)+np.power(f[1]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[2]['cascade1TaupedeDict']['X'],2)+np.power(f[2]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[3]['cascade1TaupedeDict']['X'],2)+np.power(f[3]['cascade1TaupedeDict']['Y'],2)),np.sqrt(np.power(f[4]['cascade1TaupedeDict']['X'],2)+np.power(f[4]['cascade1TaupedeDict']['Y'],2))]
inDeepCore0 = [(radiusTaupede0[0] < 150) & (f[0]['cascade0TaupedeDict']['Z'] < -150) & (f[0]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[1] < 150) & (f[1]['cascade0TaupedeDict']['Z'] < -150) & (f[1]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[2] < 150) & (f[2]['cascade0TaupedeDict']['Z'] < -150) & (f[2]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[3] < 150) & (f[3]['cascade0TaupedeDict']['Z'] < -150) & (f[3]['cascade0TaupedeDict']['Z'] > -500),(radiusTaupede0[4] < 150) & (f[4]['cascade0TaupedeDict']['Z'] < -150) & (f[4]['cascade0TaupedeDict']['Z'] > -500)]
inDeepCore1 = [(radiusTaupede1[0] < 150) & (f[0]['cascade1TaupedeDict']['Z'] < -150) & (f[0]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[1] < 150) & (f[1]['cascade1TaupedeDict']['Z'] < -150) & (f[1]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[2] < 150) & (f[2]['cascade1TaupedeDict']['Z'] < -150) & (f[2]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[3] < 150) & (f[3]['cascade1TaupedeDict']['Z'] < -150) & (f[3]['cascade1TaupedeDict']['Z'] > -500),(radiusTaupede1[4] < 150) & (f[4]['cascade1TaupedeDict']['Z'] < -150) & (f[4]['cascade1TaupedeDict']['Z'] > -500)]
cdt_DC_Reco = [inDeepCore0[0] & inDeepCore1[0],inDeepCore0[1] & inDeepCore1[1],inDeepCore0[2] & inDeepCore1[2],inDeepCore0[3] & inDeepCore1[3],inDeepCore0[4] & inDeepCore1[4]]
# Condition to get meaningful reconstructed variables 
cdt_basic = [(f[0]['cascade0TaupedeDict']['Energy'] > 0) & (f[0]['cascade1TaupedeDict']['Energy'] > 0) & (f[0]['millipedeDict']['E_cascade0_15m'] > 0) & (f[0]['millipedeDict']['E_cascade1_15m'] > 0) & (f[0]['cascade0TaupedeDict']['Energy'] < 1000) & (f[0]['cascade1TaupedeDict']['Energy'] < 1000) & (f[0]['millipedeDict']['E_tot'] < 2000) & (f[0]['taupedeDict']['bestFitLength'] > 0),
             (f[1]['cascade0TaupedeDict']['Energy'] > 0) & (f[1]['cascade1TaupedeDict']['Energy'] > 0) & (f[1]['millipedeDict']['E_cascade0_15m'] > 0) & (f[1]['millipedeDict']['E_cascade1_15m'] > 0) & (f[1]['cascade0TaupedeDict']['Energy'] < 1000) & (f[1]['cascade1TaupedeDict']['Energy'] < 1000) & (f[1]['millipedeDict']['E_tot'] < 2000) & (f[1]['taupedeDict']['bestFitLength'] > 0),
             (f[2]['cascade0TaupedeDict']['Energy'] > 0) & (f[2]['cascade1TaupedeDict']['Energy'] > 0) & (f[2]['millipedeDict']['E_cascade0_15m'] > 0) & (f[2]['millipedeDict']['E_cascade1_15m'] > 0) & (f[2]['cascade0TaupedeDict']['Energy'] < 1000) & (f[2]['cascade1TaupedeDict']['Energy'] < 1000) & (f[2]['millipedeDict']['E_tot'] < 2000) & (f[2]['taupedeDict']['bestFitLength'] > 0),
             (f[3]['cascade0TaupedeDict']['Energy'] > 0) & (f[3]['cascade1TaupedeDict']['Energy'] > 0) & (f[3]['millipedeDict']['E_cascade0_15m'] > 0) & (f[3]['millipedeDict']['E_cascade1_15m'] > 0) & (f[3]['cascade0TaupedeDict']['Energy'] < 1000) & (f[3]['cascade1TaupedeDict']['Energy'] < 1000) & (f[3]['millipedeDict']['E_tot'] < 2000) & (f[3]['taupedeDict']['bestFitLength'] > 0),
             (f[4]['cascade0TaupedeDict']['Energy'] > 0) & (f[4]['cascade1TaupedeDict']['Energy'] > 0) & (f[4]['millipedeDict']['E_cascade0_15m'] > 0) & (f[4]['millipedeDict']['E_cascade1_15m'] > 0) & (f[4]['cascade0TaupedeDict']['Energy'] < 1000) & (f[4]['cascade1TaupedeDict']['Energy'] < 1000) & (f[4]['millipedeDict']['E_tot'] < 5000) & (f[4]['taupedeDict']['bestFitLength'] > 0)]
# Cascade not nan
cdt_nan = np.logical_not(np.isnan(f[0]['cascade0Dict']['Energy'])) & np.logical_not(np.isnan(f[0]['cascade1Dict']['Energy']))
# Condition on energies
cdt_E = [(f[0]['cascade0TaupedeDict']['Energy'] > 3) & (f[0]['cascade1TaupedeDict']['Energy'] > 3),
             (f[1]['cascade0TaupedeDict']['Energy'] > 3) & (f[1]['cascade1TaupedeDict']['Energy'] > 3),
             (f[2]['cascade0TaupedeDict']['Energy'] > 3) & (f[2]['cascade1TaupedeDict']['Energy'] > 3),
             (f[3]['cascade0TaupedeDict']['Energy'] > 3) & (f[3]['cascade1TaupedeDict']['Energy'] > 3),
             (f[4]['cascade0TaupedeDict']['Energy'] > 3) & (f[4]['cascade1TaupedeDict']['Energy'] > 3)]
# Condition on BDT output
#cdt_bdt = [f[0]['bdtOutput'][:] > -0.55, f[1]['bdtOutput'][:] > -0.55, f[2]['bdtOutput'][:] > -0.55, f[3]['bdtOutput'][:] > -0.55, f[4]['bdtOutput'][:] > -0.55]
# Condition on CC/NC
cdt_CC = [f[1]['leptonDict']['Type']%2!=0,
          f[2]['leptonDict']['Type']%2!=0,
          f[3]['leptonDict']['Type']%2!=0]
cdt_NC_only = [f[1]['leptonDict']['Type']%2==0,
               f[2]['leptonDict']['Type']%2==0,
               f[3]['leptonDict']['Type']%2==0]
# Condition on Taupede best-fit zenith
cdt_zenith = [np.cos(f[0]['cascade0TaupedeDict']['Zenith']) < 0, np.cos(f[1]['cascade0TaupedeDict']['Zenith']) < 0, np.cos(f[2]['cascade0TaupedeDict']['Zenith']) < 0, np.cos(f[3]['cascade0TaupedeDict']['Zenith']) < 0, np.cos(f[4]['cascade0TaupedeDict']['Zenith']) < 0]
# Final condition (signal and CC)
#cdt = [cdt_basic[0] & cdt_DC_Reco[0] & cdt_bdt[0], cdt_basic[1] & cdt_DC_Reco[1] & cdt_bdt[1] & cdt_CC[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_bdt[2] & cdt_CC[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_bdt[3] & cdt_CC[2], cdt_basic[4] & cdt_DC_Reco[4] & cdt_bdt[4]]
cdt = [cdt_basic[0] & cdt_DC_Reco[0], cdt_basic[1] & cdt_DC_Reco[1] & cdt_CC[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_CC[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_CC[2], cdt_basic[4] & cdt_DC_Reco[4]]
#cdt = [cdt_basic[0], cdt_basic[1] & cdt_CC[0], cdt_basic[2] & cdt_CC[1], cdt_basic[3] & cdt_CC[2], cdt_basic[4]]
# Final condition (NC)
#cdt_NC = [cdt_basic[1] & cdt_DC_Reco[1] & cdt_bdt[1] & cdt_NC_only[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_bdt[2] & cdt_NC_only[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_bdt[3] & cdt_NC_only[2]]
cdt_NC = [cdt_basic[1] & cdt_DC_Reco[1] & cdt_NC_only[0], cdt_basic[2] & cdt_DC_Reco[2] & cdt_NC_only[1], cdt_basic[3] & cdt_DC_Reco[3] & cdt_NC_only[2]]
#cdt_NC = [cdt_NC_only[0], cdt_NC_only[1], cdt_NC_only[2]]

print(str(len(f[0]['neutrinoDict']['Energy'][cdt[0]]))+' ('+str(100*len(f[0]['neutrinoDict']['Energy'][cdt[0]])/len(f[0]['neutrinoDict']['Energy']))+'%) signal events pass the condition')
print(str(len(f[1]['neutrinoDict']['Energy'][cdt[1]]))+' ('+str(100*len(f[1]['neutrinoDict']['Energy'][cdt[1]])/len(f[1]['neutrinoDict']['Energy'][cdt_CC[0]]))+'%) nue CC events pass the condition')
print(str(len(f[2]['neutrinoDict']['Energy'][cdt[2]]))+' ('+str(100*len(f[2]['neutrinoDict']['Energy'][cdt[2]])/len(f[2]['neutrinoDict']['Energy'][cdt_CC[1]]))+'%) numu CC events pass the condition')
print(str(len(f[3]['neutrinoDict']['Energy'][cdt[3]]))+' ('+str(100*len(f[3]['neutrinoDict']['Energy'][cdt[3]])/len(f[3]['neutrinoDict']['Energy'][cdt_CC[2]]))+'%) nutau CC events pass the condition')
print(str(len(f[1]['neutrinoDict']['Energy'][cdt_NC[0]])+len(f[2]['neutrinoDict']['Energy'][cdt_NC[1]])+len(f[3]['neutrinoDict']['Energy'][cdt_NC[2]]))+' ('+str(100*(len(f[1]['neutrinoDict']['Energy'][cdt_NC[0]])+len(f[2]['neutrinoDict']['Energy'][cdt_NC[1]])+len(f[3]['neutrinoDict']['Energy'][cdt_NC[2]]))/(len(f[1]['neutrinoDict']['Energy'][cdt_NC_only[0]])+len(f[2]['neutrinoDict']['Energy'][cdt_NC_only[1]])+len(f[3]['neutrinoDict']['Energy'][cdt_NC_only[2]])))+'%) nu NC events pass the condition')
print(str(len(f[4]['leptonDict']['Energy'][cdt[4]]))+' ('+str(100*len(f[4]['leptonDict']['Energy'][cdt[4]])/len(f[4]['leptonDict']['Energy']))+'%) muon events pass the condition')

# Weighting (choose a value for the mixing parameter - does not apply for U-mass 2D plots)
U_tau4_sq = 1e-03
weights = [get_weight(inputfile[0],cdt[0],U_tau4_sq,mode),1000*f[1]['weight']['value'][cdt[1]]/nFiles[1],1000*f[2]['weight']['value'][cdt[2]]/nFiles[2],1000*f[3]['weight']['value'][cdt[3]]/nFiles[3],1000*f[4]['weight']['value'][cdt[4]]/nFiles[4]]
weights_NC = np.concatenate((f[1]['weight']['value'][cdt_NC[0]]/nFiles[1],f[2]['weight']['value'][cdt_NC[1]]/nFiles[2],f[3]['weight']['value'][cdt_NC[2]]/nFiles[3]))
#weights[np.isnan(weights)] = 0
#weights[weights<1e-20] = 0
#weights = np.ones(len(weights))

array1 = []
array2 = []
array3 = []
array4 = []
array5 = []
array6 = []
bins = []

#nNaN = 0
#for fi in f:
#	for e in fi['cascade1Dict']['Energy']:
#		if np.isnan(e): nNaN += 1
#print str(100*float(nNaN)/len(f['cascade1Dict']['Energy']))+' % of NaN'

if variable == 'E0':
	array1 = f[0]['cascade0TaupedeDict']['Energy'][cdt[0]]
	array2 = f[1]['cascade0TaupedeDict']['Energy'][cdt[1]]
	array3 = f[2]['cascade0TaupedeDict']['Energy'][cdt[2]]
	array4 = f[3]['cascade0TaupedeDict']['Energy'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade0TaupedeDict']['Energy'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Energy'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Energy'][cdt_NC[2]]))
	array6 = f[4]['cascade0TaupedeDict']['Energy'][cdt[4]]
	bins = np.logspace(np.log10(3),np.log10(1e+03),31)
	labelX = r'$E_{0}$ (GeV)'
elif variable == 'E1':
	array1 = f[0]['cascade1TaupedeDict']['Energy'][cdt[0]]
	array2 = f[1]['cascade1TaupedeDict']['Energy'][cdt[1]]
	array3 = f[2]['cascade1TaupedeDict']['Energy'][cdt[2]]
	array4 = f[3]['cascade1TaupedeDict']['Energy'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade1TaupedeDict']['Energy'][cdt_NC[0]],f[2]['cascade1TaupedeDict']['Energy'][cdt_NC[1]],f[3]['cascade1TaupedeDict']['Energy'][cdt_NC[2]]))
	array6 = f[4]['cascade1TaupedeDict']['Energy'][cdt[4]]
	bins = np.logspace(np.log10(3),np.log10(1e+03),31)
	labelX = r'$E_{1}$ (GeV)'
#elif variable == 'energyAsym':
#	array1 = (f[0]['cascade0TaupedeDict']['Energy'][cdt[0]]-f[0]['cascade1TaupedeDict']['Energy'][cdt[0]])/(f[0]['cascade0TaupedeDict']['Energy'][cdt[0]]+f[0]['cascade1TaupedeDict']['Energy'][cdt[0]]) 
#	array2 = (f[1]['cascade0TaupedeDict']['Energy'][cdt[1]]-f[1]['cascade1TaupedeDict']['Energy'][cdt[1]])/(f[1]['cascade0TaupedeDict']['Energy'][cdt[1]]+f[1]['cascade1TaupedeDict']['Energy'][cdt[1]])
#	array3 = (f[2]['cascade0TaupedeDict']['Energy'][cdt[2]]-f[2]['cascade1TaupedeDict']['Energy'][cdt[2]])/(f[2]['cascade0TaupedeDict']['Energy'][cdt[2]]+f[2]['cascade1TaupedeDict']['Energy'][cdt[2]])
#	array4 = (f[3]['cascade0TaupedeDict']['Energy'][cdt[3]]-f[3]['cascade1TaupedeDict']['Energy'][cdt[3]])/(f[3]['cascade0TaupedeDict']['Energy'][cdt[3]]+f[3]['cascade1TaupedeDict']['Energy'][cdt[3]])
#	array5_E0 = np.concatenate((f[1]['cascade0TaupedeDict']['Energy'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Energy'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Energy'][cdt_NC[2]]))
#	array5_E1 = np.concatenate((f[1]['cascade1TaupedeDict']['Energy'][cdt_NC[0]],f[2]['cascade1TaupedeDict']['Energy'][cdt_NC[1]],f[3]['cascade1TaupedeDict']['Energy'][cdt_NC[2]]))
#	array5 = (array5_E0-array5_E1)/(array5_E0+array5_E1)
#	bins = np.linspace(-1,1,30)
#	isLogX = False
#	labelX = r'Energy asymmetry'
elif variable == 'energyAsym':
	array1 = (f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]-f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]])/(f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]+f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]]) 
	array2 = (f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]-f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])/(f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]+f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])
	array3 = (f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]-f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])/(f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]+f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])
	array4 = (f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]-f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])/(f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]+f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])
	array5_E0 = np.concatenate((f[1]['millipedeDict']['E_cascade0_15m'][cdt_NC[0]],f[2]['millipedeDict']['E_cascade0_15m'][cdt_NC[1]],f[3]['millipedeDict']['E_cascade0_15m'][cdt_NC[2]]))
	array5_E1 = np.concatenate((f[1]['millipedeDict']['E_cascade1_15m'][cdt_NC[0]],f[2]['millipedeDict']['E_cascade1_15m'][cdt_NC[1]],f[3]['millipedeDict']['E_cascade1_15m'][cdt_NC[2]]))
	array5 = (array5_E0-array5_E1)/(array5_E0+array5_E1)
	array6 = (f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]-f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])/(f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]+f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])
	bins = np.linspace(-1,1,30)
	isLogX = False
	labelX = r'Energy asymmetry'
elif variable == 'energyConf':
	array1 = (f[0]['millipedeDict']['E_cascade0_15m'][cdt[0]]+f[0]['millipedeDict']['E_cascade1_15m'][cdt[0]])/f[0]['millipedeDict']['E_tot'][cdt[0]]
	array2 = (f[1]['millipedeDict']['E_cascade0_15m'][cdt[1]]+f[1]['millipedeDict']['E_cascade1_15m'][cdt[1]])/f[1]['millipedeDict']['E_tot'][cdt[1]]
	array3 = (f[2]['millipedeDict']['E_cascade0_15m'][cdt[2]]+f[2]['millipedeDict']['E_cascade1_15m'][cdt[2]])/f[2]['millipedeDict']['E_tot'][cdt[2]]
	array4 = (f[3]['millipedeDict']['E_cascade0_15m'][cdt[3]]+f[3]['millipedeDict']['E_cascade1_15m'][cdt[3]])/f[3]['millipedeDict']['E_tot'][cdt[3]]
	array5_E0 = np.concatenate((f[1]['millipedeDict']['E_cascade0_15m'][cdt_NC[0]],f[2]['millipedeDict']['E_cascade0_15m'][cdt_NC[1]],f[3]['millipedeDict']['E_cascade0_15m'][cdt_NC[2]]))
	array5_E1 = np.concatenate((f[1]['millipedeDict']['E_cascade1_15m'][cdt_NC[0]],f[2]['millipedeDict']['E_cascade1_15m'][cdt_NC[1]],f[3]['millipedeDict']['E_cascade1_15m'][cdt_NC[2]]))
	array5_E_tot = np.concatenate((f[1]['millipedeDict']['E_tot'][cdt_NC[0]],f[2]['millipedeDict']['E_tot'][cdt_NC[1]],f[3]['millipedeDict']['E_tot'][cdt_NC[2]]))
	array5 = (array5_E0+array5_E1)/(array5_E_tot)
	array6 = (f[4]['millipedeDict']['E_cascade0_15m'][cdt[4]]+f[4]['millipedeDict']['E_cascade1_15m'][cdt[4]])/f[4]['millipedeDict']['E_tot'][cdt[4]]
	bins = np.linspace(0,1,41)
	isLogX = False
	labelX = r'Energy confinement'
elif variable == 'cosTheta0':
	array1 = np.cos(f[0]['cascade0TaupedeDict']['Zenith'][cdt[0]])
	array2 = np.cos(f[1]['cascade0TaupedeDict']['Zenith'][cdt[1]])
	array3 = np.cos(f[2]['cascade0TaupedeDict']['Zenith'][cdt[2]])
	array4 = np.cos(f[3]['cascade0TaupedeDict']['Zenith'][cdt[3]])
	array5 = np.cos(np.concatenate((f[1]['cascade0TaupedeDict']['Zenith'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Zenith'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Zenith'][cdt_NC[2]])))
	array6 = np.cos(f[4]['cascade0TaupedeDict']['Zenith'][cdt[4]])
	bins = np.linspace(-1,0.2,30)
	labelX = r'$\cos(\theta_{0})$'
	isLogX = False
elif variable == 'TaupedeZenithMinusRetroZenith':
	array1 = f[0]['cascade0TaupedeDict']['Zenith'][cdt[0]]-f[0]['retroDict']['zenith'][cdt[0]] 
	array2 = f[1]['cascade0TaupedeDict']['Zenith'][cdt[1]]-f[1]['retroDict']['zenith'][cdt[1]]
	array3 = f[2]['cascade0TaupedeDict']['Zenith'][cdt[2]]-f[2]['retroDict']['zenith'][cdt[2]]
	array4 = f[3]['cascade0TaupedeDict']['Zenith'][cdt[3]]-f[3]['retroDict']['zenith'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade0TaupedeDict']['Zenith'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Zenith'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Zenith'][cdt_NC[2]]))-np.concatenate((f[1]['retroDict']['zenith'][cdt_NC[0]],f[2]['retroDict']['zenith'][cdt_NC[1]],f[3]['retroDict']['zenith'][cdt_NC[2]]))
	bins = np.linspace(-3,3,30)
	labelX = r'$\theta_{\mathrm{Taupede}}-\theta_{\mathrm{retro}}$'
	isLogX = False
elif variable == 'PID':
	array1 = f[0]['retroDict']['PID'][cdt[0]]
	array2 = f[1]['retroDict']['PID_FullSky'][cdt[1]]
	array3 = f[2]['retroDict']['PID_FullSky'][cdt[2]]
	array4 = f[3]['retroDict']['PID_FullSky'][cdt[3]]
	array5 = np.concatenate((f[1]['retroDict']['PID_FullSky'][cdt_NC[0]],f[2]['retroDict']['PID_FullSky'][cdt_NC[1]],f[3]['retroDict']['PID_FullSky'][cdt_NC[2]]))
	array6 = f[4]['retroDict']['PID_FullSky'][cdt[4]]
	bins = np.linspace(0,1,41)
	labelX = r'PID'
	isLogX = False
elif variable == 'bestFitLength':
	array1 = f[0]['taupedeDict']['bestFitLength'][cdt[0]]
	array2 = f[1]['taupedeDict']['bestFitLength'][cdt[1]]
	array3 = f[2]['taupedeDict']['bestFitLength'][cdt[2]]
	array4 = f[3]['taupedeDict']['bestFitLength'][cdt[3]]
	array5 = np.concatenate((f[1]['taupedeDict']['bestFitLength'][cdt_NC[0]],f[2]['taupedeDict']['bestFitLength'][cdt_NC[1]],f[3]['taupedeDict']['bestFitLength'][cdt_NC[2]]))
	array6 = f[4]['taupedeDict']['bestFitLength'][cdt[4]]
	bins = np.linspace(0,400,41)
	labelX = r'Best fit distance (m)'
	isLogX = False
elif variable == 'logL':
	array1 = f[0]['taupedeDict']['LogL'][cdt[0]]
	array2 = f[1]['taupedeDict']['LogL'][cdt[1]]
	array3 = f[2]['taupedeDict']['LogL'][cdt[2]]
	array4 = f[3]['taupedeDict']['LogL'][cdt[3]]
	array5 = np.concatenate((f[1]['taupedeDict']['LogL'][cdt_NC[0]],f[2]['taupedeDict']['LogL'][cdt_NC[1]],f[3]['taupedeDict']['LogL'][cdt_NC[2]]))
	bins = np.logspace(np.log10(1e+02),np.log10(1e+04),30)
	labelX = r'logL'
elif variable == 'rLogL_taupede':
	array1 = f[0]['taupedeDict']['rLogL'][cdt[0]]
	array2 = f[1]['taupedeDict']['rLogL'][cdt[1]]
	array3 = f[2]['taupedeDict']['rLogL'][cdt[2]]
	array4 = f[3]['taupedeDict']['rLogL'][cdt[3]]
	array5 = np.concatenate((f[1]['taupedeDict']['rLogL'][cdt_NC[0]],f[2]['taupedeDict']['rLogL'][cdt_NC[1]],f[3]['taupedeDict']['rLogL'][cdt_NC[2]]))
	array6 = f[4]['taupedeDict']['rLogL'][cdt[4]]
	#bins = np.logspace(np.log10(5),np.log10(50),30)
	bins = np.linspace(0,50,41)
	isLogX = False
	labelX = r'rLogL'
elif variable == 'rLogL_millipede':
	array1 = f[0]['millipedeDict']['rLogL'][cdt[0]]
	array2 = f[1]['millipedeDict']['rLogL'][cdt[1]]
	array3 = f[2]['millipedeDict']['rLogL'][cdt[2]]
	array4 = f[3]['millipedeDict']['rLogL'][cdt[3]]
	array5 = np.concatenate((f[1]['millipedeDict']['rLogL'][cdt_NC[0]],f[2]['millipedeDict']['rLogL'][cdt_NC[1]],f[3]['millipedeDict']['rLogL'][cdt_NC[2]]))
	array6 = f[4]['millipedeDict']['rLogL'][cdt[4]]
	#bins = np.logspace(np.log10(5),np.log10(50),30)
	bins = np.linspace(0,1,41)
	isLogX = False
	labelX = r'rLogL_millipede'
elif variable == 'rLogL_diff':
	array1 = f[0]['taupedeDict']['rLogL'][cdt[0]]-f[0]['millipedeDict']['rLogL'][cdt[0]] 
	array2 = f[1]['taupedeDict']['rLogL'][cdt[1]]-f[1]['millipedeDict']['rLogL'][cdt[1]]
	array3 = f[2]['taupedeDict']['rLogL'][cdt[2]]-f[2]['millipedeDict']['rLogL'][cdt[2]]
	array4 = f[3]['taupedeDict']['rLogL'][cdt[3]]-f[3]['millipedeDict']['rLogL'][cdt[3]]
	array5 = np.concatenate((f[1]['taupedeDict']['rLogL'][cdt_NC[0]],f[2]['taupedeDict']['rLogL'][cdt_NC[1]],f[3]['taupedeDict']['rLogL'][cdt_NC[2]]))-np.concatenate((f[1]['millipedeDict']['rLogL'][cdt_NC[0]],f[2]['millipedeDict']['rLogL'][cdt_NC[1]],f[3]['millipedeDict']['rLogL'][cdt_NC[2]]))
	#bins = np.logspace(np.log10(5),np.log10(50),30)
	bins = np.linspace(0,20,30)
	isLogX = False
	labelX = r'rLogL(taupede)-rLogL(millipede)'
elif variable == 'normedChi2':
	array1 = f[0]['taupedeDict']['chiSquared'][cdt[0]]/f[0]['taupedeDict']['chiSquared_dof'][cdt[0]] 
	array2 = f[1]['taupedeDict']['chiSquared'][cdt[1]]/f[1]['taupedeDict']['chiSquared_dof'][cdt[1]]
	array3 = f[2]['taupedeDict']['chiSquared'][cdt[2]]/f[2]['taupedeDict']['chiSquared_dof'][cdt[2]]
	array4 = f[3]['taupedeDict']['chiSquared'][cdt[3]]/f[3]['taupedeDict']['chiSquared_dof'][cdt[3]]
	array5 = np.concatenate((f[1]['taupedeDict']['chiSquared'][cdt_NC[0]],f[2]['taupedeDict']['chiSquared'][cdt_NC[1]],f[3]['taupedeDict']['chiSquared'][cdt_NC[2]]))/np.concatenate((f[1]['taupedeDict']['chiSquared_dof'][cdt_NC[0]],f[2]['taupedeDict']['chiSquared_dof'][cdt_NC[1]],f[3]['taupedeDict']['chiSquared_dof'][cdt_NC[2]]))
	bins = np.linspace(0,200,30)
	isLogX = False
	labelX = r'$\rm{\chi^{2}/N_{dof}}$'
elif variable == 'normedChi2_diff':
	array1 = (f[0]['taupedeDict']['chiSquared'][cdt[0]]/f[0]['taupedeDict']['chiSquared_dof'][cdt[0]])-(f[0]['millipedeDict']['chiSquared'][cdt[0]]/f[0]['millipedeDict']['chiSquared_dof'][cdt[0]]) 
	array2 = (f[1]['taupedeDict']['chiSquared'][cdt[1]]/f[1]['taupedeDict']['chiSquared_dof'][cdt[1]])-(f[1]['millipedeDict']['chiSquared'][cdt[1]]/f[1]['millipedeDict']['chiSquared_dof'][cdt[1]])
	array3 = (f[2]['taupedeDict']['chiSquared'][cdt[2]]/f[2]['taupedeDict']['chiSquared_dof'][cdt[2]])-(f[2]['millipedeDict']['chiSquared'][cdt[2]]/f[2]['millipedeDict']['chiSquared_dof'][cdt[2]])
	array4 = (f[3]['taupedeDict']['chiSquared'][cdt[3]]/f[3]['taupedeDict']['chiSquared_dof'][cdt[3]])-(f[3]['millipedeDict']['chiSquared'][cdt[3]]/f[3]['millipedeDict']['chiSquared_dof'][cdt[3]])
	array5 = (np.concatenate((f[1]['taupedeDict']['chiSquared'][cdt_NC[0]],f[2]['taupedeDict']['chiSquared'][cdt_NC[1]],f[3]['taupedeDict']['chiSquared'][cdt_NC[2]]))/np.concatenate((f[1]['taupedeDict']['chiSquared_dof'][cdt_NC[0]],f[2]['taupedeDict']['chiSquared_dof'][cdt_NC[1]],f[3]['taupedeDict']['chiSquared_dof'][cdt_NC[2]])))-(np.concatenate((f[1]['millipedeDict']['chiSquared'][cdt_NC[0]],f[2]['millipedeDict']['chiSquared'][cdt_NC[1]],f[3]['millipedeDict']['chiSquared'][cdt_NC[2]]))/np.concatenate((f[1]['millipedeDict']['chiSquared_dof'][cdt_NC[0]],f[2]['millipedeDict']['chiSquared_dof'][cdt_NC[1]],f[3]['millipedeDict']['chiSquared_dof'][cdt_NC[2]])))
	bins = np.linspace(0,200,30)
	isLogX = False
	labelX = r'$\rm{\chi^{2}/N_{dof}}(taupede) - \rm{\chi^{2}/N_{dof}}(millipede)$'
elif variable == 'cascade0X':
	array1 = f[0]['cascade0TaupedeDict']['X'][cdt[0]]
	array2 = f[1]['cascade0TaupedeDict']['X'][cdt[1]]
	array3 = f[2]['cascade0TaupedeDict']['X'][cdt[2]]
	array4 = f[3]['cascade0TaupedeDict']['X'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade0TaupedeDict']['X'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['X'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['X'][cdt_NC[2]]))
	array6 = f[4]['cascade0TaupedeDict']['X'][cdt[4]]
	bins = np.linspace(-600,600,31)
	isLogX = False
	labelX = 'Taupede first cascade X (m)'
elif variable == 'cascade0Y':
	array1 = f[0]['cascade0TaupedeDict']['Y'][cdt[0]]
	array2 = f[1]['cascade0TaupedeDict']['Y'][cdt[1]]
	array3 = f[2]['cascade0TaupedeDict']['Y'][cdt[2]]
	array4 = f[3]['cascade0TaupedeDict']['Y'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade0TaupedeDict']['Y'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Y'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Y'][cdt_NC[2]]))
	array6 = f[4]['cascade0TaupedeDict']['Y'][cdt[4]]
	bins = np.linspace(-600,600,31)
	isLogX = False
	labelX = 'Taupede first cascade Y (m)'
elif variable == 'cascade0Z':
	array1 = f[0]['cascade0TaupedeDict']['Z'][cdt[0]]
	array2 = f[1]['cascade0TaupedeDict']['Z'][cdt[1]]
	array3 = f[2]['cascade0TaupedeDict']['Z'][cdt[2]]
	array4 = f[3]['cascade0TaupedeDict']['Z'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade0TaupedeDict']['Z'][cdt_NC[0]],f[2]['cascade0TaupedeDict']['Z'][cdt_NC[1]],f[3]['cascade0TaupedeDict']['Z'][cdt_NC[2]]))
	array6 = f[4]['cascade0TaupedeDict']['Z'][cdt[4]]
	bins = np.linspace(-500,200,36)
	isLogX = False
	labelX = 'Taupede first cascade Z (m)'
elif variable == 'cascade1Z':
	array1 = f[0]['cascade1TaupedeDict']['Z'][cdt[0]]
	array2 = f[1]['cascade1TaupedeDict']['Z'][cdt[1]]
	array3 = f[2]['cascade1TaupedeDict']['Z'][cdt[2]]
	array4 = f[3]['cascade1TaupedeDict']['Z'][cdt[3]]
	array5 = np.concatenate((f[1]['cascade1TaupedeDict']['Z'][cdt_NC[0]],f[2]['cascade1TaupedeDict']['Z'][cdt_NC[1]],f[3]['cascade1TaupedeDict']['Z'][cdt_NC[2]]))
	array6 = f[4]['cascade1TaupedeDict']['Z'][cdt[4]]
	bins = np.linspace(-500,500,36)
	isLogX = False
	labelX = 'Taupede second cascade Z (m)'
elif variable == 'bdtOutput':
	array1 = f[0]['bdtOutput'][:][cdt[0]] 
	array2 = f[1]['bdtOutput'][:][cdt[1]]
	array3 = f[2]['bdtOutput'][:][cdt[2]]
	array4 = f[3]['bdtOutput'][:][cdt[3]]
	array5 = np.concatenate((f[1]['bdtOutput'][:][cdt_NC[0]],f[2]['bdtOutput'][:][cdt_NC[1]],f[3]['bdtOutput'][:][cdt_NC[2]]))
	array6 = f[4]['bdtOutput'][:][cdt[4]]
	bins = np.linspace(-1,0,51)
	isLogX = False
	labelX = 'BDT output'

if variable != 'U_mass':
	print('Calculating errors...')
	yerr1 = histogram_errors(array1,weights[0],bins)
	print('yerr1 done...')
	yerr2 = histogram_errors(array2,weights[1],bins)
	print('yerr2 done...')
	yerr3 = histogram_errors(array3,weights[2],bins)
	print('yerr3 done...')
	yerr4 = histogram_errors(array4,weights[3],bins)
	print('yerr4 done...')
	yerr5 = histogram_errors(array5,weights_NC,bins)
	print('yerr5 done...')
	yerr6 = histogram_errors(array6,weights[4],bins)
	print('yerr6 done...')
	#plot_5hist_1D(array1,array2,array3,array4,array5,labelX,labelY,'OscNext level '+str(level),bins,weights[0],weights[1],weights[2],weights[3],weights_NC,yerr1,yerr2,yerr3,yerr4,yerr5,isLogX,isLogY,variable)
	plot_6hist_1DStack(array1,array2,array3,array4,array5,array6,labelX,labelY,'OscNext '+str(level),bins,weights[0],weights[1],weights[2],weights[3],weights_NC,weights[4],yerr1,yerr2,yerr3,yerr4,yerr5,yerr6,isLogX,isLogY,variable)
	#plot_5hist_1DStack(array1,array2,array3,array4,array5,labelX,labelY,'OscNext '+str(level),bins,weights[0],weights[1],weights[2],weights[3],weights[4],yerr1,yerr2,yerr3,yerr4,yerr5,isLogX,isLogY,variable)
