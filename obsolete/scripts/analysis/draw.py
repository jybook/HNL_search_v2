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
parser.add_argument('-i', '--inputfile', type=str, default='190606_L7_995.h5', help='Name of the file')
parser.add_argument('-p', '--particle', type=str, default='neutrino', help='Name of the particle (neutrino, hnl, cascade0, cascade1)')
parser.add_argument('-v', '--variable', type=str, default='energy', help='Name of the variable (see below)')
parser.add_argument('-m', '--mode', type=str, default='rate', help='rate or event')
args = parser.parse_args()

inputfile = args.inputfile
particle = args.particle
variable = args.variable
mode = args.mode

f = h5py.File(inputfile, 'r')
nFiles = float(inputfile.split('.')[0].split('_')[-1])
level = inputfile.split('.')[0].split('_')[-2]

labelY = 'Rate (mHz/bin)'
if variable == 'U_mass': labelY = 'Rate (mHz)'
if mode == 'event':
	labelY = 'Events/bin'
	if '_' in variable: labelY = 'Events (8 years)'
labelX = ''
isLogX = True
isLogY = True

# Kinematic cut: uncomment the desired condition
c = 299792458.
mass = f['hnlDict']['Mass']
gamma = (f['hnlDict']['Energy']+mass)/mass
speed = c*np.sqrt(1-np.power(1./gamma,2))
distance = f['propertiesDict']['distance']
radius0 = np.sqrt(np.power(f['cascade0Dict']['X'],2)+np.power(f['cascade0Dict']['Y'],2))
radius1 = np.sqrt(np.power(f['cascade1Dict']['X'],2)+np.power(f['cascade1Dict']['Y'],2))
inDeepCore0 = (radius0 < 150) & (f['cascade0Dict']['Z'] < -150) & (f['cascade0Dict']['Z'] > -500)
inIceCube0 = (radius0 < 600) & (f['cascade0Dict']['Z'] < 500) & (f['cascade0Dict']['Z'] > -500)
inDeepCore1 = (radius1 < 150) & (f['cascade1Dict']['Z'] < -150) & (f['cascade1Dict']['Z'] > -500)
inIceCube1 = (radius1 < 600) & (f['cascade1Dict']['Z'] < 500) & (f['cascade1Dict']['Z'] > -500)
# "Reconstructible" double cascade condition
cdt_DC = (distance > 20) & ((inDeepCore0 & (f['cascade0Dict']['Energy'] > 5) & inDeepCore1 & (f['cascade1Dict']['Energy'] > 5)) | (inDeepCore0 & (f['cascade0Dict']['Energy'] > 5) & inIceCube1 & (f['cascade1Dict']['Energy'] > 20)) | (inDeepCore1 & (f['cascade1Dict']['Energy'] > 5) & inIceCube0 & (f['cascade0Dict']['Energy'] > 20)) | (inIceCube0 & (f['cascade0Dict']['Energy'] > 20) & inIceCube1 & (f['cascade1Dict']['Energy'] > 20)))
# Single cascade condition
cdt_SC = (distance < 20) | ((distance > 20) & ((radius0 > 600) | (radius1 > 600) | (f['cascade0Dict']['Z'] < -500) | (f['cascade1Dict']['Z'] < -500) | (f['cascade0Dict']['Z'] > 500) | (f['cascade1Dict']['Z'] > 500)))
# One cascade in I3 condition
cdt_1C = (radius0 > 600) | (radius1 > 600) | (f['cascade0Dict']['Z'] < -500) | (f['cascade1Dict']['Z'] < -500) | (f['cascade0Dict']['Z'] > 500) | (f['cascade1Dict']['Z'] > 500)
# Track condition
cdt_T = f['retroDict']['PID'] > 0.8
# Cascade energies condition
cdt_E = (f['cascade0Dict']['Energy'] > 5) & (f['cascade1Dict']['Energy'] > 5)
# Condition on the mass range
#condition = (f['hnlDict']['Mass'] > 1.7) & (f['hnlDict']['Mass'] < 1.8)
# This is "no condition" (always True)
cdt0 = np.ones(len(f['neutrinoDict']['Energy']), dtype=bool)
# Final condition
condition = cdt_E

print(str(len(f['neutrinoDict']['Energy'][cdt_DC]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_DC]))/len(f['neutrinoDict']['Energy']))+'%) pass the DC condition')
print(str(len(f['neutrinoDict']['Energy'][cdt_SC]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_SC]))/len(f['neutrinoDict']['Energy']))+'%) pass the SC condition')
print(str(len(f['neutrinoDict']['Energy'][cdt_1C]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_1C]))/len(f['neutrinoDict']['Energy']))+'%) pass the 1C condition')
print(str(len(f['neutrinoDict']['Energy'][cdt_DC & cdt_T]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_DC & cdt_T]))/len(f['neutrinoDict']['Energy']))+'%) pass the DC and the T conditions')
print(str(len(f['neutrinoDict']['Energy'][cdt_SC & cdt_T]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_SC & cdt_T]))/len(f['neutrinoDict']['Energy']))+'%) pass the SC and the T conditions')
print(str(len(f['neutrinoDict']['Energy'][condition]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][condition]))/len(f['neutrinoDict']['Energy']))+'%) pass the final condition')

# Weighting (choose a value for the mixing parameter - does not apply for U-mass 2D plots)
U_tau4_sq = 1e-03
weights = get_weight(inputfile,condition,U_tau4_sq,mode)
weights[np.isnan(weights)] = 0
#weights[weights<1e-20] = 0
weights = np.ones(len(weights))

array = []
bins = []

nNaN = 0
for e in f['cascade1Dict']['Energy']:
	if np.isnan(e): nNaN += 1
print(str(100*float(nNaN)/len(f['cascade1Dict']['Energy']))+' % of NaN')

if variable == 'mass':
	array = f['hnlDict']['Mass'][condition]
	bins = np.logspace(np.log10(0.1),np.log10(3.),30)
	labelX = 'HNL mass'
elif variable == 'energy':
	array = f[particle+'Dict']['Energy'][condition]
	if particle== 'neutrino': bins = np.logspace(np.log10(2.),np.log10(1e+04),30)
	else: bins = np.logspace(np.log10(1e-03),np.log10(1e+04),30)
	labelX = 'Energy (GeV)'
elif variable == 'zenith':
	array = np.cos(f[particle+'Dict']['Zenith'][condition])
	bins = np.linspace(-1.,0.2,30)
	labelX = r'$\cos(\theta)$'
	isLogX = False
elif variable == 'radius':
	array = np.sqrt(np.power(f[particle+'Dict']['X'][condition],2)+np.power(f[particle+'Dict']['Y'][condition],2))
	bins = np.linspace(0,1000,30)
	labelX = 'Radius (m)'
	isLogX = False
elif variable == 'z':
	array = f[particle+'Dict']['Z'][condition]
	bins = np.linspace(-700,-100,30)
	labelX = 'Z (m)'
	isLogX = False
elif variable == 'lifetime':
	array = f['propertiesDict']['lifetime'][condition]
	bins = np.logspace(np.log10(1e-02),np.log10(1e+05),30)
	labelX = 'Lifetime (ns)'
elif variable == 'distance':
	array = f['propertiesDict']['distance'][condition]
	bins = np.logspace(np.log10(1e-01),np.log10(1000),30)
	labelX = 'Distance between cascades (m)'
elif variable == 'bestFitLength':
	array = f['taupedeDict']['bestFitLength'][condition]
	bins = np.logspace(np.log10(1e-01),np.log10(1000),30)
	labelX = 'Taupede best fit distance between cascades (m)'
elif variable == 'distanceResolution':
	array = (f['propertiesDict']['distance'][condition]-f['taupedeDict']['bestFitLength'][condition])/f['propertiesDict']['distance'][condition]
	bins = np.linspace(-10,10,30)
	labelX = r'$\mathrm{(L_{True}-L_{Best \ fit})/L_{True}}$'
	isLogX = False
elif variable == 'E0Resolution':
	array = (f['cascade0Dict']['Energy'][condition]-f['cascade0TaupedeDict']['Energy'][condition])/f['cascade0Dict']['Energy'][condition]
	bins = np.linspace(-10,10,30)
	labelX = r'$\mathrm{(E_{0,True}-E_{0,Best \ fit})/E_{0,True}}$'
	isLogX = False
elif variable == 'E1Resolution':
	array = (f['cascade1Dict']['Energy'][condition]-f['cascade1TaupedeDict']['Energy'][condition])/f['cascade1Dict']['Energy'][condition]
	bins = np.linspace(-10,10,30)
	labelX = r'$\mathrm{(E_{1,True}-E_{1,Best \ fit})/E_{1,True}}$'
	isLogX = False
elif variable == 'normedChiSquared':
	array = f['taupedeDict']['chiSquared'][condition]/f['taupedeDict']['chiSquared_dof'][condition]
	bins = np.linspace(0,1000,30)
	labelX = r'$\chi^{2}/\mathrm{N_{d.o.f.}}$'
	isLogX = False
elif variable == 'PID':
	array = f['retroDict']['PID'][condition]
	bins = np.linspace(0,1,40)
	labelX = r'PID'
	isLogX = False
elif variable == 'pulseCharge':
	event = np.random.choice(f['event']['value'],1)
	height = []
	x = []
	iPulse = 0
	for charge in f['pulseDict']['charge']:
		if f['pulseDict']['event'][iPulse] == event:
			height.append(charge)
			x.append(f['pulseDict']['time'][iPulse])
		iPulse += 1
	isLogX = False
	isLogY = False
	labelX = 'Pulse time'
	labelY = 'Pulse charge'
	plot_bar(x,height,labelX,labelY,'',isLogX,isLogY,'pulse')
elif variable == 'pulseTotalCharge':
	array = f['total_charge']['value'][condition]
	#bins = np.logspace(np.log10(5),np.log10(1000),30)
	bins = np.linspace(5,200,30)
	isLogX = False
	labelX = 'Pulse total charge'
elif variable == 'domCharge':
	array = []
	weights = []
	total_charge = []
	for i in range(np.count_nonzero(condition)):
		condition_event = f['pulseDict']['event']==f['event']['value'][i]
		DOM_list = []
		string_list = []
		total_charge_tmp = []
		weights_tmp = []
		for j in range(len(f['pulseDict']['charge'][condition_event])):
			DOM = f['pulseDict']['om'][condition_event][j]
			string = f['pulseDict']['stringID'][condition_event][j]
			indexDOM = np.where(DOM_list==DOM)
			indexString = np.where(string_list==string)
			index = np.intersect1d(indexDOM,indexString)
			if index.size == 0:
				DOM_list.append(DOM)
				string_list.append(string)
				total_charge_tmp.append(f['pulseDict']['charge'][condition_event][j])
				weights_tmp.append(get_weight(inputfile,condition,U_tau4_sq,mode)[i])
			else:
				total_charge_tmp[index[0]] += f['pulseDict']['charge'][condition_event][j]
		total_charge = total_charge + total_charge_tmp
		weights = weights + weights_tmp
	array = total_charge
	weights = np.ones(len(array))
	bins = np.linspace(0,10,30)
	isLogX = False
	#isLogY = False
	labelX = 'DOM total charge'
elif variable == 'domXY':
	arrayX = f['pulseDict']['omX']
	arrayY = f['pulseDict']['omY']
	bins = np.linspace(-600,600,100)
	h, binsX, binsY, _ = plt.hist2d(arrayX,arrayY,bins)
	isLogX = False
	isLogY = False
	labelX = 'DOM X'
	labelY = 'DOM Y'
	#plot_2D_scatter(arrayX,arrayY,labelX,labelY,'',isLogX,isLogY,'domX','domY')
	plot_2D(h.T,labelX,labelY,'',r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,isLogX,isLogY,False,'domX','domY')
elif variable == 'domX':
	array = f['pulseDict']['omX']
	weights = np.ones(len(array))
	bins = np.arange(-600,600,30)
	isLogX = False
	#isLogY = False
	labelX = 'DOM X (m)'
elif variable == 'domZ':
	array = f['pulseDict']['omZ']
	weights = np.ones(len(array))
	bins = np.arange(-500,500,25)
	isLogX = False
	#isLogY = False
	labelX = 'DOM Z (m)'
elif variable == 'cosalpha-cascades':
	cosTheta0 = np.cos(f['cascade0Dict']['Zenith'][condition])
	cosPhi0 = np.cos(f['cascade0Dict']['Azimuth'][condition])
	sinTheta0 = np.sin(f['cascade0Dict']['Zenith'][condition])
	sinPhi0 = np.sin(f['cascade0Dict']['Azimuth'][condition])
	cosTheta1 = np.cos(f['cascade1Dict']['Zenith'][condition])
	cosPhi1 = np.cos(f['cascade1Dict']['Azimuth'][condition])
	sinTheta1 = np.sin(f['cascade1Dict']['Zenith'][condition])
	sinPhi1 = np.sin(f['cascade1Dict']['Azimuth'][condition])
	array = (sinTheta0*cosPhi0*sinTheta1*cosPhi1) + (sinTheta0*sinPhi0*sinTheta1*sinPhi1) + (cosTheta0*cosTheta1)
	bins = np.linspace(-1.,1.,30)
	labelX = r'$\cos(\alpha)$ between the two cascades'
	isLogX = False
elif variable == 'U_mass':
	hbar = 6.582119569E-25
	mass = f['hnlDict']['Mass'][condition]
	bins_mass = np.logspace(np.log10(0.1),np.log10(3.),30)
	bins_U = np.logspace(np.log10(1e-05),np.log10(1e-01),30)
	h = []
	for U_tau4_sq in tqdm(np.delete(bins_U,-1)):
		lifetime_proper = hbar/(FullWidth_vec(mass)*U_tau4_sq)
		weights = get_weight(inputfile,condition,U_tau4_sq,mode)
		n, bins, _ = plt.hist(mass,bins_mass,weights=weights)
		h.append(n)
	h = np.asarray(h)
	#plot_2D_contour(h,'HNL mass (GeV)',r'$|U_{\tau4}|^{2}$',labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins_mass,bins_U,True,True,True,'HNL-mass','decay-length')
	plot_2D(h,'HNL mass (GeV)',r'$|U_{\tau4}|^{2}$',labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins_mass,bins_U,True,True,False,'HNL-mass','decay-length')
elif variable == 'E0_E1':
	arrayX = f['cascade0Dict']['Energy'][condition]
	arrayY = f['cascade1Dict']['Energy'][condition]
	labelX = r'$E_{0}$ (GeV)'
	labelYbis = r'$E_{1}$ (GeV)'
	bins = np.logspace(np.log10(1e-03),np.log10(1e+04),30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,bins,weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,True,True,True,'E0','E1')
elif variable == 'E0Taupede_E1Taupede':
	arrayX = f['cascade0TaupedeDict']['Energy'][condition]
	arrayY = f['cascade1TaupedeDict']['Energy'][condition]
	labelX = r'$E_{0}$ (GeV)'
	labelYbis = r'$E_{1}$ (GeV)'
	bins = np.logspace(np.log10(1),np.log10(1e+03),30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,bins,weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,True,True,True,'E0','E1')
elif variable == 'E0True_E0BF':
	arrayX = f['cascade0Dict']['Energy'][condition]
	arrayY = f['cascade0TaupedeDict']['Energy'][condition]
	labelX = r'$\mathrm{E_{0,True}}$ (GeV)'
	labelYbis = r'$\mathrm{E_{0,Best \ fit}}$ (GeV)'
	bins = np.logspace(np.log10(1),np.log10(3e+02),30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,bins,weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,True,True,False,'E0True','E0BF')
elif variable == 'E1True_E1BF':
	arrayX = f['cascade1Dict']['Energy'][condition]
	arrayY = f['cascade1TaupedeDict']['Energy'][condition]
	labelX = r'$\mathrm{E_{1,True}}$ (GeV)'
	labelYbis = r'$\mathrm{E_{1,Best \ fit}}$ (GeV)'
	bins = np.logspace(np.log10(3),np.log10(3e+02),30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,bins,weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,True,True,False,'E1True','E1BF')
elif variable == 'LTrue_LFit':
	arrayX = f['propertiesDict']['distance'][condition]
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	labelX = r'$\mathrm{L_{True}}$ (m)'
	labelYbis = r'$\mathrm{L_{Best \ fit}}$ (m)'
	#bins = np.logspace(np.log10(1),np.log10(1e+03),30)
	bins = np.linspace(0,1000,30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,bins,weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins,bins,False,False,True,'L-True','L-bestFit')
elif variable == 'cosTheta_E':
	arrayY = np.cos(f['hnlDict']['Zenith'][condition])
	arrayX = f['hnlDict']['Energy'][condition]
	labelYbis = r'$\cos(\theta_{\mathrm{HNL}})$'
	labelX = r'$\mathrm{E_{HNL}}$'
	binsX = np.logspace(np.log10(1),np.log10(1e+03),30)
	binsY = np.linspace(-1.0,0.2,30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,True,False,False,'cosZenith-HNL','E-HNL')
elif variable == 'retro_cosTheta_E':
	arrayY = np.cos(f['retroDict']['zenith'][condition])
	arrayX = f['retroDict']['cascadeEnergy'][condition]
	labelYbis = r'$\cos(\theta_{\mathrm{retro}})$'
	labelX = r'$\mathrm{E_{retro}^{track}}$'
	binsX = np.logspace(np.log10(1),np.log10(3e+02),30)
	binsY = np.linspace(-1.0,0.2,30)
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	plot_2D(h.T,labelX,labelYbis,labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,True,False,False,'cosZenith-retro','E-retro')
elif variable == 'PID_U':
	hbar = 6.582119569E-25
	arrayX = f['retroDict']['PID'][condition]
	binsX = np.linspace(0,1,40)
	bins_U = np.logspace(np.log10(1e-05),np.log10(1e-01),30)
	h = []
	for U_tau4_sq in tqdm(np.delete(bins_U,-1)):
		weights = get_weight(inputfile,condition,U_tau4_sq,mode)
		weights[np.isnan(weights)] = 0
		n, bins, _ = plt.hist(arrayX,binsX,weights=weights)
		h.append(n)
	h = np.asarray(h)
	plot_2D(h,'PID',r'$|U_{\tau4}|^{2}$',labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,bins_U,False,True,True,'PID','U')
elif variable == 'PID_mass':
	hbar = 6.582119569E-25
	arrayX = f['retroDict']['PID'][condition]
	binsX = np.linspace(0,1,40)
	bins_mass = np.logspace(np.log10(0.1),np.log10(3),30)
	h = []
	for m in tqdm(np.delete(bins_mass,-1)):
		weights = get_weight_fixedM(inputfile,condition,1e-01,m,mode)
		weights[np.isnan(weights)] = 0
		n, bins, _ = plt.hist(arrayX,binsX,weights=weights)
		h.append(n)
	h = np.asarray(h)
	plot_2D(h,'PID',r'HNL mass (GeV)',labelY,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,bins_mass,False,True,True,'PID','mass')
elif variable == 'PID_U_mass':
	hbar = 6.582119569E-25
	bins_U = np.logspace(np.log10(1e-05),np.log10(1e-01),30)
	bins_mass = np.logspace(np.log10(0.1),np.log10(3),30)
	h = np.zeros((len(bins_mass)-1,len(bins_U)-1))
	i = 0
	for i in tqdm(xrange(len(bins_mass)-1)):
		condition_mass = (f['hnlDict']['Mass'] > bins_mass[i]) & (f['hnlDict']['Mass'] <= bins_mass[i+1])
		condition_tot = condition & condition_mass
		j = 0
		for j in xrange(len(bins_U)-1):
			weights = get_weight_fixedM(inputfile,condition_tot,bins_U[j],bins_mass[i],mode)
			h[i][j] = np.sum(weights*f['retroDict']['PID'][condition_tot])/np.sum(weights)
			j += 1
		i += 1
	h = np.asarray(h)
	h[np.isnan(h)] = 0
	plot_2D(h.T,'HNL mass (GeV)',r'$|U_{\tau4}|^{2}$','Mean PID',r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',bins_mass,bins_U,True,True,False,'mass','U')


if variable != 'U_mass':
	yerr = histogram_errors(array,weights,bins)
	plot_hist_1D(array,labelX,labelY,'',bins,weights,yerr,True,isLogX,isLogY,variable)
