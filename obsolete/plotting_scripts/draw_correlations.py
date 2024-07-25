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
	#print(pdf_inverse,pdf_exp1,pdf_exp2,pdf_exp)
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

def getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights):
	result = []
	mean = np.zeros(len(binsX)-1)
	total = np.zeros(len(binsX)-1)
	correlation = 0
	sumTot = 0
	sumx = 0
	sumx2 = 0
	sumy = 0
	sumy2 = 0
	sumxy = 0
	for i in range(len(arrayY)):
		if (arrayY[i] < binsY[0]) or (arrayY[i] > binsY[len(binsY)-1]) or (arrayX[i] < binsX[0]) or (arrayX[i] > binsX[len(binsX)-1]): continue
		sumTot += weights[i]
		sumx  += arrayX[i]*weights[i]
		sumx2  += arrayX[i]*arrayX[i]*weights[i]
		sumy  += arrayY[i]*weights[i]
		sumy2  += arrayY[i]*arrayY[i]*weights[i]
		sumxy += arrayX[i]*arrayY[i]*weights[i]
		for j in range(len(binsX)-1):
			if (arrayX[i] >= binsX[j]) and (arrayX[i] < binsX[j+1]):
				mean[j] += arrayY[i]*weights[i]
				total[j] += weights[i]
	mean = mean/total
	sumxNorm = sumx/sumTot
	sumyNorm = sumy/sumTot
	rmsX = np.sqrt(np.abs((sumx2/sumTot)-(sumxNorm*sumxNorm)))
	rmsY = np.sqrt(np.abs((sumy2/sumTot)-(sumyNorm*sumyNorm)))
	correlation = ((sumxy/sumTot) - (sumx/sumTot)*(sumy/sumTot))/(rmsX*rmsY)
	result.append(mean)
	result.append(correlation)
	return result

def getMeanAndCorrelation2(arrayY,binsX,binsY):
	result = []
	mean = np.zeros(len(binsX)-1)
	total = np.zeros(len(binsX)-1)
	correlation = 0
	sumTot = 0
	sumx = 0
	sumx2 = 0
	sumy = 0
	sumy2 = 0
	sumxy = 0
	for j in tqdm(range(len(binsX)-1)):
		weights = get_weight(inputfile,condition,binsX[j],mode)
		weights[np.isnan(weights)] = 0
		weights[weights<1e-20] = 0
		for i in range(len(arrayY)):
			if (arrayY[i] < binsY[0]) or (arrayY[i] > binsY[len(binsY)-1]): continue
			sumTot += weights[i]
			sumx  += binsX[j]*weights[i]
			sumx2  += binsX[j]*binsX[j]*weights[i]
			sumy  += arrayY[i]*weights[i]
			sumy2  += arrayY[i]*arrayY[i]*weights[i]
			sumxy += binsX[j]*arrayY[i]*weights[i]
			mean[j] += arrayY[i]*weights[i]
	mean = mean/sumTot
	sumxNorm = sumx/sumTot
	sumyNorm = sumy/sumTot
	rmsX = np.sqrt(np.abs((sumx2/sumTot)-(sumxNorm*sumxNorm)))
	rmsY = np.sqrt(np.abs((sumy2/sumTot)-(sumyNorm*sumyNorm)))
	correlation = ((sumxy/sumTot) - (sumx/sumTot)*(sumy/sumTot))/(rmsX*rmsY)
	result.append(mean)
	result.append(correlation)
	return result

# Vectorization
weight_lifetime_vec = np.vectorize(weight_lifetime)
FullWidth_vec = np.vectorize(lambda x: FullWidth(x))

# Define input arguments
parser = argparse.ArgumentParser(description='Plot the selected quantity for the selected particle')
parser.add_argument('-i', '--inputfile', type=str, default='190606_L7_995.h5', help='Name of the file')
parser.add_argument('-v', '--variable', type=str, default='energy', help='Name of the variable (see below)')
parser.add_argument('-m', '--mode', type=str, default='rate', help='rate or event')
args = parser.parse_args()

inputfile = args.inputfile
variable = args.variable
mode = args.mode

f = h5py.File(inputfile, 'r')
nFiles = float(inputfile.split('.')[0].split('_')[-1])
level = inputfile.split('.')[0].split('_')[-2]

labelZ = 'Rate (mHz)'
if mode == 'event':
	labelZ = 'Events (8 years)'
isLogX = True
isLogY = True

# Kinematic cut: uncomment the desired condition
if '1906' in inputfile:
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
	cdt_DC = (distance > 50) & ((inDeepCore0 & (f['cascade0Dict']['Energy'] > 5) & inDeepCore1 & (f['cascade1Dict']['Energy'] > 5)) | (inDeepCore0 & (f['cascade0Dict']['Energy'] > 5) & inIceCube1 & (f['cascade1Dict']['Energy'] > 20)) | (inDeepCore1 & (f['cascade1Dict']['Energy'] > 5) & inIceCube0 & (f['cascade0Dict']['Energy'] > 20)) | (inIceCube0 & (f['cascade0Dict']['Energy'] > 20) & inIceCube1 & (f['cascade1Dict']['Energy'] > 20)))
	# Strong "reconstructible" double cascade condition
	cdt_strongDC = (distance > 30) & ((inDeepCore0 & (f['cascade0Dict']['Energy'] > 10) & inDeepCore1 & (f['cascade1Dict']['Energy'] > 10)) | (inDeepCore0 & (f['cascade0Dict']['Energy'] > 10) & inIceCube1 & (f['cascade1Dict']['Energy'] > 50)) | (inDeepCore1 & (f['cascade1Dict']['Energy'] > 10) & inIceCube0 & (f['cascade0Dict']['Energy'] > 50)) | (inIceCube0 & (f['cascade0Dict']['Energy'] > 50) & inIceCube1 & (f['cascade1Dict']['Energy'] > 50)))
	# Single cascade condition
	cdt_SC = (distance < 20) | ((distance > 20) & ((radius0 > 600) | (radius1 > 600) | (f['cascade0Dict']['Z'] < -500) | (f['cascade1Dict']['Z'] < -500) | (f['cascade0Dict']['Z'] > 500) | (f['cascade1Dict']['Z'] > 500)))
	# One cascade in I3 condition
	cdt_1C = (radius0 > 600) | (radius1 > 600) | (f['cascade0Dict']['Z'] < -500) | (f['cascade1Dict']['Z'] < -500) | (f['cascade0Dict']['Z'] > 500) | (f['cascade1Dict']['Z'] > 500)
	# Both cascades in DeepCore condition
	cdt_DeepCore = inDeepCore0 & inDeepCore1
	# Condition on the mass range
	cdt_m = (f['hnlDict']['Mass'] > 0.5)# & (f['hnlDict']['Mass'] < 1.8)
# Track condition
#cdt_T = f['retroDict']['PID'] > 0.8
# Basic condition on energies
cdt_E = (f['cascade0TaupedeDict']['Energy'] > 5) & (f['cascade1TaupedeDict']['Energy'] > 5)
# This is "no condition" (always True)
cdt0 = np.ones(len(f['leptonDict']['Energy']), dtype=bool)
# Final condition
condition = cdt0

#print(str(len(f['neutrinoDict']['Energy'][cdt_DC]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_DC]))/len(f['neutrinoDict']['Energy']))+'%) pass the DC condition')
#print(str(len(f['neutrinoDict']['Energy'][cdt_SC]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_SC]))/len(f['neutrinoDict']['Energy']))+'%) pass the SC condition')
#print(str(len(f['neutrinoDict']['Energy'][cdt_1C]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_1C]))/len(f['neutrinoDict']['Energy']))+'%) pass the 1C condition')
#print(str(len(f['neutrinoDict']['Energy'][cdt_DC & cdt_T]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_DC & cdt_T]))/len(f['neutrinoDict']['Energy']))+'%) pass the DC and the T conditions')
#print(str(len(f['neutrinoDict']['Energy'][cdt_SC & cdt_T]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][cdt_SC & cdt_T]))/len(f['neutrinoDict']['Energy']))+'%) pass the SC and the T conditions')
#print(str(len(f['neutrinoDict']['Energy'][condition]))+' events '+'('+str(100*float(len(f['neutrinoDict']['Energy'][condition]))/len(f['neutrinoDict']['Energy']))+'%) pass the final condition')

# Weighting (choose a value for the mixing parameter - does not apply for U-mass 2D plots)
U_tau4_sq = 1e-03
#weights = get_weight(inputfile,condition,U_tau4_sq,mode)
if not '1906' in inputfile:
	nFiles = float(inputfile.split('.')[0].split('_')[-1])
	weights = 1000*f['weight']['value']/nFiles
#for i in range(len(mass)):
#	if mass[i]>2: print(f['weight']['value'][condition][i],weights[i])
weights[np.isnan(weights)] = 0
#weights[weights<1e-20] = 0
#weights = np.ones(len(weights))

array = []
bins = []

#nNaN = 0
#for e in f['cascade1Dict']['Energy']:
#	if np.isnan(e): nNaN += 1
#print(str(100*float(nNaN)/len(f['cascade1Dict']['Energy']))+' % of NaN')

if variable == 'cosTheta_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	#binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	binsX = np.arange(0.1,3,0.1)
	isLogX = False
	labelX = 'HNL mass (GeV)'
	arrayY = np.cos(f['cascade0TaupedeDict']['Zenith'][condition])
	binsY = np.linspace(-1,1,30)
	labelY = r'Taupede best fit $\cos(\theta)$'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','cosTheta')
elif variable == 'cosTheta0True_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	#binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	binsX = np.arange(0.1,3,0.1)
	isLogX = False
	labelX = 'HNL mass (GeV)'
	arrayY = np.cos(f['cascade0Dict']['Zenith'][condition])
	binsY = np.linspace(-1,1,30)
	labelY = r'First cascade true $\cos(\theta)$'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','cosTheta0True')
elif variable == 'cosThetaTaupede_cosThetaTrue':
	#arrayX = np.cos(f['cascade0Dict']['Zenith'][condition])
	arrayX = np.cos(f['leptonDict']['Zenith'][condition])
	binsX = np.linspace(-1,1,30)
	isLogX = False
	labelX = r'$\cos(\theta_{\mathrm{True}})$'
	arrayY = np.cos(f['cascade0TaupedeDict']['Zenith'][condition])
	#arrayY = np.cos(f['retroDict']['zenith'][condition])
	binsY = np.linspace(-1,1,30)
	labelY = r'$\cos(\theta_{\mathrm{Taupede}})$'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'cosTheta0True','cosThetaTaupede')
elif variable == 'cosThetaRetro_cosThetaTrue':
	arrayX = np.cos(f['cascade1Dict']['Zenith'][condition])
	binsX = np.linspace(-1,1,30)
	isLogX = False
	labelX = r'$\cos(\theta_{\mathrm{True}})$'
	arrayY = np.cos(f['retroDict']['zenith'][condition])
	binsY = np.linspace(-1,1,30)
	labelY = r'$\cos(\theta_{\mathrm{retro}})$'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'cosTheta0True','cosThetaRetro')
elif variable == 'cosTheta_U':
	binsX = np.logspace(np.log10(1e-05),np.log10(1),30)
	labelX = r'$|U_{\tau4}|^{2}$'
	arrayY = np.cos(f['cascade0TaupedeDict']['Zenith'][condition])
	binsY = np.linspace(-1,1,30)
	labelY = r'Taupede best fit $\cos(\theta)$'
	isLogY = False
	h = []
	for U_tau4_sq in tqdm(np.delete(binsX,-1)):
		weights = get_weight(inputfile,condition,U_tau4_sq,mode)
		weights[np.isnan(weights)] = 0
		weights[weights<1e-20] = 0
		n, bins, _ = plt.hist(arrayY,binsY,weights=weights)
		h.append(n)
	h = np.asarray(h)
	meanCorr = getMeanAndCorrelation2(arrayY,binsX,binsY)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'U','cosTheta')
elif variable == 'E0_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	labelX = 'HNL mass (GeV)'
	arrayY = f['cascade0TaupedeDict']['Energy'][condition]
	binsY = np.logspace(np.log10(5),np.log10(1000),30)
	labelY = r'Taupede best fit $E_{0}$'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','E0')
elif variable == 'E0_U':
	binsX = np.logspace(np.log10(1e-05),np.log10(1),30)
	labelX = r'$|U_{\tau4}|^{2}$'
	arrayY = f['cascade0TaupedeDict']['Energy'][condition]
	binsY = np.logspace(np.log10(5),np.log10(1000),30)
	labelY = r'Taupede best fit $E_{0}$'
	h = []
	for U_tau4_sq in tqdm(np.delete(binsX,-1)):
		weights = get_weight(inputfile,condition,U_tau4_sq,mode)
		weights[np.isnan(weights)] = 0
		weights[weights<1e-20] = 0
		n, bins, _ = plt.hist(arrayY,binsY,weights=weights)
		h.append(n)
	h = np.asarray(h)
	meanCorr = getMeanAndCorrelation2(arrayY,binsX,binsY)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'U','E0')
elif variable == 'E1_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	labelX = 'HNL mass (GeV)'
	arrayY = f['cascade1TaupedeDict']['Energy'][condition]
	binsY = np.logspace(np.log10(5),np.log10(1000),30)
	labelY = r'Taupede best fit $E_{1}$'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,True,True,True,'mass','E1')
elif variable == 'L_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	#binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	binsX = np.arange(0.1,3,0.1)
	isLogX = False
	labelX = 'HNL mass (GeV)'
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	#binsY = np.logspace(np.log10(1),np.log10(1000),30)
	binsY = np.arange(0,800,20)
	isLogY = False
	labelY = r'Best fit distance (m)'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','BF-distance')
elif variable == 'L_ETot':
	arrayX = f['cascade0TaupedeDict']['Energy'][condition]+f['cascade1TaupedeDict']['Energy'][condition]
	binsX = np.logspace(np.log10(5),np.log10(5e+03),30)
	labelX = 'Total energy (Taupede best-fit) (GeV)'
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	#binsY = np.logspace(np.log10(1),np.log10(1000),30)
	binsY = np.arange(0,500,30)
	isLogY = False
	labelY = r'Cascade separation (Taupede best-fit) (m)'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'ETot','BF-distance')
elif variable == 'L_cosTheta':
	arrayX = np.cos(f['cascade0TaupedeDict']['Zenith'][condition])
	binsX = np.linspace(-1,1,30)
	isLogX = False
	labelX = r'$\cos(\theta)$ (Taupede best-fit)'
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	#binsY = np.logspace(np.log10(1),np.log10(1000),30)
	binsY = np.arange(0,500,30)
	isLogY = False
	labelY = r'Cascade separation (Taupede best-fit) (m)'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'cosTheta','BF-distance')
if variable == 'L_U':
	binsX = np.logspace(np.log10(1e-05),np.log10(1),30)
	labelX = r'$|U_{\tau4}|^{2}$'
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	binsY = np.arange(0,400,10)
	isLogY = False
	labelY = r'Best fit distance (m)'
	h = []
	for U_tau4_sq in tqdm(np.delete(binsX,-1)):
		weights = get_weight(inputfile,condition,U_tau4_sq,mode)
		weights[np.isnan(weights)] = 0
		weights[weights<1e-20] = 0
		n, bins, _ = plt.hist(arrayY,binsY,weights=weights)
		h.append(n)
	h = np.asarray(h)
	meanCorr = getMeanAndCorrelation2(arrayY,binsX,binsY)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'U','L')
elif variable == 'LTrue_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	#binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	binsX = np.arange(0.1,3,0.1)
	isLogX = False
	labelX = 'HNL mass (GeV)'
	arrayY = f['propertiesDict']['distance'][condition]
	#binsY = np.logspace(np.log10(1),np.log10(1000),30)
	binsY = np.arange(0,800,20)
	isLogY = False
	labelY = r'True distance (m)'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','distance')
elif variable == 'L_LTrue':
	arrayX = f['propertiesDict']['distance'][condition]
	binsX = np.arange(50,400,10)
	isLogX = False
	labelX = 'True distance (m)'
	arrayY = f['taupedeDict']['bestFitLength'][condition]
	binsY = np.arange(50,400,10)
	isLogY = False
	labelY = r'Best fit distance (m)'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	meanCorr = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'distance','BF-distance')
elif variable == 'LRes_mass':
	arrayX = f['hnlDict']['Mass'][condition]
	#binsX = np.logspace(np.log10(0.1),np.log10(3),30)
	binsX = np.arange(0.1,3,0.1)
	isLogX = False
	labelX = 'HNL mass (GeV)'
	arrayY = (f['propertiesDict']['distance'][condition]-f['taupedeDict']['bestFitLength'][condition])/f['propertiesDict']['distance'][condition]
	binsY = np.linspace(-10,10,25)
	isLogY = False
	labelY = r'Separation resolution'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	meanCorr = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','separation-res')
elif variable == 'LRes_LTrue':
	arrayX = f['propertiesDict']['distance'][condition]
	binsX = np.arange(0,400,10)
	isLogX = False
	labelX = 'True distance (m)'
	arrayY = (f['propertiesDict']['distance'][condition]-f['taupedeDict']['bestFitLength'][condition])/f['propertiesDict']['distance'][condition]
	binsY = np.linspace(-10,10,25)
	isLogY = False
	labelY = r'Separation resolution'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	meanCorr = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'mass','separation-res')
elif variable == 'LRes_NuETrue':
	arrayX = f['neutrinoDict']['Energy'][condition]
	binsX = np.logspace(np.log10(1),np.log10(1000),30)
	labelX = 'True neutrino energy (GeV)'
	arrayY = (f['propertiesDict']['distance'][condition]-f['taupedeDict']['bestFitLength'][condition])/f['propertiesDict']['distance'][condition]
	binsY = np.linspace(-10,10,25)
	isLogY = False
	labelY = r'Separation resolution'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	meanCorr = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'neurtino-energy','separation-res')
elif variable == 'LRes_ETot':
	arrayX = f['cascade0Dict']['Energy'][condition]+f['cascade1Dict']['Energy'][condition]
	binsX = np.logspace(np.log10(1),np.log10(10000),30)
	labelX = 'Total energy (GeV)'
	arrayY = (f['propertiesDict']['distance'][condition]-f['taupedeDict']['bestFitLength'][condition])/f['propertiesDict']['distance'][condition]
	binsY = np.linspace(-10,10,25)
	isLogY = False
	labelY = r'Separation resolution'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	meanCorr = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)
	print('Correlation = '+str(meanCorr[1]))
	plot_2D(h.T,meanCorr[0],labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'ETot','separation-res')
elif variable == 'E0_E0True':
	arrayX = f['cascade0Dict']['Energy'][condition]
	binsX = np.logspace(np.log10(5),np.log10(1000),30)
	labelX = r'True $E_{0}$'
	arrayY = f['cascade0TaupedeDict']['Energy'][condition]
	binsY = np.logspace(np.log10(5),np.log10(1000),30)
	labelY = r'Taupede best fit $E_{0}$'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,True,True,True,'E0True','E0')
elif variable == 'nuE_conf':
	arrayX = f['neutrinoDict']['Energy'][condition]
	binsX = np.logspace(np.log10(1),np.log10(1000),30)
	labelX = r'Neutrino true energy (GeV)$'
	arrayY = (f['cascade0TaupedeDict']['Energy'][condition]+f['cascade1TaupedeDict']['Energy'][condition])/f['millipedeDict']['E_tot'][condition]
	binsY = np.linspace(0,2,30)
	isLogY = False
	labelY = r'Energy confinement'
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'nuE','conf')
elif variable == 'x0True_x0Reco':
	arrayX = f['cascade0Dict']['X'][condition]
	binsX = np.linspace(-600,600,31)
	labelX = r'First cascade true X (m)'
	isLogX = False
	arrayY = f['cascade0TaupedeDict']['X'][condition]
	binsY = np.linspace(-600,600,31)
	labelY = r'First cascade reconstructed X (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'x0True','x0Reco')
elif variable == 'x1True_x1Reco':
	arrayX = f['cascade1Dict']['X'][condition]
	binsX = np.linspace(-600,600,31)
	labelX = r'Second cascade true X (m)'
	isLogX = False
	arrayY = f['cascade1TaupedeDict']['X'][condition]
	binsY = np.linspace(-600,600,31)
	labelY = r'Second cascade reconstructed X (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'x1True','x1Reco')
elif variable == 'y0True_y0Reco':
	arrayX = f['cascade0Dict']['Y'][condition]
	binsX = np.linspace(-600,600,31)
	labelX = r'First cascade true Y (m)'
	isLogX = False
	arrayY = f['cascade0TaupedeDict']['Y'][condition]
	binsY = np.linspace(-600,600,31)
	labelY = r'First cascade reconstructed Y (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'y0True','y0Reco')
elif variable == 'y1True_y1Reco':
	arrayX = f['cascade1Dict']['Y'][condition]
	binsX = np.linspace(-600,600,31)
	labelX = r'Second cascade true Y (m)'
	isLogX = False
	arrayY = f['cascade1TaupedeDict']['Y'][condition]
	binsY = np.linspace(-600,600,31)
	labelY = r'Second cascade reconstructed Y (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'y1True','y1Reco')
elif variable == 'z0True_z0Reco':
	arrayX = f['cascade0Dict']['Z'][condition]
	binsX = np.linspace(-500,500,31)
	labelX = r'First cascade true Z (m)'
	isLogX = False
	arrayY = f['cascade0TaupedeDict']['Z'][condition]
	binsY = np.linspace(-500,500,31)
	labelY = r'First cascade reconstructed Z (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'z0True','z0Reco')
elif variable == 'z1True_z1Reco':
	arrayX = f['cascade1Dict']['Z'][condition]
	binsX = np.linspace(-500,500,31)
	labelX = r'Second cascade true Z (m)'
	isLogX = False
	arrayY = f['cascade1TaupedeDict']['Z'][condition]
	binsY = np.linspace(-500,500,31)
	labelY = r'Second cascade reconstructed Z (m)'
	isLogY = False
	h, xedges, yedges, _ = plt.hist2d(arrayX,arrayY,[binsX,binsY],weights=weights)
	mean = getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[0]
	print('Correlation = '+str(getMeanAndCorrelation(arrayX,arrayY,binsX,binsY,weights)[1]))
	plot_2D(h.T,mean,labelX,labelY,labelZ,r'$\nu_{\tau}$ NC up-scattering events (LI, '+level+')',binsX,binsY,isLogX,isLogY,True,'z1True','z1Reco')
