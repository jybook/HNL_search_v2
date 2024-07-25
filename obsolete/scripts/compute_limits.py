import h5py
import matplotlib.pyplot as plt
import numpy as np
from modules.plot import *
from modules.utils import *
from modules.DecayWidths import FullWidth
import math
from tqdm import tqdm
import sys, argparse
from skhep.dataset.numpydataset import *
import uproot
from skhep.dataset.selection import Selection
import ROOT
#from Utilities.utilities import destruct_objects
#from Utilities.RooFit import RooDataset, RemoveEmptyBins
#from PyLHCb.Root.RooFitUtils import ResidualPlot
#import probfit
import iminuit
from root_numpy import array2tree, array2hist
import pandas as pd

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

f = [h5py.File(inputfile[0], 'r'),h5py.File(inputfile[1], 'r'),h5py.File(inputfile[2], 'r'),h5py.File(inputfile[3], 'r'),h5py.File(inputfile[4], 'r')]
nFiles = [float(inputfile[0].split('.')[0].split('_')[-1]),float(inputfile[1].split('.')[0].split('_')[-1]),float(inputfile[2].split('.')[0].split('_')[-1]),float(inputfile[3].split('.')[0].split('_')[-1]),float(inputfile[4].split('.')[0].split('_')[-1])]
level = inputfile[0].split('.')[0].split('_')[-2]

# Conditions
# "No condition" (always True)
cdt0 = [np.ones(len(f[0]['neutrinoDict']['Energy']), dtype=bool),
        np.ones(len(f[1]['neutrinoDict']['Energy']), dtype=bool),
        np.ones(len(f[2]['neutrinoDict']['Energy']), dtype=bool),
        np.ones(len(f[3]['neutrinoDict']['Energy']), dtype=bool),
        np.ones(len(f[4]['leptonDict']['Energy']), dtype=bool)]
# Final condition
cdt = cdt0#[cdt_M,cdt0[1],cdt0[2],cdt0[3],cdt0[4]]

# Weighting
U_tau4_sq = 1e-01
weights = [get_weight(inputfile[0],cdt[0],U_tau4_sq,''),1000*f[1]['weight']['value'][cdt[1]]/nFiles[1],1000*f[2]['weight']['value'][cdt[2]]/nFiles[2],1000*f[3]['weight']['value'][cdt[3]]/nFiles[3],1000*f[4]['weight']['value'][cdt[4]]/nFiles[4]]

# Determine data (observed), bg, and signal
signal = f[0]['hnlDict']['Mass']
bg = np.sum(weights[1])+sum(weights[2])+sum(weights[3])+sum(weights[4])

# Constant mass model parameter, we will manuallly set it to scan over it.
mass = ROOT.RooRealVar("mass", "mass", 1.0, 0.1, 3.0)
mass.setConstant(True)

# Dummy observable for the RooUniforms.
x = ROOT.RooRealVar("x", "x", 0.0)

# The number of points to scan over for our mass model parameter, and the actual values.
n_points = 30
mass_vals = np.logspace(np.log10(0.1), np.log10(3.0), n_points+1)

# Signal model (histogram PDF) 
hist_sig_temp = ROOT.TH1D("sigHist","sigHist",n_points,mass_vals)
array_binned_sig, _ = np.histogram(signal,mass_vals,weights=weights[0])
hist_sig = array2hist(array_binned_sig,hist_sig_temp)
sigDataHist = ROOT.RooDataHist("sigDataHist","sigDataHist",ROOT.RooArgSet(mass),hist_sig)
nsig_bsm = ROOT.RooHistPdf("sigHistPdf","sigHistPdf",mass,sigDataHist,0)

# The parameter of interest, aka signal strength parameter.
mu = ROOT.RooRealVar("mu", "mu", 1, 0.0, 100000000.0)

# The number of signal is the standard model prediction times the signal strength.
nsig = ROOT.RooProduct("nsig", "nsig", ROOT.RooArgList(mu, nsig_bsm))
# The number of background has a gaussian uncertainty
nbkg = ROOT.RooRealVar("nbkg", "nbkg", bg, 0.0, 1.1*bg)

# Signal and background pdfs (uniforms PDFs if no systematics)
sig = ROOT.RooUniform("sig", "sig", x)
bkg = ROOT.RooUniform("bkg", "bkg", x)
#mean_bg = ROOT.RooFit.RooConst(bg)
#sigma_bg = ROOT.RooFit.RooConst(0.1*bg)
#bkg = ROOT.RooGaussian("bkg", "bkg", x, mean_bg, sigma_bg)

# The final model for the counting experiment.
model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(sig, bkg), ROOT.RooArgList(nsig, nbkg))

# A dictionary to store the expected upper limits plus +/- 1 and 2 sigma confidence bands.
result_dict = {}
for i in [-2, -1, 0, 1, 2]:
    result_dict[i] = np.zeros(n_points + 1)

# Now, do the actual scan over model parameters.
for i_point in range(n_points + 1):
	print(i_point)
	mass.setVal(mass_vals[i_point])
	
	# Create a workspace for the model for given mass value.
	ws = ROOT.RooWorkspace()
	ws.Import(model)
	
	conf = ROOT.RooStats.ModelConfig("model" + str(i_point), ws)
	conf.SetPdf(model)
	conf.SetParametersOfInterest(mu)
	conf.SetObservables(x)
	conf.SetNuisanceParameters(nbkg)
	
	# S+B model
	model_sb = conf
	model_sb.SetName("MODEL_SB" + str(i_point))
	
	# BKG only
	model_b = conf.Clone()
	model_b.SetName("MODEL_B" + str(i_point))
	
	mu.setVal(1.0)
	model_sb.SetSnapshot(ROOT.RooArgSet(mu))
	
	mu.setVal(0.0)
	model_b.SetSnapshot(ROOT.RooArgSet(mu))
	
	# Generate simulated observed data according to null hypothesis.
	mu.setVal(0.0)
	data = model.generate(x)
	
	# Test statistic calculator.
	calc = ROOT.RooStats.AsymptoticCalculator(data, model_b, model_sb)
	calc.SetOneSided(True)
	calc.SetQTilde(False)
	
	# Create, configure and run hypothesis test inverted.
	test = ROOT.RooStats.HypoTestInverter(calc)
	test.SetConfidenceLevel(0.9)
	test.UseCLs(True)
	test.SetFixedScan(5000, 0, 200000.0)
	r = test.GetInterval()

	# Store the results in our dictionary.
	for i in [-2, -1, 0, 1, 2]:
		result_dict[i][i_point] = r.GetExpectedUpperLimit(i)

# Plotting

plt.clf()
fig, ax = plt.subplots()
plt.fill_between(
    mass_vals, result_dict[-2], result_dict[2], color="#ffff00", step="pre", label="expected limit (+/- 2 sig)"
)
plt.fill_between(
    mass_vals, result_dict[-1], result_dict[1], color="#00ff00", step="pre", label="expected limit (+/- 1 sig)"
)
plt.step(mass_vals, result_dict[0], "k--", label="expected limit (median)")
ax.text(0.15,0.7,r'$|U_{\tau4}|^{2} = 10^{-1}$',transform=ax.transAxes)
plt.legend(loc='best')
plt.xlim(0.1, 3.)
plt.ylim(1e+01, 1e+08)
plt.xlabel('HNL mass (GeV)')
plt.xscale('log')
plt.ylabel('Signal strength upper limit (90% CLs)')
plt.yscale('log')
fig.tight_layout()
plt.savefig('brazil.png')
plt.show()
