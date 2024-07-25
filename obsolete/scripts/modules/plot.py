import h5py
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.lines as lines
#import matplotlib.image.pcolormesh
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
#from utils import MidpointNormalize

# Function to plot one 1D histograms
def plot_hist_1D(data, label_X, label_Y, title, bins, weights, yerr, drawErr, isXLog, isYLog, variable):
	# Define canvas geometry
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	# Compute bins and x-axis error bars
	xerr = np.zeros(len(yerr))
	for i in range(len(xerr)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	# Define the histogram
	n, bins, _ = plt.hist(data,bins,weights=weights,histtype='step',label='',color='white')
	# Print total rates
	print('Total rate = '+str(np.sum(n))+' mHz')
	# Compute the "per bin" histogram
	#color = 'blue'
	#if drawErr: color = 'white'
	n_perBin = np.divide(n,2*xerr)
	yerr_perBin = np.divide(yerr,2*xerr)
	if drawErr: plt.errorbar(bin_centers,n_perBin,xerr=xerr,yerr=yerr_perBin,fmt='.b')
	#if drawErr: plt.errorbar(bin_centers,n_perBin,xerr=xerr,yerr=0,fmt='.b')
	# Compute mean
	mean = 0
	total = 0
	for i in range(len(data)):
		if(data[i] < bins[0] or data[i] > bins[-1]): continue
		mean += data[i]*weights[i]
		total += weights[i]
	mean = mean/total
	print('Mean = '+str(mean))
	# X axis 
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	#plt.ylim(1e-13,1E-05)
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	# Other settings
	plt.tick_params(labelsize=18)
	plt.title(title)
	plt.title(title,fontsize=16,pad=20)
	fig.tight_layout()
	#plt.legend(loc='best',frameon=False)
	#ax.text(0.6,0.9,'GENIE L7',transform=ax.transAxes,fontsize=15)
	#ax.text(0.7,0.9,'HNL (post-L7)',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.85,'7.5 years',transform=ax.transAxes,fontsize=15)
	# Save and display figure
	plt.savefig(variable+'.pdf')
	plt.savefig(variable+'.png')
	plt.show()

# Function to plot one 1D graph
def plot_graph_1D(n, err, bins, label_X, label_Y, title, isXLog, isYLog, particle, quantity_X, quantity_Y, isRatio):
	plt.clf()
	fig, ax = plt.subplots()
	#x = np.zeros(len(n))
	xerr = np.zeros(len(n))
	for i in range(len(n)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	plt.errorbar(bin_centers,n,xerr=xerr,yerr=err,fmt='o',markersize=4)
	plt.xlabel(label_X)
	plt.title(title)
	if isYLog: plt.yscale('log')
	if isXLog: plt.xscale('log')
	if isRatio:
		# Plot line at y mean
		yMean = 0
		total = 0
		for i in range(len(n)):
			yMean += pow(n[i]/err[i],2)*n[i]
			total += pow(n[i]/err[i],2)
		yMean = yMean/total
		#line = plt.axhline(yMean,color='red')
		line = plt.axhline(1,color='red')
		# Zoom on +/- 20% around the line
		#plt.ylim(0.9*yMean,1.1*yMean)
		#plt.ylim(0.75,1.0)
		#plt.ylim(0.9,1.4)
		#ax.legend([line],['Mean ratio'],loc='upper left',frameon=False)
		#ax.text(0.5,0.9,'Nominal (David) / Nominal (Spencer)',transform=ax.transAxes)
		#ax.text(0.7,0.92,'+200 m / Nominal',transform=ax.transAxes)
		ax.text(0.15,0.92,'+200 m / Nominal',transform=ax.transAxes)
		if quantity_Y == '':
			plt.ylabel('Ratio of rates')
			plt.savefig(particle+'-'+quantity_X+'_ratio-of-rates.pdf')
		else:
			plt.ylabel('Ratio of mean energies')
			plt.savefig('mean-'+particle+'-'+quantity_Y+'_as_'+quantity_X+'_ratio.pdf')
	else:
		plt.ylabel('Mean '+label_Y)
		plt.savefig('mean-'+particle+'-'+quantity_Y+'_as_'+quantity_X+'.pdf')
	plt.show()

# Function to plot one 1D graph
def plot_bar(x, height, label_X, label_Y, title, isXLog, isYLog, quantity):
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)

	plt.bar(x,height)

	# X axis 
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	# Other settings
	plt.tick_params(labelsize=18)
	plt.title(title)
	plt.title(title,fontsize=16,pad=20)
	fig.tight_layout()

	plt.savefig('bar_'+quantity+'.png')
	plt.savefig('bar_'+quantity+'.pdf')
	plt.show()

# Function to plot two 1D histograms
def plot_2hist_1D(data1, data2, label_X, label_Y, title, bins, weights1, weights2, yerr1, yerr2, drawErr, isXLog, isYLog, variable):
	# Define canvas geometry
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	# Compute bins and x-axis error bars
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	color1 = 'black'
	color2 = 'red'
	if drawErr:
		color1 = 'white'
		color2 = 'white'
	n1, bins, _ = plt.hist(data1,bins,weights=weights1,histtype='step',density=False,color=color1,fill=False)
	n2, bins, _ = plt.hist(data2,bins,weights=weights2,histtype='step',density=False,color=color2,fill=False)
	n1_perBin = np.divide(n1,2*xerr)
	n2_perBin = np.divide(n2,2*xerr)
	# Compute total rates
	integral1 = 0
	integral2 = 0
	for i in range(len(n1)):
		integral1 += n1[i]
		integral2 += n2[i]
	print('Total rate 1 = '+str(integral1))
	print('Total rate 2 = '+str(integral2))

	if drawErr:
		yerr1_perBin = np.divide(yerr1,2*xerr)
		yerr2_perBin = np.divide(yerr2,2*xerr)
		errorbar1 = ax.errorbar(bin_centers,n1,xerr=xerr,yerr=yerr1,fmt='bo',markersize=4)
		errorbar2 = ax.errorbar(bin_centers,n2,xerr=xerr,yerr=yerr2,fmt='ro',markersize=4)

	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	#plt.ylim(1E-06,1E-01)
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(which='both',labelsize=18)
	fig.tight_layout()

	#if drawErr: ax.legend([errorbar2,errorbar1],['NuMu events','HNL events'],loc=loc,frameon=False,fontsize=15)
	##else: plt.legend([n2,n1],['NuMu events','HNL events'],loc=loc,frameon=False,fontsize=15)
	#else: plt.legend(loc='best',frameon=False,fontsize=15)
	if drawErr: ax.legend([errorbar1,errorbar2],['Gen-level (rate = '+str(round(1000*integral1,1))+' mHz)','Det-level (rate = '+str(round(1000*integral2,1))+' mHz)'],loc='best',frameon=False,fontsize=15)
	else: ax.legend(loc='best',frameon=False,fontsize=15)

	plt.savefig(variable+'_2hist.png')
	plt.savefig(variable+'_2hist.pdf')

	plt.show()

# Function to plot two 1D histograms and their ratio
def plot_2hist_ratio(data1, data2, label_X, label_Y, title, bins, weights1, weights2, yerr1, yerr2, drawErr, isXLog, isYLog, variable):
	### Define canvas geometry
	plt.clf()
	fig = plt.figure()
	#plt.subplots(1,2,sharex=True)
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	#spec = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[3,1])
	spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3,1])
	### Upper plot
	ax = fig.add_subplot(spec[0,0])
	### Compute bins and x-axis error bars
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	color1 = 'black'
	color2 = 'red'
	if drawErr:
		color1 = 'white'
		color2 = 'white'
	density = False
	n1, bins, _ = plt.hist(data1,bins,weights=weights1,histtype='step',color=color1,fill=False)
	n2, bins, _ = plt.hist(data2,bins,weights=weights2,histtype='step',color=color2,fill=False)
	n1_perBin = np.divide(n1,2*xerr)
	n2_perBin = np.divide(n2,2*xerr)
	### Compute total rates
	integral1 = np.sum(n1)
	integral2 = np.sum(n2)
	print('Total rate 1 = '+str(integral1)+' mHz')
	print('Total rate 2 = '+str(integral2)+' mHz')
	if density:
		n1_perBin = np.divide(n1_perBin,integral1)
		n2_perBin = np.divide(n2_perBin,integral2)
	### Plot error bars if asked for
	if drawErr:
		yerr1_perBin = np.divide(yerr1,2*xerr)
		yerr2_perBin = np.divide(yerr2,2*xerr)
		if density:
			yerr1_perBin = np.divide(yerr1_perBin,integral1)
			yerr2_perBin = np.divide(yerr2_perBin,integral2)
		if 'zenith' in variable:
			errorbar1 = ax.errorbar(bin_centers,n1,xerr=xerr,yerr=yerr1,fmt='bo',markersize=4)
			errorbar2 = ax.errorbar(bin_centers,n2,xerr=xerr,yerr=yerr2,fmt='ro',markersize=4)
		else:
			errorbar1 = ax.errorbar(bin_centers,n1_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='bo',markersize=4)
			errorbar2 = ax.errorbar(bin_centers,n2_perBin,xerr=xerr,yerr=yerr2_perBin,fmt='ro',markersize=4)
	### X axis
	if isXLog: plt.xscale('log')
	plt.tick_params(axis='x',which='both',labelsize=0)
	### Y axis
	#plt.ylim(1E-08*n2[4],2E+04*n1[4])
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	plt.tick_params(axis='y',which='both',labelsize=14)
	### Legend 
	loc = 'best'
	#leg_HNL = r'LI $\nu_{\tau} \rightarrow$HNL NC'
	#leg_HNL += '\n'
	#leg_HNL += r'(M = 1 GeV, $|U_{\tau4}|^{2}$ = $10^{-3}$)'
	#if drawErr: ax.legend([errorbar1,errorbar2],[r'GENIE $\nu_{\mu}$ CC',leg_HNL],loc=loc,frameon=False,fontsize=15)
	#if drawErr: ax.legend([errorbar1,errorbar2],[r'LI $\nu_{\mu}$ NC',r'GENIE $\nu_{\mu}$ NC'],loc=loc,frameon=False,fontsize=15)
	if drawErr: ax.legend([errorbar1,errorbar2],['Icetray LI','ICOS LI'],loc=loc,frameon=False,fontsize=15)
	#else: plt.legend([n2,n1],['NuMu events','HNL events'],loc=loc,frameon=False,fontsize=15)
	else: plt.legend(loc=loc,frameon=False,fontsize=15)
	#ax.text(0.15,0.9,'OscNext L7',transform=ax.transAxes,fontsize=15)
	#ax.text(0.15,0.85,r'$\cos\theta \leq 0$',transform=ax.transAxes,fontsize=15)
	plt.title(title)
	### Ratio plot
	ax = fig.add_subplot(spec[1,0])
	### Define histograms and errors
	n_1, bins_X, _ = plt.hist(data1,bins,weights=weights1,color='white')
	n_2, bins_X, _ = plt.hist(data2,bins,weights=weights2,color='white')
	ratio = np.divide(n1_perBin,n2_perBin,out=np.zeros_like(n1_perBin),where=n2_perBin!=0)
	yerr = np.multiply(ratio,np.sqrt(np.add(np.square(np.divide(yerr1,n_1,out=np.zeros_like(yerr1),where=n_1!=0)),np.square(np.divide(yerr2,n_2,out=np.zeros_like(yerr2),where=n_2!=0)))))
	xerr = np.zeros(len(ratio))
	for i in range(len(ratio)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	errorbar_ratio = plt.errorbar(bin_centers,ratio,xerr=xerr,yerr=yerr,fmt='o',markersize=4)
#	### Plot line at y mean
#	yMean = 0
#	total = 0
#	for i in range(len(ratio)):
#		weight_ratio = 0
#		#if yerr[i] != 0: weight_ratio = ratio[i]/yerr[i]
#		if yerr[i] != 0: weight_ratio = 1/math.sqrt(yerr[i])
#		yMean += pow(weight_ratio,2)*ratio[i]
#		total += pow(weight_ratio,2)
#	yMean = yMean/total
#	print(yMean)
	line = plt.axhline(1,color='red')
	#line = plt.axhline(1,color='red')
	### X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	plt.tick_params(axis='x',which='both',labelsize=18)
	### Y axis
	#plt.ylim(0.5*yMean,1.5*yMean)
	plt.ylim(0.1*np.min(ratio[np.nonzero(ratio)]),10*np.max(ratio[np.nonzero(ratio)]))
	plt.ylabel('Ratio')
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	#ax.yaxis.grid(True, which='minor')
	plt.minorticks_on()
	plt.tick_params(axis='y',which='both',labelsize=12)
	### Legend
	ax.legend([line],['1'],loc='best',frameon=False,fontsize=12)
	### Aesthetics  
	fig.tight_layout()
	### Save and plot
	plt.savefig(variable+'_2hist_ratio.pdf')
	plt.savefig(variable+'_2hist_ratio.png')
	plt.show()

# Function to plot three 1D histograms
def plot_3hist_1D(data1, data2, data3, label_X, label_Y, title, bins, weights1, weights2, weights3, yerr1, yerr2, yerr3, drawErr, isXLog, isYLog, variable):
	# Define canvas geometry
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	# Compute bins and x-axis error bars
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	color1 = 'black'
	color2 = 'red'
	color3 = 'blue'
	if drawErr:
		color1 = 'white'
		color2 = 'white'
		color3 = 'white'
	n1, bins, _ = plt.hist(data1,bins,weights=weights1,histtype='step',density=True,color=color1,fill=False,label=r'genie $\nu_\mu$ CC')
	n2, bins, _ = plt.hist(data2,bins,weights=weights2,histtype='step',density=True,color=color2,fill=False,label=r'LI $\nu_\mu$ CC')
	n3, bins, _ = plt.hist(data3,bins,weights=weights3,histtype='step',density=True,color=color3,fill=False,label=r'LI HNL NC')
	n1_perBin = np.divide(n1,2*xerr)
	n2_perBin = np.divide(n2,2*xerr)
	n3_perBin = np.divide(n3,2*xerr)
	# Compute total rates
	integral1 = 0
	integral2 = 0
	integral3 = 0
	for i in range(len(n1)):
		integral1 += n1[i]
		integral2 += n2[i]
		integral3 += n2[i]
	print('Total rate 1 = '+str(integral1))
	print('Total rate 2 = '+str(integral2))
	print('Total rate 3 = '+str(integral3))

	if drawErr:
		yerr1_perBin = np.divide(yerr1,2*xerr)
		yerr2_perBin = np.divide(yerr2,2*xerr)
		yerr3_perBin = np.divide(yerr3,2*xerr)
		errorbar1 = ax.errorbar(bin_centers,n1,xerr=xerr,yerr=yerr1,fmt='bo',markersize=4)
		errorbar2 = ax.errorbar(bin_centers,n2,xerr=xerr,yerr=yerr2,fmt='ro',markersize=4)
		errorbar3 = ax.errorbar(bin_centers,n3,xerr=xerr,yerr=yerr3,fmt='ko',markersize=4)

	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	#plt.ylim(1E-06,1E-01)
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(which='both',labelsize=18)
	fig.tight_layout()

	loc = 'upper right'
	if 'zenith' in variable: loc = 'upper left'
	if drawErr: ax.legend([errorbar3,errorbar2,errorbar1],['NuMu events','HNL events'],loc=loc,frameon=False,fontsize=15)
	#else: plt.legend([n2,n1],['NuMu events','HNL events'],loc=loc,frameon=False,fontsize=15)
	else: plt.legend(loc=loc,frameon=False,fontsize=15)

	plt.savefig(variable+'_3hist.pdf')
	plt.savefig(variable+'_3hist.png')

	plt.show()

# Function to plot 5 1D histograms
def plot_5hist_1D(data1, data2, data3, data4, data5, label_X, label_Y, title, bins, weights1, weights2, weights3, weights4, weights5, yerr1, yerr2, yerr3, yerr4, yerr5, isXLog, isYLog, quantity):
	plt.clf()
	#fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
	fig = plt.figure()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)

	# create grid for different subplots
	#spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[2, 1], height_ratios=[1, 2])
	spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3,1])
	### Upper plot
	ax1 = fig.add_subplot(spec[0])
	
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)):
		xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	n1, bins, _ = plt.hist(data1,bins,weights=weights1,color='white')
	n2, bins, _ = plt.hist(data2,bins,weights=weights2,color='white')
	n3, bins, _ = plt.hist(data3,bins,weights=weights3,color='white')
	n4, bins, _ = plt.hist(data4,bins,weights=weights4,color='white')
	n5, bins, _ = plt.hist(data5,bins,weights=weights5,color='white')
	n1_perBin = np.divide(n1,2*xerr)
	n2_perBin = np.divide(n2,2*xerr)
	n3_perBin = np.divide(n3,2*xerr)
	n4_perBin = np.divide(n4,2*xerr)
	n5_perBin = np.divide(n5,2*xerr)
	yerr1_perBin = np.divide(yerr1,2*xerr)
	yerr2_perBin = np.divide(yerr2,2*xerr)
	yerr3_perBin = np.divide(yerr3,2*xerr)
	yerr4_perBin = np.divide(yerr4,2*xerr)
	yerr5_perBin = np.divide(yerr5,2*xerr)
	# Compute total rates
	integral1 = np.sum(n1)
	integral2 = np.sum(n2)
	integral3 = np.sum(n3)
	integral4 = np.sum(n4)
	integral5 = np.sum(n5)
	print('Total rate 1 = '+str(integral1)+' mHz')
	print('Total rate 2 = '+str(integral2)+' mHz')
	print('Total rate 3 = '+str(integral3)+' mHz')
	print('Total rate 4 = '+str(integral4)+' mHz')
	print('Total rate 5 = '+str(integral5)+' mHz')
	errorbar1 = plt.errorbar(bin_centers,n1_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='ko',markersize=4,label=r'HNL $|U_{\tau4}|^{2} = 10^{-3}$')
	errorbar2 = plt.errorbar(bin_centers,n2_perBin,xerr=xerr,yerr=yerr2_perBin,fmt='ro',markersize=4,label=r'$\nu_{e}$ CC')
	errorbar3 = plt.errorbar(bin_centers,n3_perBin,xerr=xerr,yerr=yerr3_perBin,fmt='bo',markersize=4,label=r'$\nu_{\mu}$ CC')
	errorbar4 = plt.errorbar(bin_centers,n4_perBin,xerr=xerr,yerr=yerr4_perBin,fmt='mo',markersize=4,label=r'$\nu_{\tau}$ CC')
	errorbar5 = plt.errorbar(bin_centers,n5_perBin,xerr=xerr,yerr=yerr5_perBin,fmt='go',markersize=4,label=r'NC')
	
	### X axis
	if isXLog: plt.xscale('log')
	#plt.tick_params(axis='x',which='both',labelsize=0)
	### Y axis
	plt.ylabel(label_Y,labelpad=10)
	ax1.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	plt.tick_params(axis='y',which='both',labelsize=18)

	plt.title(title,fontsize=16,pad=20)

	### Legend
	plt.legend(loc='best',frameon=False,fontsize=12)
	#plt.legend([line],['1'],loc='best',frameon=False,fontsize=12)
	
	### Ratio plot
	#plt = fig.add_subplot(spec[1,0])
	ax2 = fig.add_subplot(spec[1], sharex=ax1)

    ### Define histograms and errors
	n_den = n2+n3+n4+n5
	yerr_den = np.sqrt(np.power(yerr2,2)+np.power(yerr3,2)+np.power(yerr4,2)+np.power(yerr5,2))
	ratio = np.divide(n1,n_den,out=np.zeros_like(n1),where=n_den!=0)
	yerr = np.multiply(ratio,np.sqrt(np.add(np.square(np.divide(yerr1,n1,out=np.zeros_like(yerr1),where=n1!=0)),np.square(np.divide(yerr_den,n_den,out=np.zeros_like(yerr_den),where=n_den!=0)))))
	xerr = np.zeros(len(ratio))
	for i in range(len(ratio)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	errorbar_ratio = plt.errorbar(bin_centers,ratio,xerr=xerr,yerr=yerr,fmt='o',markersize=4)
	
	#line = plt.plthline(1,color='red')

	### X axis
	plt.xlabel(label_X,labelpad=10)
	ax2.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	plt.tick_params(axis='x',which='both',labelsize=18)
	### Y axis
	plt.ylabel('S/B',labelpad=10)
	ax2.yaxis.label.set_size(18)
	#plt.yaxis.grid(True, which='minor')
	plt.minorticks_on()
	#plt.tick_params(axis='y',which='both',labelsize=12)
	plt.tick_params(axis='y',which='both',labelsize=14)

	### Aesthetics
	plt.setp(ax1.get_xticklabels(), visible=False)
	fig.tight_layout()
	plt.savefig(quantity+'_signalVSbg.pdf')
	plt.savefig(quantity+'_signalVSbg.png')
	plt.show()

# Function to plot 5 1D histograms
def plot_5hist_1DStack(data1, data2, data3, data4, data5, label_X, label_Y, title, bins, weights1, weights2, weights3, weights4, weights5, yerr1, yerr2, yerr3, yerr4, yerr5, isXLog, isYLog, quantity):
	plt.clf()
	#fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
	fig = plt.figure()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)

	# create grid for different subplots
	#spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[2, 1], height_ratios=[1, 2])
	spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3,1])
	### Upper plot
	ax1 = fig.add_subplot(spec[0])
	
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)):
		xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])

	nSig, bins, _ = plt.hist(data1,bins,weights=weights1,color='white')
	nBG1, bins, _ = plt.hist(data2,bins,weights=weights2,color='white')
	nBG2, bins, _ = plt.hist(data3,bins,weights=weights3,color='white')
	nBG3, bins, _ = plt.hist(data4,bins,weights=weights4,color='white')
	nBG4, bins, _ = plt.hist(data5,bins,weights=weights5,color='white')
	nSig_perBin = np.divide(nSig,2*xerr)
	nBG1_perBin = np.divide(nBG1,2*xerr)
	nBG2_perBin = np.divide(nBG2,2*xerr)
	nBG3_perBin = np.divide(nBG3,2*xerr)
	nBG4_perBin = np.divide(nBG4,2*xerr)
	yerr1_perBin = np.divide(yerr1,2*xerr)
	yerr2_perBin = np.divide(yerr2,2*xerr)
	yerr3_perBin = np.divide(yerr3,2*xerr)
	yerr4_perBin = np.divide(yerr4,2*xerr)
	yerr5_perBin = np.divide(yerr5,2*xerr)
	# Compute total rates
	integral1 = np.sum(nSig)
	integral2 = np.sum(nBG1)
	integral3 = np.sum(nBG2)
	integral4 = np.sum(nBG3)
	integral5 = np.sum(nBG4)
	print('Total rate 1 = '+str(integral1)+' mHz')
	print('Total rate 2 = '+str(integral2)+' mHz')
	print('Total rate 3 = '+str(integral3)+' mHz')
	print('Total rate 4 = '+str(integral4)+' mHz')
	print('Total rate 5 = '+str(integral5)+' mHz')
	errorbar1 = plt.errorbar(bin_centers,nSig_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='ko',markersize=4,label=r'HNL $|U_{\tau4}|^{2} = 10^{-3}$')
	colors = ['r','b','m','g']
	barStack2 = plt.bar(np.delete(bins,-1),nBG2_perBin,width=2*xerr,align='edge',color=colors[1],bottom=nBG3_perBin+nBG1_perBin+nBG4_perBin,label=r'$\nu_{\mu}$')
	barStack4 = plt.bar(np.delete(bins,-1),nBG4_perBin,width=2*xerr,align='edge',color=colors[3],bottom=nBG3_perBin+nBG1_perBin,label=r'Atmospheric muons')
	barStack1 = plt.bar(np.delete(bins,-1),nBG1_perBin,width=2*xerr,align='edge',color=colors[0],bottom=nBG3_perBin,label=r'$\nu_{e}$')
	barStack3 = plt.bar(np.delete(bins,-1),nBG3_perBin,width=2*xerr,align='edge',color=colors[2],label=r'$\nu_{\tau}$')

	# Get max and min
	minSig = np.min(nSig_perBin[np.nonzero(nSig_perBin)])
	maxBG = max(np.concatenate((nBG1_perBin,nBG2_perBin,nBG3_perBin,nBG4_perBin)))
	print(minSig,maxBG)

	### X axis
	if isXLog: plt.xscale('log')
	#plt.tick_params(axis='x',which='both',labelsize=0)
	### Y axis
	plt.ylabel(label_Y,labelpad=10)
	ax1.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	plt.tick_params(axis='y',which='both',labelsize=18)
	plt.ylim(0.1*minSig,500*maxBG)

	plt.title(title,fontsize=16,pad=20)

	### Legend
	plt.legend(loc='best',frameon=False,fontsize=12)
	#frame.set_edgecolor('black')
	#plt.legend([line],['1'],loc='best',frameon=False,fontsize=12)
	
	### Ratio plot
	#plt = fig.add_subplot(spec[1,0])
	ax2 = fig.add_subplot(spec[1], sharex=ax1)

    ### Define histograms and errors
	n_den = nBG1+nBG2+nBG3+nBG4
	yerr_den = np.sqrt(np.power(yerr2,2)+np.power(yerr3,2)+np.power(yerr4,2))
	ratio = np.divide(nSig,n_den,out=np.zeros_like(nSig),where=n_den!=0)
	yerr = np.multiply(ratio,np.sqrt(np.add(np.square(np.divide(yerr1,nSig,out=np.zeros_like(yerr1),where=nSig!=0)),np.square(np.divide(yerr_den,n_den,out=np.zeros_like(yerr_den),where=n_den!=0)))))
	xerr = np.zeros(len(ratio))
	for i in range(len(ratio)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	errorbar_ratio = plt.errorbar(bin_centers,ratio,xerr=xerr,yerr=yerr,fmt='o',markersize=4)
	
	#line = plt.plthline(1,color='red')

	### X axis
	plt.xlabel(label_X,labelpad=10)
	ax2.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	plt.tick_params(axis='x',which='both',labelsize=18)
	### Y axis
	plt.ylabel('S/B',labelpad=10)
	ax2.yaxis.label.set_size(18)
	#plt.yaxis.grid(True, which='minor')
	plt.minorticks_on()
	#plt.tick_params(axis='y',which='both',labelsize=12)
	plt.tick_params(axis='y',which='both',labelsize=14)

	### Aesthetics
	plt.setp(ax1.get_xticklabels(), visible=False)
	fig.tight_layout()
	plt.savefig(quantity+'_signalVSbg.pdf')
	plt.savefig(quantity+'_signalVSbg.png')
	plt.show()

# Function to plot 5 1D histograms
def plot_6hist_1DStack(data1, data2, data3, data4, data5, data6, label_X, label_Y, title, bins, weights1, weights2, weights3, weights4, weights5, weights6, yerr1, yerr2, yerr3, yerr4, yerr5, yerr6, isXLog, isYLog, quantity):
	plt.clf()
	#fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
	fig = plt.figure()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)

	# create grid for different subplots
	#spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[2, 1], height_ratios=[1, 2])
	spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3,1])
	### Upper plot
	ax1 = fig.add_subplot(spec[0])
	
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)):
		xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])

	nSig, bins, _ = plt.hist(data1,bins,weights=weights1,color='white')
	nBG1, bins, _ = plt.hist(data2,bins,weights=weights2,color='white')
	nBG2, bins, _ = plt.hist(data3,bins,weights=weights3,color='white')
	nBG3, bins, _ = plt.hist(data4,bins,weights=weights4,color='white')
	nBG4, bins, _ = plt.hist(data5,bins,weights=weights5,color='white')
	nBG5, bins, _ = plt.hist(data6,bins,weights=weights6,color='white')
	nSig_perBin = np.divide(nSig,2*xerr)
	nBG1_perBin = np.divide(nBG1,2*xerr)
	nBG2_perBin = np.divide(nBG2,2*xerr)
	nBG3_perBin = np.divide(nBG3,2*xerr)
	nBG4_perBin = np.divide(nBG4,2*xerr)
	nBG5_perBin = np.divide(nBG5,2*xerr)
	yerr1_perBin = np.divide(yerr1,2*xerr)
	yerr2_perBin = np.divide(yerr2,2*xerr)
	yerr3_perBin = np.divide(yerr3,2*xerr)
	yerr4_perBin = np.divide(yerr4,2*xerr)
	yerr5_perBin = np.divide(yerr5,2*xerr)
	yerr6_perBin = np.divide(yerr6,2*xerr)
	# Compute total rates
	integral1 = np.sum(nSig)
	integral2 = np.sum(nBG1)
	integral3 = np.sum(nBG2)
	integral4 = np.sum(nBG3)
	integral5 = np.sum(nBG4)
	integral6 = np.sum(nBG5)
	print('Total rate 1 = '+str(integral1)+' mHz')
	print('Total rate 2 = '+str(integral2)+' mHz')
	print('Total rate 3 = '+str(integral3)+' mHz')
	print('Total rate 4 = '+str(integral4)+' mHz')
	print('Total rate 5 = '+str(integral5)+' mHz')
	print('Total rate 6 = '+str(integral6)+' mHz')
	errorbar1 = plt.errorbar(bin_centers,nSig_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='ko',markersize=4,label=r'HNL $|U_{\tau4}|^{2} = 10^{-3}$')
	colors = ['r','b','m','g','y']
	barStack2 = plt.bar(np.delete(bins,-1),nBG2_perBin,width=2*xerr,align='edge',color=colors[1],bottom=nBG4_perBin+nBG3_perBin+nBG1_perBin+nBG5_perBin,label=r'$\nu_{\mu}$ CC')
	barStack5 = plt.bar(np.delete(bins,-1),nBG5_perBin,width=2*xerr,align='edge',color=colors[4],bottom=nBG4_perBin+nBG3_perBin+nBG1_perBin,label=r'Atmospheric muons')
	barStack1 = plt.bar(np.delete(bins,-1),nBG1_perBin,width=2*xerr,align='edge',color=colors[0],bottom=nBG4_perBin+nBG3_perBin,label=r'$\nu_{e}$ CC')
	barStack3 = plt.bar(np.delete(bins,-1),nBG3_perBin,width=2*xerr,align='edge',color=colors[2],bottom=nBG4_perBin,label=r'$\nu_{\tau}$ CC')
	barStack4 = plt.bar(np.delete(bins,-1),nBG4_perBin,width=2*xerr,align='edge',color=colors[3],label=r'NC')

	# Get max and min
	minSig = np.min(nSig_perBin[np.nonzero(nSig_perBin)])
	maxBG = max(np.concatenate((nBG1_perBin,nBG2_perBin,nBG3_perBin,nBG4_perBin,nBG5_perBin)))
	print(minSig,maxBG)

	### X axis
	if isXLog: plt.xscale('log')
	#plt.tick_params(axis='x',which='both',labelsize=0)
	### Y axis
	plt.ylabel(label_Y,labelpad=10)
	ax1.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	plt.tick_params(axis='y',which='both',labelsize=18)
	plt.ylim(0.1*minSig,500*maxBG)

	plt.title(title,fontsize=16,pad=20)

	### Legend
	plt.legend(loc='best',frameon=False,fontsize=12)
	#frame.set_edgecolor('black')
	#plt.legend([line],['1'],loc='best',frameon=False,fontsize=12)
	
	### Ratio plot
	#plt = fig.add_subplot(spec[1,0])
	ax2 = fig.add_subplot(spec[1], sharex=ax1)

    ### Define histograms and errors
	n_den = nBG1+nBG2+nBG3+nBG4+nBG5
	yerr_den = np.sqrt(np.power(yerr2,2)+np.power(yerr3,2)+np.power(yerr4,2)+np.power(yerr5,2))
	ratio = np.divide(nSig,n_den,out=np.zeros_like(nSig),where=n_den!=0)
	yerr = np.multiply(ratio,np.sqrt(np.add(np.square(np.divide(yerr1,nSig,out=np.zeros_like(yerr1),where=nSig!=0)),np.square(np.divide(yerr_den,n_den,out=np.zeros_like(yerr_den),where=n_den!=0)))))
	xerr = np.zeros(len(ratio))
	for i in range(len(ratio)): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	errorbar_ratio = plt.errorbar(bin_centers,ratio,xerr=xerr,yerr=yerr,fmt='o',markersize=4)
	
	#line = plt.plthline(1,color='red')

	### X axis
	plt.xlabel(label_X,labelpad=10)
	ax2.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	plt.tick_params(axis='x',which='both',labelsize=18)
	### Y axis
	plt.ylabel('S/B',labelpad=10)
	ax2.yaxis.label.set_size(18)
	#plt.yaxis.grid(True, which='minor')
	plt.minorticks_on()
	#plt.tick_params(axis='y',which='both',labelsize=12)
	plt.tick_params(axis='y',which='both',labelsize=14)

	### Aesthetics
	plt.setp(ax1.get_xticklabels(), visible=False)
	fig.tight_layout()
	plt.savefig(quantity+'_signalVSbg.pdf')
	plt.savefig(quantity+'_signalVSbg.png')
	plt.show()

# Function to plot 6 1D histograms
def plot_4hist_1D(n1, n2, n3, n4, label_X, label_Y, title, bins, yerr1, yerr2, yerr3, yerr4, isXLog, isYLog, particle, quantity):
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)):
		xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])

	n1_perBin = np.divide(n1,2*xerr)
	n2_perBin = np.divide(n2,2*xerr)
	n3_perBin = np.divide(n3,2*xerr)
	n4_perBin = np.divide(n4,2*xerr)
	yerr1_perBin = np.divide(yerr1,2*xerr)
	yerr2_perBin = np.divide(yerr2,2*xerr)
	yerr3_perBin = np.divide(yerr3,2*xerr)
	yerr4_perBin = np.divide(yerr4,2*xerr)

	errorbar1 = ax.errorbar(bin_centers,n1_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='ko',markersize=4)
	errorbar2 = ax.errorbar(bin_centers,n2_perBin,xerr=xerr,yerr=yerr2_perBin,fmt='ro',markersize=4)
	errorbar3 = ax.errorbar(bin_centers,n3_perBin,xerr=xerr,yerr=yerr3_perBin,fmt='bo',markersize=4)
	errorbar4 = ax.errorbar(bin_centers,n4_perBin,xerr=xerr,yerr=yerr4_perBin,fmt='mo',markersize=4)
	
	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	h_max_list = [n1.max(),n2.max(),n3.max(),n4.max()]
	h_min_list = [n1.min(),n2.min(),n3.min(),n4.min()]
	#plt.ylim(5e-07,1e+02)
	#plt.ylim(5e-02,1e+01)
	#if quantity == 'Energy': plt.ylim(1E-13,5)
	#elif quantity == 'Zenith': plt.ylim(1E-06,5E-02)
	#if quantity == 'Energy': plt.ylim(5E-12,5E-02)
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(labelsize=18)
	plt.title(title)
	plt.title(title,fontsize=16,pad=20)
	plt.legend([errorbar1,errorbar2,errorbar3,errorbar4],['True-True','True-False','False-True','False-False'],loc='best',frameon=False,fontsize=12,ncol=2)
	
	#ax.text(0.1,0.9,'LeptonInjector',transform=ax.transAxes,fontsize=15)
	
	fig.tight_layout()
	plt.savefig(particle+'_'+quantity+'_4histos.pdf')
	plt.savefig(particle+'_'+quantity+'_4histos.png')
	plt.show()

# Function to plot 6 1D histograms
def plot_9hist_1D(data1, data2, data3, data4, data5, data6, data7, data8, data9, label_X, label_Y, title, bins, weights1, weights2, weights3, weights4, weights5, weights6, weights7, weights8, weights9, yerr1, yerr2, yerr3, yerr4, yerr5, yerr6, yerr7, yerr8, yerr9, isXLog, isYLog, particle, quantity):
	plt.clf()
	fig, ax = plt.subplots()
	fig.set_figheight(7)
	fig.set_figwidth(7.5)
	
	xerr = np.zeros(len(bins)-1)
	for i in range(len(xerr)):
		xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	n1, bins, _ = plt.hist(data1,bins,weights=weights1,label='',color='white')
	n2, bins, _ = plt.hist(data2,bins,weights=weights2,label='',color='white')
	n3, bins, _ = plt.hist(data3,bins,weights=weights3,label='',color='white')
	n4, bins, _ = plt.hist(data4,bins,weights=weights4,label='',color='white')
	n5, bins, _ = plt.hist(data5,bins,weights=weights5,label='',color='white')
	n6, bins, _ = plt.hist(data6,bins,weights=weights6,label='',color='white')
	n7, bins, _ = plt.hist(data7,bins,weights=weights7,label='',color='white')
	n8, bins, _ = plt.hist(data8,bins,weights=weights8,label='',color='white')
	n9, bins, _ = plt.hist(data9,bins,weights=weights9,label='',color='white')
	#n1_perBin = np.divide(n1,2*xerr)
	#n2_perBin = np.divide(n2,2*xerr)
	#n3_perBin = np.divide(n3,2*xerr)
	#n4_perBin = np.divide(n4,2*xerr)
	#n5_perBin = np.divide(n5,2*xerr)
	#n6_perBin = np.divide(n6,2*xerr)
	#n7_perBin = np.divide(n7,2*xerr)
	#n8_perBin = np.divide(n8,2*xerr)
	#n9_perBin = np.divide(n9,2*xerr)
	#yerr1_perBin = np.divide(yerr1,2*xerr)
	#yerr2_perBin = np.divide(yerr2,2*xerr)
	#yerr3_perBin = np.divide(yerr3,2*xerr)
	#yerr4_perBin = np.divide(yerr4,2*xerr)
	#yerr5_perBin = np.divide(yerr5,2*xerr)
	#yerr6_perBin = np.divide(yerr6,2*xerr)
	#yerr7_perBin = np.divide(yerr7,2*xerr)
	#yerr8_perBin = np.divide(yerr8,2*xerr)
	#yerr9_perBin = np.divide(yerr9,2*xerr)
	n1_perBin = n1/n1
	n2_perBin = n2/n1 
	n3_perBin = n3/n1 
	n4_perBin = n4/n1 
	n5_perBin = n5/n1 
	n6_perBin = n6/n1 
	n7_perBin = n7/n1 
	n8_perBin = n8/n1 
	n9_perBin = n9/n1 
	yerr1_perBin = n1_perBin*np.sqrt(pow(yerr1/n1,2)+pow(yerr1/n1,2))
	yerr2_perBin = n2_perBin*np.sqrt(pow(yerr2/n2,2)+pow(yerr1/n1,2))
	yerr3_perBin = n3_perBin*np.sqrt(pow(yerr3/n3,2)+pow(yerr1/n1,2))
	yerr4_perBin = n4_perBin*np.sqrt(pow(yerr4/n4,2)+pow(yerr1/n1,2))
	yerr5_perBin = n5_perBin*np.sqrt(pow(yerr5/n5,2)+pow(yerr1/n1,2))
	yerr6_perBin = n6_perBin*np.sqrt(pow(yerr6/n6,2)+pow(yerr1/n1,2))
	yerr7_perBin = n7_perBin*np.sqrt(pow(yerr7/n7,2)+pow(yerr1/n1,2))
	yerr8_perBin = n8_perBin*np.sqrt(pow(yerr8/n8,2)+pow(yerr1/n1,2))
	yerr9_perBin = n9_perBin*np.sqrt(pow(yerr9/n9,2)+pow(yerr1/n1,2))
	# Compute total rates
	integral1 = np.sum(n1_perBin)/(len(bins)-1) 
	integral2 = np.sum(n2_perBin)/(len(bins)-1)
	integral3 = np.sum(n3_perBin)/(len(bins)-1)
	integral4 = np.sum(n4_perBin)/(len(bins)-1)
	integral5 = np.sum(n5_perBin)/(len(bins)-1)
	integral6 = np.sum(n6_perBin)/(len(bins)-1)
	integral7 = np.sum(n7_perBin)/(len(bins)-1)
	integral8 = np.sum(n8_perBin)/(len(bins)-1)
	integral9 = np.sum(n9_perBin)/(len(bins)-1)
	print('Total Gen rate = '+str(integral1))
	print('Total Phot rate = '+str(integral2))
	print('Total Det rate = '+str(integral3))
	print('Total L2 rate = '+str(integral4))
	print('Total L3 rate = '+str(integral5))
	print('Total L4 rate = '+str(integral6))
	print('Total L5 rate = '+str(integral7))
	print('Total L6 rate = '+str(integral8))
	print('Total L7 rate = '+str(integral9))
	errorbar1 = ax.errorbar(bin_centers,n1_perBin,xerr=xerr,yerr=yerr1_perBin,fmt='ko',markersize=4)
	errorbar2 = ax.errorbar(bin_centers,n2_perBin,xerr=xerr,yerr=yerr2_perBin,fmt='ro',markersize=4)
	errorbar3 = ax.errorbar(bin_centers,n3_perBin,xerr=xerr,yerr=yerr3_perBin,fmt='bo',markersize=4)
	errorbar4 = ax.errorbar(bin_centers,n4_perBin,xerr=xerr,yerr=yerr4_perBin,fmt='mo',markersize=4)
	errorbar5 = ax.errorbar(bin_centers,n5_perBin,xerr=xerr,yerr=yerr5_perBin,fmt='go',markersize=4)
	errorbar6 = ax.errorbar(bin_centers,n6_perBin,xerr=xerr,yerr=yerr6_perBin,fmt='co',markersize=4)
	errorbar7 = ax.errorbar(bin_centers,n7_perBin,xerr=xerr,yerr=yerr7_perBin,fmt='yo',markersize=4)
	errorbar8 = ax.errorbar(bin_centers,n8_perBin,xerr=xerr,yerr=yerr8_perBin,fmt='bo',markersize=4)
	errorbar9 = ax.errorbar(bin_centers,n9_perBin,xerr=xerr,yerr=yerr9_perBin,fmt='ro',markersize=4)
	
	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	h_max_list = [n1.max(),n2.max(),n3.max(),n4.max()]
	h_min_list = [n1.min(),n2.min(),n3.min(),n4.min()]
	plt.ylim(5e-07,1e+02)
	#plt.ylim(5e-02,1e+01)
	#if quantity == 'Energy': plt.ylim(1E-13,5)
	#elif quantity == 'Zenith': plt.ylim(1E-06,5E-02)
	#if quantity == 'Energy': plt.ylim(5E-12,5E-02)
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(labelsize=18)
	plt.title(title)
	plt.title(title,fontsize=16,pad=20)
	#plt.legend([errorbar1,errorbar2,errorbar3,errorbar4,errorbar5,errorbar6,errorbar7,errorbar8,errorbar9],['L2 (rate = '+str(round(1000*integral1,2))+' mHz)','L3 (rate = '+str(round(1000*integral2,2))+' mHz)','L4 (rate = '+str(round(1000*integral3,2))+' mHz)','L5 (rate = '+str(round(1000*integral4,2))+' mHz)','L6 (rate = '+str(round(1000*integral5,2))+' mHz)','L7 (rate = '+str(round(1000*integral6,2))+' mHz)'],loc=loc,frameon=False,fontsize=12)
	plt.legend([errorbar1,errorbar2,errorbar3,errorbar4,errorbar5,errorbar6,errorbar7,errorbar8,errorbar9],['Gen','Phot','Det','L2','L3','L4','L5','L6','L7'],loc='best',frameon=False,fontsize=12,ncol=3)
	
	#ax.text(0.1,0.9,'LeptonInjector',transform=ax.transAxes,fontsize=15)
	#ax.text(0.1,0.9,r'GENIE $\nu_{\mu}$ CC',transform=ax.transAxes,fontsize=15)
	#ax.text(0.1,0.85,'OscNext selection',transform=ax.transAxes,fontsize=15)
	
	fig.tight_layout()
	plt.savefig(particle+'_'+quantity+'_9histos.pdf')
	plt.savefig(particle+'_'+quantity+'_9histos.png')
	plt.show()

def plot_2graphs_1D(n1, n2, err1, err2, bins, label_X, title, isYLog, particle, quantity_X, quantity_Y):
	plt.clf()
	xerr = np.zeros(len(bins)-1)
	for i in range(len(bins)-1): xerr[i] = 0.5*(bins[i+1] - bins[i])
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	plt.errorbar(bin_centers,n1,xerr=xerr,yerr=err1,fmt='bo')
	plt.errorbar(bin_centers,n2,xerr=xerr,yerr=err2,fmt='ro')
	plt.xlabel(label_X)
	plt.ylabel('Mean '+particle+' '+quantity_Y)
	plt.title(title)
	if isYLog: plt.yscale('log')
	plt.savefig('mean-'+particle+'-'+quantity_Y+'_as_'+quantity_X+'_overlap.pdf')
	plt.show()

# Function to plot 2D histograms
def plot_2D(h, mean, label_X, label_Y, label_Z, title, bins_X, bins_Y, isXLog, isYLog, isZLog, quantity_1, quantity_2):
#def plot_2D(h, label_X, label_Y, label_Z, title, bins_X, bins_Y, isXLog, isYLog, isZLog, quantity_1, quantity_2):
	plt.clf()
	fig = plt.figure()
	fig.set_figheight(7.5)
	fig.set_figwidth(9.5)
	ax = fig.subplots()
	# Remove bins with content < 1
	#h_cleaned = []
	#for i in range(len(h)):
	#	h_cleaned.append(h[i])
	#	for j in range(len(h[i])):
	#		if h[i][j] >= 1: h_cleaned[i][j] = h[i][j]
	#		else: h_cleaned[i][j] = 0

	#h = np.asarray(h_cleaned)

	#h_min = 0
	h_min = np.min(h[np.nonzero(h)])
	h_max = h.max()

	print(h_min, h_max)

	# Mask zeros
	h = np.ma.masked_where(h==0,h)
	#h = np.ma.masked_where(h<1,h)

	if isZLog: norm = colors.LogNorm(vmin=h_min, vmax=h_max)
	else: norm = colors.Normalize(vmin=h_min,vmax=h_max)
	#norm = MidpointNormalize(midpoint=1,vmin=h_min, vmax=h_max)
	pcm = ax.pcolormesh(bins_X, bins_Y, h, norm=norm, cmap='jet')#cmap='Greens')#, cmap='Wistia')
	#pcm = ax.pcolormesh(bins_X, bins_Y, h, cmap='Greens')#, cmap='Wistia')
	# To change colormap: pcm = ax.pcolormesh(bins_X, bins_Y, h, vmin=h.min(), vmax=h.max(), cmap='PuBu_r')
	# Show Z value on each X-Y bin
	#for i in range(len(bins_Y)-1):
	#    for j in range(len(bins_X)-1): ax.text(bins_X[j]+0.025,(bins_Y[i]+bins_Y[i+1])/2, round(h[i,j],2), color="k", ha="center", va="center", fontweight="bold", fontsize=7)

	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	#plt.xticks(np.arange(0,2,1),['False','True','True'],horizontalalignment='left')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	#plt.yticks(np.arange(0,2,1),['False','True','True'],horizontalalignment='left',rotation=90)
	# Z axis
	cbar = plt.colorbar(pcm, ax=ax)
	cbar.ax.set_ylabel(label_Z,fontsize=18,labelpad=12)
	cbar.ax.tick_params(labelsize=15)

	#ax.text(0.6,0.9,'GENIE (post-L2)',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.85,'7.5 years',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.3,'LeptonInjector (post-L2)',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.25,'7.5 years',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.2,r'$M = 1$ GeV',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.15,r'$|U_{\tau4}|^{2} = 10^{-3}$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.9,r'$|Cascade_{1}^{Z}| < 600 m$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.7,0.9,r'$|U_{\tau4}|^{2} = 10^{-1}$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.7,0.8,r'PID$>$0.8',transform=ax.transAxes,fontsize=15)

	## Mean line
	ax.stairs(mean,bins_X,baseline=None,color='k',lw=3)
	
	plt.tick_params(labelsize=18)
	plt.title(title,fontsize=18,pad=20)

	fig.tight_layout()

	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D.png')
	plt.show()

#def plot_2D_scatter(data1, data2, mean, binsX, label_X, label_Y, title, isXLog, isYLog, quantity_1, quantity_2):
def plot_2D_scatter(data1, data2, label_X, label_Y, title, isXLog, isYLog, quantity_1, quantity_2):
	plt.clf()
	fig = plt.figure()
	fig.set_figheight(8)
	fig.set_figwidth(9.5)
	ax = fig.subplots()

	# Scatter plot
	ax.scatter(data1,data2,s=1)
	## Mean line
	#ax.stairs(mean,binsX,baseline=None,color='r',lw=2)

	# X axis
	plt.xlabel(label_X,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	#plt.xticks(np.arange(0,2,1),['False','True','True'],horizontalalignment='left')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	#plt.yticks(np.arange(0,2,1),['False','True','True'],horizontalalignment='left',rotation=90)

	#ax.text(0.6,0.9,'GENIE (post-L2)',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.85,'7.5 years',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.3,'LeptonInjector (post-L2)',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.25,'7.5 years',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.2,r'$M = 1$ GeV',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.15,r'$|U_{\tau4}|^{2} = 10^{-3}$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.6,0.9,r'$|Cascade_{1}^{Z}| < 600 m$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.7,0.9,r'$|U_{\tau4}|^{2} = 10^{-1}$',transform=ax.transAxes,fontsize=15)
	#ax.text(0.7,0.8,r'PID$>$0.8',transform=ax.transAxes,fontsize=15)
	
	plt.tick_params(labelsize=18)
	plt.title(title,fontsize=18,pad=20)

	fig.tight_layout()

	plt.savefig(quantity_1+'_'+quantity_2+'_scatter.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_scatter.png')
	plt.show()

def plot_2D_text(h, title, bins_X, bins_Y, quantity_1, quantity_2):
	fig, ax = plt.subplots()

	h_min = np.min(h[np.nonzero(h)])
	h_max = h.max()
	norm = colors.Normalize(vmin=h_min,vmax=h_max)
	norm = colors.LogNorm(vmin=h_min, vmax=h_max)
	im = ax.imshow(h, norm=norm, cmap='Reds')
	# We want to show all ticks...
	ax.set_xticks(np.arange(2))
	ax.set_yticks(np.arange(2))
	# ... and label them with the respective list entries
	ax.set_xticklabels(bins_X)
	ax.set_yticklabels(bins_Y)

	plt.xlabel('Cascade 0 in DeepCore',labelpad=10)
	ax.xaxis.label.set_size(15)
	plt.ylabel('Cascade 1 in DeepCore',labelpad=10)
	ax.yaxis.label.set_size(15)
	
	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
	         rotation_mode="anchor")
	
	# Loop over data dimensions and create text annotations.
	for i in range(len(bins_X)):
	    for j in range(len(bins_Y)):
	        text = ax.text(j, i, np.format_float_scientific(h[i, j],precision=2),
	                       ha="center", va="center", color="k", fontsize=15)
	
	ax.set_title(title)
	fig.tight_layout()

	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D-text.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D-text.png')
	plt.show()

# Function to plot 2D histograms with contour
def plot_2D_contour(h, label_X, label_Y, label_Z, title, bins_X, bins_Y, isXLog, isYLog, isZLog, quantity_1, quantity_2):
	plt.clf()
	fig = plt.figure()
	fig.set_figheight(8)
	fig.set_figwidth(9.5)
	ax = fig.subplots()

	h_min = np.min(h[np.nonzero(h)])
	h_max = h.max()

	print(h_min, h_max)
	
	# Mask zeros
	h = np.ma.masked_where(h==0,h)
	#h = np.ma.masked_where(h<1,h)

	if isZLog: norm = colors.LogNorm(vmin=h_min, vmax=h_max)
	else: norm = colors.Normalize(vmin=h_min,vmax=h_max)
	#norm = MidpointNormalize(midpoint=1,vmin=h_min, vmax=h_max)
	pcm = ax.pcolormesh(bins_X, bins_Y, h, norm=norm, cmap='jet')

	#levels = np.logspace(np.log10(h_min),np.log10(h_max),40)#[1e-07,1e-06,1e-05,1e-04,1e-03,1e-02]
	#levels = np.logspace(np.log10(4e-06),np.log10(1e-02),50)#[1e-07,1e-06,1e-05,1e-04,1e-03,1e-02]
	levels_exp = np.arange(np.floor(np.log10(h_min)-1),np.ceil(np.log10(h_max)+1))
	levels = np.power(10, levels_exp)
	#CS = ax.contourf(np.delete(bins_X,-1),np.delete(bins_Y,-1),h,levels,cmap='Greens',norm=colors.LogNorm())
	#CS = ax.contourf(np.delete(bins_X,-1),np.delete(bins_Y,-1),h,cmap='Greens',norm=colors.LogNorm(vmin=h_min, vmax=h_max))

	#CS2 = ax.contour(np.delete(bins_X,-1),np.delete(bins_Y,-1),h,levels=1.0,colors='b',norm=colors.LogNorm())
	CS2 = ax.contour(np.delete(bins_X,-1),np.delete(bins_Y,-1),h,levels=[4e-06],colors='b',norm=norm)
	#CS2 = ax.contour(np.delete(bins_X,-1),np.delete(bins_Y,-1),h,levels=5.3e-06,colors='b',norm=colors.LogNorm())
	
	# X axis
	plt.xlabel(label_X)#,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	# Z axis
	cbar = plt.colorbar(pcm, ax=ax)
	cbar.ax.set_ylabel(label_Z,fontsize=18,labelpad=12)
	cbar.ax.tick_params(labelsize=15)
	#cbar = fig.colorbar(CS)
	#cbar.ax.set_ylabel(label_Z,fontsize=18)
	#cbar.ax.tick_params(labelsize=15)

	cbar.add_lines(CS2)
	#ax.text(0.3,0.85,'> 1 event in 8 years',transform=ax.transAxes,fontsize=15)
	
	plt.tick_params(which='major',labelsize=18)
	plt.tick_params(axis='x',which='minor',labelbottom=False,labelsize=0)
	plt.title(title,fontsize=18,pad=20)
	
	fig.tight_layout()
	
	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D_contour.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_hist2D_contour.png')
	plt.show()

# Function to plot 2D histograms with contour for sensitivity
def plot_2D_contour_5files(h, label_X, label_Y, title, bins_X, bins_Y, isXLog, isYLog, quantity_1, quantity_2):
	plt.clf()
	fig = plt.figure()
	fig.set_figheight(8)
	fig.set_figwidth(8.5)
	ax = fig.subplots()

	CS = []
	labels = ['L3','L4','L5','L6','L7']
	color = ['k','b','r','g','c']
	i = 0
	for hist in h:
		CSi = ax.contour(np.delete(bins_X,-1),np.delete(bins_Y,-1),hist,levels=1.0,norm=colors.LogNorm(),colors=color[i])
		CS.append(CSi)
		i += 1

	for i in range(len(CS)):
		CS[i].collections[0].set_label(labels[i])

	# X axis
	plt.xlabel(label_X)#,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(which='major',labelsize=18)
	plt.tick_params(axis='x',which='minor',labelbottom=False,labelsize=0)
	plt.title(title,fontsize=18,pad=20)

	leg = plt.legend(loc='lower right',frameon=False,fontsize=15)
	leg.set_title('> 1 event/8 years contour',prop={'size':'xx-large'})
	
	fig.tight_layout()
	
	plt.savefig(quantity_1+'_'+quantity_2+'_countour-5files.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_countour-5files.png')
	plt.show()

# Function to plot 2D histograms with contour for sensitivity
def plot_2D_sensitivity(h, label_X, label_Y, title, bins_X, bins_Y, bg, isXLog, isYLog, quantity_1, quantity_2):
	plt.clf()
	fig = plt.figure()
	fig.set_figheight(8)
	fig.set_figwidth(8.5)
	ax = fig.subplots()

	CS = []
	labels = []
	color = ['k','b','r']
	i = 0
	for hist in h:
		CSi = ax.contour(np.delete(bins_X,-1),np.delete(bins_Y,-1),hist,levels=0.05,norm=colors.LogNorm(),colors=color[i])
		CS.append(CSi)
		labels.append('BG = '+str(bg[i])+' event(s)')
		i += 1

	for i in range(len(CS)):
		CS[i].collections[0].set_label(labels[i])

	# X axis
	plt.xlabel(label_X)#,labelpad=10)
	ax.xaxis.label.set_size(18)
	if isXLog: plt.xscale('log')
	# Y axis
	plt.ylabel(label_Y)
	ax.yaxis.label.set_size(18)
	if isYLog: plt.yscale('log')
	
	plt.tick_params(which='major',labelsize=18)
	plt.tick_params(axis='x',which='minor',labelbottom=False,labelsize=0)
	plt.title(title,fontsize=18,pad=20)

	leg = plt.legend(loc='lower right',frameon=False,fontsize=15)
	leg.set_title(r'95% CL',prop={'size':'xx-large'})
	
	fig.tight_layout()
	
	plt.savefig(quantity_1+'_'+quantity_2+'_sensitivity.pdf')
	plt.savefig(quantity_1+'_'+quantity_2+'_sensitivity.png')
	plt.show()
