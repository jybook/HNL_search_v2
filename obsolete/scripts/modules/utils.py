import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import tables

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# Function to extract column from 2D array
def column(matrix, i):
    return [row[i] for row in matrix]

# Function to calculate array of errors for histogram with weighted events
def histogram_errors(values, weights, bins):
    errors = np.zeros(len(bins)-1)
    for i in range(len(values)):
        for j in range(len(bins)-1):
            if values[i] >= bins[j] and values[i] < bins[j+1]: errors[j] += pow(weights[i],2)
    return np.sqrt(errors)

class import_MC(object):
    def __init__(self, infiles, MC, desired_keys='', scale=1.): 
        self.scale = scale
        self.infiles    = infiles
        self.MC         = MC
        if 'Multisim' in infiles[0]:
            print('Looks like a Multisim file.')
            self.Multisim = True
        else:
            self.Multisim = False
        self.Multisim = False
        self.desired_keys = desired_keys
        self.allocate_variables(infiles,desired_keys,self.Multisim)
    def allocate_variables(self, infiles, desired_keys, Multisim):
        variable_dict = {}
        f = tables.open_file(infiles[0])
        for i,j in f.root.keysDict.coldescrs.iteritems():
            if desired_keys:
                    if i in desired_keys:
                        variable_dict[i] = []
            else:
                variable_dict[i] = []
        MultisimAmplitudesArray = []
        MultisimPhasesArray = []
        MultisimModesArray = []
        for infile in infiles:
            if self.MC == 'True':
                scale_factor = float(infile.split('_')[-1].split('.h5')[0])/1000./self.scale
                self.livetime = scale_factor
            else:
                scale_factor = float(infile.split('_')[-1].split('.h5')[0])/1000./self.scale
                self.livetime = scale_factor
            f = tables.open_file(infile) 
            try:
                all_charges  = f.root.charges.cols.item[:]
                charge_array = [all_charges[x:x+60] for x in xrange(0, len(all_charges), 60)]
                self.charges = np.sum(charge_array, axis=0)
            except:
                pass
            if Multisim:
                ##print(f.root.MultisimAmplitudes.cols.item[:])
                #print(np.shape(f.root.MultisimAmplitudes.cols.item[:]))
                MultisimAmplitudesArray.extend(f.root.MultisimAmplitudes.cols.item[:])
                MultisimPhasesArray.extend(f.root.MultisimPhases.cols.item[:])
                MultisimModesArray.extend(f.root.MultisimModes.cols.item[:])
            for i,j in f.root.keysDict.coldescrs.iteritems():
                if desired_keys:
                    if i == 'weights':
                        variable_dict['weights'] += np.ndarray.tolist(f.root.keysDict[:]['weights']/scale_factor)
                    elif i in desired_keys:
                        variable_dict[i] += np.ndarray.tolist(f.root.keysDict[:][i])
                else:
                    if i == 'weights':
                        variable_dict['weights'] += np.ndarray.tolist(f.root.keysDict[:]['weights']/scale_factor)
                    else:
                        variable_dict[i] += np.ndarray.tolist(f.root.keysDict[:][i])
            for i,j in f.root.keysDict.coldescrs.iteritems():
                if desired_keys:
                    if i == 'weights':
                        if self.MC == False:
                            print('Uploading Data')
                            setattr(self, 'weights', np.ones(len(variable_dict['weights'])))
                        else:
                            setattr(self, 'weights', np.asarray(variable_dict['weights']))
                    elif i in desired_keys:
                        setattr(self, i, np.asarray(variable_dict[i]))
                else:
                    if i == 'weights':
                        setattr(self, 'weights', np.asarray(variable_dict['weights']))
                    else:
                        setattr(self, i, np.asarray(variable_dict[i]))
        if self.Multisim:
            len_of_array = int(len(MultisimAmplitudesArray)/len(self.weights))
            self.MultisimAmplitudes = np.asarray([MultisimAmplitudesArray[x:x+len_of_array] for x in xrange(0, len(MultisimAmplitudesArray), len_of_array)])
            self.MultisimPhases     = np.asarray([MultisimPhasesArray[x:x+len_of_array] for x in xrange(0, len(MultisimPhasesArray), len_of_array)])
            self.MultisimModes      = np.asarray([MultisimModesArray[x:x+len_of_array] for x in xrange(0, len(MultisimModesArray), len_of_array)])
        self.rate = sum(self.weights)
        print('Rate: '+ str(sum(self.weights))+' mHz +/- '+ str(np.sqrt(sum(np.power(self.weights,2))))+' mHz')
