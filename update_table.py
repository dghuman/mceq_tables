import matplotlib.pyplot as plt
import numpy as np
import pickle
import os, sys
from colour import Color

#import solver related modules
from MCEq.core import MCEqRun
import mceq_config as config
from MCEq.geometry.density_profiles import GeneralizedTarget
#import primary model choices
import crflux.models as pm

# Load file using pickle
with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_table_dtheta_10.pkl', 'rb') as input:
    table = pickle.load(input)

with open('/data/p-one/dghuman/simulation/MCEq_tables/UofA_tables/fluxtable_as_list_1460.pkl', 'rb') as input4:
    uofa_table = pickle.load(input4)

# Save file function using pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)

# Convert to rate functions
def Area(r, theta):
    return np.pi*r*r*np.cos(theta)

def SolidAngle(theta, dtheta):
    return 2*np.pi*(np.cos(theta - dtheta/2) - np.cos(theta + dtheta/2))

rate_table = {}
uofa_rate_table = {}
uofa_e_grid = uofa_table['e_grid']
e_grid = table['e_grid']
l = [x for x in table.keys() if not isinstance(x, str)]
thetas = np.array(l)*np.pi/180

for key in l:
    area_ = Area(700E2, key*np.pi/180.) #in cm2 
    solidangle_ = SolidAngle(key*np.pi/180., 10*np.pi/180.)
    rate_table[key] = area_*solidangle_*np.array(table[key])
    uofa_rate_table[key] = area_*solidangle_*np.array(uofa_table[key])

rate_table['e_grid'] = e_grid
rate_table['e_bins'] = table['e_bins']
rate_table['e_widths'] = table['e_widths']
uofa_rate_table['e_grid'] = uofa_e_grid

save_object(rate_table, "/data/p-one/dghuman/simulation/MCEq_tables/rate_muon_table_dtheta_10.pkl")
save_object(uofa_rate_table, "/data/p-one/dghuman/simulation/MCEq_tables/mute_rate_muon_table_dtheta_10.pkl")

