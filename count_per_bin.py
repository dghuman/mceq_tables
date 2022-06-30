import matplotlib.pyplot as plt
import numpy as np
import pickle
import os, sys
from colour import Color
import logging

#logging.basicConfig(filename='livetime.log', encoding='utf-8', level=logging.DEBUG)

# Load files using pickle
with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_dict_0006.pkl', 'rb') as input1:
    muon_dict = pickle.load(input1)

with open('/data/p-one/dghuman/simulation/MCEq_tables/rate_muon_table_dtheta_10.pkl', 'rb') as input2:
    rate = pickle.load(input2)

with open('/data/p-one/dghuman/simulation/MCEq_tables/mute_rate_muon_table_dtheta_10.pkl', 'rb') as input3:
    mute_rate = pickle.load(input3)
    
# Save file function using pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)
    
def find_closest(arr, val):
    index = np.argmin(np.abs(np.array(arr)-val))
    return index

theta_bin_diff = 10    
e_grid = rate['e_grid']
e_bins = rate['e_bins']
e_widths = rate['e_widths']
e_diff = e_grid[1:] - e_grid[:-1]

mute_e_grid = mute_rate['e_grid']
bin_thetas = [x for x in rate.keys() if not isinstance(x, str)]
live_time = np.zeros(len(bin_thetas))
counts = np.zeros(len(bin_thetas))

muon_Energy = np.array(muon_dict['Energy'])
muon_thetas = np.array(muon_dict['theta'])
muon_phis = np.array(muon_dict['phi'])
new_theta = []
time_weight = []
zeros = []

new_dict = {}
new_dict['Energy'] = []
new_dict['theta'] = []
new_dict['phi'] = []
new_dict['run_time'] = []
new_dict['mute_run_time'] = []
new_dict['time_length'] = []
new_dict['mute_time_length'] = []
new_dict['E_bin'] = []
new_dict['theta_bin'] = []
new_dict['e_grid'] = e_grid
new_dict['e_bin_width'] = []


run_time = 0
mute_run_time = 0

'''
    if energy_ == 0:
        energy_bin_width = e_grid[energy_+1] - e_grid[energy_]
    elif energy_ == len(muon_Energy)-1:
        energy_bin_width = e_grid[energy_] - e_grid[energy_-1]
    else:
        energy_bin_width = (e_grid[energy_+1] - e_grid[energy_-1])/2
'''


for i in range(len(muon_Energy)):
    theta_bin = find_closest(bin_thetas, muon_thetas[i])
    theta_ = bin_thetas[theta_bin]
    energy_ = find_closest(e_grid, muon_Energy[i])
    mute_energy_ = find_closest(mute_e_grid, muon_Energy[i])
    rate_ = rate[theta_][energy_]
    mute_rate_ = mute_rate[theta_][mute_energy_]
    energy_bin_width = e_widths[energy_]
    if rate_ < 1E-17:
        #time_ = 0
        zeros.append(muon_Energy[i])
    else:
        time_ = 1/(energy_bin_width*rate_) #1/(muon_Energy[i]*rate_)
    if mute_rate_ < 1E-17:
        mute_time_ = 0
    else:
        mute_time_ = 1/(energy_bin_width*mute_rate_) #(muon_Energy[i]*mute_rate_)
    time_ = mute_time_                                             # <-----------Remove this if not using mute result
    run_time += time_
    mute_run_time += mute_time_
    new_dict['e_bin_width'].append(energy_bin_width)
    new_dict['run_time'].append(run_time)
    new_dict['mute_run_time'].append(mute_run_time)
    new_dict['time_length'].append(time_)
    new_dict['mute_time_length'].append(mute_time_)
    new_dict['Energy'].append(muon_Energy[i])
    new_dict['E_bin'].append(energy_)
    new_dict['theta_bin'].append(theta_bin)
    new_dict['theta'].append(muon_thetas[i])
    new_dict['phi'].append(muon_phis[i])
    time_weight.append(time_)
    new_theta.append(muon_thetas[i])
    live_time[theta_bin] += time_
    counts[theta_bin] += 1

new_dict['zeros'] = zeros    
info = {}
info['live_time'] = live_time
info['counts'] = counts
info['thetas'] = bin_thetas
info['raw_thetas'] = new_theta
info['time_weights'] = time_weight

print('Avg livetime for first bin from mceq: ' + str(sum(new_dict['time_length'])/len(new_dict['time_length'])))
save_object(info, "/data/p-one/dghuman/simulation/MCEq_tables/live_time_zenith.pkl")
save_object(new_dict, "/data/p-one/dghuman/simulation/MCEq_tables/new_muon_dict.pkl")
