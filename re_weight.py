## add weights based off of livetime
import numpy as np
import pickle
import matplotlib.pyplot as plt

import logging

logging.basicConfig(filename='/home/users/ghuman/log/weighting.log', level=logging.INFO)

# Save file function using pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)

def main():
    # Load files using pickle
    with open('/data/p-one/dghuman/simulation/MCEq_tables/new_muon_dict.pkl', 'rb') as input1:
        muon_dict = pickle.load(input1)
        
    e_grid = np.array(muon_dict['e_grid'])
    e_bins = np.array(muon_dict['E_bin'])
    theta = np.array(muon_dict['theta'])
    theta_bins = np.array(muon_dict['theta_bin'])

    weights = np.zeros(len(theta))

    data = {}
    min_times = []
    count = 0

    for theta_i in range(10):
        indices = np.where(theta_bins == theta_i)[0]
        logging.info(f'------------------ Starting theta bin {theta_i} ---------------------')
        if len(indices) == 0:
            data[theta_i] = {}
            data[theta_i]['mintime'] = 0
            data[theta_i]['e_weight'] = []
            data[theta_i]['indices'] = []
            min_times.append(0)
            continue
        
        hist_E, bin_edges_E = np.histogram(e_bins[indices], bins=range(len(e_grid)))
        live_time_hist = []
        E_index = []
        min_time = 0

        logging.info(f'Energy spectrum count (per energy bin) = {hist_E}')
        
        for index, val in enumerate(hist_E):
            i_val = np.where(e_bins[indices]==index)[0]
            if len(i_val) != val:
                count += 1
            if len(i_val) == 0:
                live_time_hist.append(-1)
                E_index.append(index)
                continue
            tscale = np.array(muon_dict['time_length'])[indices][i_val[0]]  # Get first occurance of timescale (they should all be the same for this energybin)
            live_time_hist.append(tscale*val)
            E_index.append(index)
            if tscale == 0 or val == 0:
                live_time_hist[-1] = -1.
                continue
            elif min_time == 0:
                min_time = tscale*val
            elif tscale*val < min_time:
                min_time = tscale*val
        
        logging.info(f'Livetime per energy bin = {live_time_hist}')
        #min_time = 1496.9781737661062 # TEMPORARY CHECK
        e_weight = min_time/np.array(live_time_hist)
        e_weight[np.where(np.array(live_time_hist) < 0)] = 0
        min_times.append(min_time)
        live_time = np.array(live_time_hist)
        live_time[np.where(live_time < 0)] = 0
        data[theta_i] = {}
        data[theta_i]['E_index'] = np.array(E_index)
        data[theta_i]['livetime'] = live_time
        data[theta_i]['mintime'] = min_time
        data[theta_i]['e_weight'] = e_weight
        data[theta_i]['indices'] = indices
        

    min_times = np.array(min_times)
    non_zero = np.where(min_times != 0)
    zero_ind = np.where(min_times == 0)
    data['mintimes'] = min_times[:-1]

    sorted_mintimes = np.sort(min_times[non_zero])
    #min_mintimes = sorted_mintimes[0]
    min_mintimes = sorted_mintimes[2] #1496.9781737661062 #    # can change to manually increase livetime and reduce number of computations (don't change elsewhere as this should follow through!)
    logging.info(f'minimum time is {min_mintimes}')
    logging.info(f'Sorted livetimes: {sorted_mintimes}')        
    min_times[zero_ind] = 1
    angle_weight = min_mintimes/min_times
    angle_weight[zero_ind] = 0
    logging.info(f'Angle weights are {angle_weight}')    

    data['full_weight'] = np.zeros(len(theta))
    data['min_mintimes'] = min_mintimes
    data['angle_weight'] = angle_weight[:-1]
    logging.info(f'Count is {count}')

    for j in range(9):
        e_weight = np.array(data[j]['e_weight'])
        data[j]['full_weight'] = e_weight*(angle_weight[j])  # find final weights down to the energy bin
        data['full_weight'][np.where(theta_bins == j)] = (e_weight[e_bins[np.where(theta_bins == j)]])*(angle_weight[j])
    
    save_object(data, '/data/p-one/dghuman/simulation/MCEq_tables/muon_weights.pkl')

if __name__ == '__main__':
    main()
    
