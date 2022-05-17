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

def find_closest(arr, val):
    index = np.argmin(np.abs(np.array(arr)-val))
    return index
        
def main():
    # Load files using pickle
    with open('/data/p-one/dghuman/simulation/MCEq_tables/new_muon_dict.pkl', 'rb') as input1:
        e_dict = pickle.load(input1)

    with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_dict_0006_triggered.pkl', 'rb') as input2:
        muon_dict = pickle.load(input2)

    with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_weights.pkl', 'rb') as input3:
        muon_weights = pickle.load(input3)

    with open('/data/p-one/dghuman/simulation/MCEq_tables/rate_muon_table_dtheta_10.pkl', 'rb') as input4:
        rate = pickle.load(input4)
        

    e_bins = np.zeros(len(muon_dict['theta']))
    e_grid = np.array(e_dict['e_grid'])

    weights = np.zeros(len(muon_dict['theta']))

    theta_bins = np.zeros(len(weights))
    bin_thetas = [x for x in rate.keys() if not isinstance(x, str)]

    data = {}

    for index, theta in enumerate(muon_dict['theta']):
        theta_bins[index] = int(find_closest(bin_thetas, theta))
        e_bins[index] = int(find_closest(e_grid, muon_dict['Energy'][index]))
        weights[index] = muon_weights[theta_bins[index]]['full_weight'][int(e_bins[index])]
        

    data['theta_bins'] = theta_bins
    data['e_bins'] = e_bins
    data['weights'] = weights
    
    save_object(data, '/data/p-one/dghuman/simulation/MCEq_tables/triggered_muon_weights.pkl')

if __name__ == '__main__':
    main()
    
