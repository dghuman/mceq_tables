# /data/p-one/mceq_db_lext_dpm191_v131.h5

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os, sys

#import solver related modules
from MCEq.core import MCEqRun
import mceq_config as config
from MCEq.geometry.density_profiles import GeneralizedTarget
#import primary model choices
import crflux.models as pm

# Save file function using pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)


# configure
config.floatlen = "float64"
config.kernel_config = "numpy"
#config.e_min = 10.
config.enable_default_tracking = False
config.muon_helicity_dependence = False
config.dXmax = 1.
config.max_density = 0.001225

config.debug_level = 2


mceq_surface = MCEqRun(
    #provide the string of the interaction model
    interaction_model='SIBYLL2.3c',
    #primary cosmic ray flux model
    #support a tuple (primary model class (not instance!), arguments)
    primary_model=(pm.GlobalSplineFitBeta, None),
    # Zenith angle in degrees. 0=vertical, 90=horizontal
    theta_deg=0.0
)

# Create a homogeneous target material
target_water = GeneralizedTarget(len_target=1550e2) # 2km
density = 1.03
target_water.add_material(0.,density,'Water')
depths = np.arange(10., 1550.1 , 10)*100 #read out at steps of 10 meters

#mceqconfig.adv_set["exclude_from_mixing"] = [13, -13]
config.leading_process = 'auto'
config.max_density = density
mceq_water = MCEqRun(
    
    interaction_model='SIBYLL2.3c',
    
    primary_model = None,
    
    theta_deg = None,
    
    medium='water',
    
    particle_list=[(13,0),(-13,0)],
    
    density_model = target_water
    
)

#Power of energy to scale the flux (the results will be returned as E**mag * flux)
mag = 0

#obtain energy grid (fixed) of the solution for the x-axis of the plots
#e_grid = mceq_surface.e_grid

#Dictionary for results
flux = {}

#Step size of angles in degrees for tables
dtheta = 10

#Define a zenith angle, counted positively from vertical direction. Theta = 0. means vertical, theta = 90. horizontal
for theta in np.arange(0, 90, dtheta):
    #Start at surface
    #Set the zenith angle
    mceq_surface.set_theta_deg(theta)
    #Run the solver
    mceq_surface.solve()
    
    #Switch to water
    mceq_water._set_state_vector(*mceq_surface._get_state_vector(), only_available=True)
    #Set length of water
    target_water.set_length((2000/np.cos(theta*np.pi/180))*100)
    mceq_water._calculate_integration_path(int_grid=depths, grid_var='X', force=True)
    mceq_water.solve()

    # total means conventional + prompt
    mu_total = (mceq_water.get_solution('total_mu+', mag, integrate=False)
                + mceq_water.get_solution('total_mu-', mag, integrate=False))

    flux[theta] = mu_total

e_grid = mceq_water.e_grid
e_widths = mceq_water.e_widths
e_bins = mceq_water.e_bins
flux['e_grid'] = e_grid
flux['e_widths'] = e_widths
flux['e_bins'] = e_bins
flux['depths'] = depths
print("Writing to /data/p-one/dghuman/simulation/MCEq_tables/muon_table_dtheta_" + str(dtheta) + ".pkl")
save_object(flux, "/data/p-one/dghuman/simulation/MCEq_tables/muon_table_dtheta_" + str(dtheta) + ".pkl")
