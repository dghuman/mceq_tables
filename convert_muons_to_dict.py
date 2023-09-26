#!/usr/bin/env python                                                                                                

# Import some useful ICECUBE modules                                                                                  
from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame  
from icecube.dataclasses import I3Particle 
import numpy as np              
import matplotlib
import pickle
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import sys
import math as m
import os, sys

# Save file function using pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)

#directory = '/data/p-one/sim/muon/sim0006/eventfiles'
directory = '/data/p-one/sim/muon/sim0006/daq6_6hit'   #atm_daq       # Files for events that were triggered
filelist_aux = os.listdir(directory)
filelist = [x for x in filelist_aux if '.i3.zst' in x ] #and os.path.getsize(directory + '/' + x) > 2000000]

muon_dict = {}
muon_dict['Energy'] = []
muon_dict['theta'] = []
muon_dict['phi'] = []

for filel in filelist:
    infile = dataio.I3File(directory + '/' + filel)

    for frame in infile:
#        if not frame.Has('I3EventHeader'):
#            continue

#        MMCTrackList = frame['MMCTrackList']
        try:
#            Muon = MMCTrackList[0].GetI3Particle()
            mctree = frame['I3MCTree'] #frame['I3MCTree_preMuonProp'] #frame['MuonGeneratorI3MCTree']
        except:
            continue
        simplecount = len(frame['singleDOMTrigger_3PMT_2DOM'])
        if simplecount < 2:
            continue
        Muon = mctree[0]
        muon_dir = Muon.dir
        muon_theta = np.abs(muon_dir.theta - np.pi)*180./np.pi # converted to degrees and shifted to point to source
        muon_phi = (muon_dir.phi)*180./np.pi
        muon_E = (Muon.energy)
        muon_dict['theta'].append(muon_theta)
        muon_dict['phi'].append(muon_phi)
        muon_dict['Energy'].append(muon_E)

'''
        # Check if 6 unique DOMs are hit
        pmtsplit = frame['I3Photons_pmtsplit']
        elements = set()                        # Sets only preserve unique elements
        for OM in pmtsplit.keys():
            elements.add(str(OM.string) + '_' + str(OM.om))   # string guarentees the comparison will be consistent
        if len(elements) < 6:
            continue
'''
        
save_object(muon_dict, "/data/p-one/dghuman/simulation/MCEq_tables/muon_dict_0006_triggered_6string_6hit.pkl")
