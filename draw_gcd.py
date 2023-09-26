#!/usr/bin/env python                                                                                                

# Import some useful ICECUBE modules                                                                                  
from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
from icecube.dataclasses import I3Particle 
import numpy as np              
import matplotlib
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import pickle
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import sys
import math as m
import os, sys

import argparse

def fill_OMs(_geo,_dom):
    for omkey in _geo.keys():
        pos = _geo[omkey].position
        _dom[0].append(pos.theta*180/np.pi)
        _dom[1].append(pos.phi*180/np.pi)
        _dom[2].append(pos.x)
        _dom[3].append(pos.y)
        _dom[4].append(pos.z)

def plot_OMs(ax, _dom):
    ax.scatter(_dom[2], _dom[3], _dom[4], marker='o', color='blue')
    #plt.grid()
    ax.set_xlabel('x (m)')
    ax.set_ylabel("y (m)")
    ax.set_zlabel('z (m)')
    #ax.title(r'Location of DOMs')
    #plt.savefig('/home/users/ghuman/simAnalysis/output/plots/llhPDF/DOM_Locations_3d.eps')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gcd",default='/data/p-one/sim/muon/sim0006/jobfiles/PONE_10String.i3.gz', help="Read in GCD file")
    parser.add_argument("-o", "--outdir", type=str, default="/data/p-one/dghuman/simulation/MCEq_tables/")
    args = parser.parse_args()

    gcd = dataio.I3File(args.gcd)
    frame = gcd.pop_frame()
#    geo = frame["I3Geometry"]
    geo = frame["I3OMGeoMap"]

    dom = [[], [], [], [], []]
    muon_vert = []
    muon_dir = []

    fill_OMs(geo,dom)
    fig = plt.figure(figsize=(10,10))
    fig.clf()
    ax = plt.axes(projection="3d")
    plot_OMs(ax,dom)
    ax.view_init(30, 45) #Top Down view (90,90)
    plt.savefig(args.outdir + 'GCD_Plot.png', dpi=300)

if __name__ == '__main__':
    main()

'''    
    dom_thetas = []
    dom_phis = []
    dom_x = []
    dom_y = []
    dom_z = []
'''
    
