#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:58:20 2023

code to turn the FAFB skeletons into neurons

@author: jordankalai
"""

#@##Install dependencies
#STEP ONE: Import packages needed for code to run
import pymaid 
import matplotlib.pyplot as plt
#import pandas as pd
#from PIL import Image
import numpy 
import csv
import sys
from cloudvolume import CloudVolume
#from tifffile import imwrite 
import navis


output_folder_name = '/Users/jordankalai/Documents/Barnhart_Lab/MCFO/HS/SWC_skeletons/FAFB/'

cell_ID = ['Putative HSE 1088755 GA', 'Putative HSE 827035 ECS ZSG', 'Putative HSN 1088633 GA', 'Putative HSN 830794 ZSG', 'Putative HSS 3796481 PG', 'Putative HSS 408972 SJH ZSG']
cells = ['HSE_Left', 'HSE_Right', 'HSN_Left', 'HSN_Right', 'HSS_Left', 'HSS_Right']

"""
Sampling HSE Right for the figure because of asymmetry
"""



"""
STEP FOUR: Open Catmaid and select neuron of interest
"""
#initiate catmaid instance
rm = pymaid.CatmaidInstance('https://garden.catmaid.org/#',
                            api_token='26cf58fabcaa16cfeb6ce627d3b6276a05b418fa',
                            http_user='fly', 
                            http_password='superfly')

#Retreive Neuron of Interest 
for name, ID in zip(cells, cell_ID):
    n=pymaid.get_neuron(ID)
    backbone = navis.prune_by_strahler(n, to_prune=-1, inplace=False)

    pruneawaytwigs = navis.prune_twigs(n,size=5,recursive=float('inf'),inplace=False)

    despiked = navis.despike_skeleton(n) #no longer despike_neuron()

    n_pr = navis.prune_twigs(despiked, size= "5 microns", recursive=float('inf') , inplace=False) #recursive=float('inf')
    #take the neuron and turn it into an SWC (do not be concerned with seperating axon and dendrite right now)
    navis.write_swc(n_pr, output_folder_name+name+'.swc')