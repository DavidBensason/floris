#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:38:16 2020

@author: dbensaso
"""
import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.optimization import YawOptimizationWindRose
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import WakeSteering_US.WakeSteering_US.namingfarm as nf
import WakeSteering_US.WakeSteering_US.cp_for_any_turb as cturb
import pdb
import os

# Instantiate the FLORIS object
file_dir = os.path.dirname(os.path.abspath(__file__))
fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, '../../example_input.json')
)


      
# Instantiate the FLORIS object
#file_dir = os.path.dirname(os.path.abspath(__file__))
#fi = wfct.floris_interface.FlorisInterface(
#    os.path.join(file_dir, '../../example_input.json')
#)

#Define Wake Model and input parameters 
fi.floris.farm.wake.velocity_model.use_yaw_added_recovery = False
fi.floris.farm.wake.deflection_model.use_secondary_steering = False

fi.floris.farm.flow_field.wake.deflection_model.deflection_multiplier = 1.2
fi.floris.farm.flow_field.wake.deflection_model.ka = 0.3

# Define wind farm coordinates and layout
#farm_name = "Jericho Mountain"
mf =pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
kf = (mf.loc[mf['p_name'] == "Apple Blossom"])
wf_coordinate = [kf["ylat"].mean(),kf["xlong"].mean()]

# Set wind farm to N_row x N_row grid with constant spacing 
# (2 x 2 grid, 5 D spacing)
D = fi.floris.farm.turbines[0].rotor_diameter
lat_y = kf['ylat'].values
long_x = kf['xlong'].values

layout_x, layout_y= nf.longlat_to_utm(lat_y, long_x)
#layout_x=layout_x1.tolist()
#layout_y=layout_y1.tolist()
coord = pd.DataFrame({'x':layout_x, 'y':layout_y})


plt.figure()
plt.gca().invert_yaxis()

x = np.linspace(-600,4000,5)
x = [-600,4000]
y0 = [-3000,-3000]
y1 = [-2000,-2000]
y2 = [-1000,-1000]
y3 = [0,0]
y4 = [1000,1000]
y5 = [2000,2000]
plt.plot(x,y0,x,y1,x,y2,x,y3,x,y4,x,y5)
plt.plot(layout_x,layout_y, 'r.', markersize=14)

plt.show()

group1 = pd.DataFrame([])
group2 = pd.DataFrame([])
group3 = pd.DataFrame([])
group4 = pd.DataFrame([])
group5 = pd.DataFrame([])
#for i in range(coord):
for i, row in coord.iterrows():
#Check if point lies in the line 

    x = layout_x[i]

    y = layout_y[i]
    c3 = 1* y + 0.3* x
    
    if y0[1]<= y < y1[1]: 
        group1= group1.append(pd.DataFrame({'x': x, 'y': y},index=[0]), ignore_index=True)
    if y1[1] <= y < y2[1]: 
        group2= group2.append(pd.DataFrame({'x': x, 'y': y},index=[0]), ignore_index=True)
    if y2[1]<= y < y3[1]: 
        group3= group3.append(pd.DataFrame({'x': x, 'y': y},index=[0]), ignore_index=True)
    if y3[1]<= y < y4[1]: 
        group4= group4.append(pd.DataFrame({'x': x, 'y': y},index=[0]), ignore_index=True)
    if y4[1] <= y < y5[1]: 
        group5= group5.append(pd.DataFrame({'x': x, 'y': y},index=[0]), ignore_index=True)
