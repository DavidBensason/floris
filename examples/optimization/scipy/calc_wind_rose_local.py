#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:05:24 2020

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

#Define Wake Model and input parameters 
fi.floris.farm.wake.velocity_model.use_yaw_added_recovery = True
fi.floris.farm.wake.deflection_model.use_secondary_steering = True

fi.floris.farm.flow_field.wake.deflection_model.deflection_multiplier = 1.2
fi.floris.farm.flow_field.wake.deflection_model.ka = 0.3

# Define wind farm coordinates and layout
#farm_name = "Jericho Mountain"
mf =pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
kf = (mf.loc[mf['p_name'] == "Wildcat Ranch"])
wf_coordinate = [kf["ylat"].mean(),kf["xlong"].mean()]

# Set wind farm to N_row x N_row grid with constant spacing 
# (2 x 2 grid, 5 D spacing)
D = fi.floris.farm.turbines[0].rotor_diameter
lat_y = kf['ylat'].values
long_x = kf['xlong'].values

layout_x, layout_y= nf.longlat_to_utm(lat_y, long_x)
#layout_x=layout_x1.tolist()
#layout_y=layout_y1.tolist()

N_turb = len(layout_x)

fi.reinitialize_flow_field(layout_array=(layout_x, layout_y), wind_direction=[270.0],wind_speed=[8.0])
fi.calculate_wake()

#Diameter and Rated power based on the wind farm
D = kf["t_rd"]
P_r = kf["t_cap"]
hub_h = kf["t_hh"]

C_p_rated = 0.43003137
C_t_rated = 0.70701647

#Normalized wind speed for any turbine
tf= pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_Cp_look_up_table')

## Enumerate so for each turbine 
for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
        turbine.rotor_diameter = D.iloc[count]
        turbine.hub_height = hub_h.iloc[count]
        T_Area = (np.pi* (D.iloc[count]**2)) /4
        U_turb_rated= (2* P_r.iloc[count]*(10**3)/ (C_p_rated * 1.225* T_Area))**(1/3)
        U_turb_norm =  tf.iloc[:,0] / U_turb_rated
        cp_new = cturb.cp_for_any_turb(U_turb_norm)
        ct_new = cturb.ct_for_any_turb(U_turb_norm)
        turbine.power_thrust_table["power"] = cp_new
        turbine.power_thrust_table["thrust"] = ct_new

# set min and max yaw offsets for optimization 
min_yaw = 0.0
max_yaw = 25.0

# Define minimum and maximum wind speed for optimizing power. 
# Below minimum wind speed, assumes power is zero.
minimum_ws = 3.0
maximum_ws = 15.0

#UNC OPTION 
unc_options={'std_wd': 4.95, 'std_yaw': 0.0,'pmf_res': 1.0, 'pdf_cutoff': 0.95}

# ================================================================================
print('Importing wind rose data...')
# ================================================================================

# Create wind rose object and import wind rose dataframe using WIND Toolkit HSDS API.
# Alternatively, load existing .csv fi.le with wind rose information.
calculate_wind_rose = True

wind_rose = rose.WindRose()

if calculate_wind_rose:

	wd_list = np.arange(0,360,5)
	ws_list = np.arange(0,26,1)

	df = wind_rose.import_from_wind_toolkit_hsds(wf_coordinate[0],
	                                                    wf_coordinate[1],
	                                                    ht = int(hub_h.iloc[0]), #should this be a float 
	                                                    wd = wd_list,
	                                                    ws = ws_list,
                                                        include_ti=True,
                                                        limit_month = None,
	                                                    st_date = None,
	                                                    en_date = None)

   
else:
    file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
    df = wind_rose.load(os.path.join(final_path, file_name))
    
# plot wind rose and save plots. File name is that of the first farm in the list
file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
wind_rose.save(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle/{}'.format(file_name))