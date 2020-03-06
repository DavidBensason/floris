#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:13:20 2020

@author: dbensaso
"""

# Copyright 2019 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

# See read the https://floris.readthedocs.io for documentation


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

def WakeSteering_function(farm_name="Jericho Mountain", 
                          calculate_wind_rose= True, 
                          use_yaw_added_recovery= False,
                          use_secondary_steering= False,
                          deflection_multiplier= 1.2,
                          ka= 0.3,
                          Optimization_case= "Just_Base"):
        
    # Instantiate the FLORIS object
    #file_dir = os.path.dirname(os.path.abspath(__file__))
    #fi = wfct.floris_interface.FlorisInterface(
    #    os.path.join(file_dir, '../../example_input.json')
    #)
    
    #Define Wake Model and input parameters 
    fi.floris.farm.wake.velocity_model.use_yaw_added_recovery = use_yaw_added_recovery
    fi.floris.farm.wake.deflection_model.use_secondary_steering = use_secondary_steering
    
    fi.floris.farm.flow_field.wake.deflection_model.deflection_multiplier = deflection_multiplier
    fi.floris.farm.flow_field.wake.deflection_model.ka = ka
    
    # Define wind farm coordinates and layout
    #farm_name = "Jericho Mountain"
    mf =pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
    kf = (mf.loc[mf['p_name'] == farm_name])
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
    print('Plotting the FLORIS flowfield...')
    # ================================================================================
    
    # Initialize the horizontal cut
    #hor_plane = wfct.cut_plane.HorPlane(
    #    fi.get_flow_data(),
    #    fi.floris.farm.turbines[0].hub_height
    #)
    
    # Initialize the horizontal cut
    hor_plane = fi.get_hor_plane(
        height=fi.floris.farm.turbines[0].hub_height
    )
    
    # Plot and show
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
    ax.set_title("Baseline flow for U= 8 m/s  ,  Wind Direction= 270 $^\circ$")
    #layout_name = str(kf['p_name'].iloc[0]) + "_layout.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\layout_farms\{}'.format(layout_name))
    #ax.set_title('Baseline flow for U = 8 m/s, Wind Direction = 270$^\circ$')
    #str(fi.reinitialize_flow_field.wind_speed)
    #str(fi.reinitialize_flow_field.wind_direction) #took these out from axis label 
    # ================================================================================
    print('Importing wind rose data...')
    # ================================================================================
    
    # Create wind rose object and import wind rose dataframe using WIND Toolkit HSDS API.
    # Alternatively, load existing .csv fi.le with wind rose information.
    #calculate_wind_rose = True
    
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
        df = wind_rose.load(os.path.join(file_dir, 'windtoolkit_geo_center_us.p'))
        
        #    file_name = str(kf['p_name'].iloc[0]) + "_Wind Farm.p"
        #    df = wind_rose.load(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\pickle_files\{}'.format(file_name))
    
    # plot wind rose and save plots
    #file_name = str(kf['p_name'].iloc[0]) + "_Wind Farm.p"
    #wind_rose.save(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\pickle_files\{}'.format(file_name))
    
    #Plot Wind Rose
    wind_rose.plot_wind_rose()
    
    
    #windrose_name = str(kf['p_name'].iloc[0]) + "_Wind_rose.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\windrose_farms\{}'.format(windrose_name))
    
    wind_rose.plot_wind_rose_ti() ##NOT WORKING CHECK THIS
    
    wind_rose.plot_ti_ws() ## ALSO NOT WORKINg
    
    #ti_ws_name = str(kf['p_name'].iloc[0]) + "_ti_ws.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\ti_ws_plots_farm\{}'.format(ti_ws_name))
    
    #ti_wd_name = str(kf['p_name'].iloc[0]) + "_ti_wd.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\ti_wd_plots_farm\{}'.format(ti_wd_name))
    #wind_rose.ti_plot_wd(kf)
    print("Pause")
    input("PRESS ENTER TO CONTINUE.")
    # =============================================================================
    print('Finding optimal yaw angles in FLORIS...')
    # =============================================================================
    
    #Three cases: 
    #1) Unc_and_base : Run optimization for both robust and non-robust
    #2) Just_Unc : Run for just robost 
    #3) Just_Base : Run for just non-robust
    
    #Optimization_case= "Unc_and_base"
    
    if Optimization_case == "Unc_and_base":
    
    
        yaw_opt = YawOptimizationWindRose(fi, df.wd, df.ws,df.ti,
                                   minimum_yaw_angle=min_yaw,
                                   maximum_yaw_angle=max_yaw,
                                   minimum_ws=minimum_ws,
                                   maximum_ws=maximum_ws,
                                   include_unc=True,
                                   unc_options=unc_options)
        
        # Determine baseline power with and without wakes
        df_base_unc = yaw_opt.calc_baseline_power()
        # Perform optimization
        df_opt_unc = yaw_opt.optimize()
        
        
        # Instantiate the Optimization object FOR NOW ASSUME TI WORKS
        yaw_opt = YawOptimizationWindRose(fi, df.wd, df.ws, df.ti,
                                       minimum_yaw_angle=min_yaw, 
                                       maximum_yaw_angle=max_yaw,
                                       minimum_ws=minimum_ws,
                                       maximum_ws=maximum_ws)
        
        # Determine baseline power with and without wakes
        df_base = yaw_opt.calc_baseline_power()
        # Perform optimization
        df_opt = yaw_opt.optimize()
        
        ##For cases WITHOUT UNC.
        # combine wind farm-level power into one dataframe
        df_power = pd.DataFrame({'ws':df.ws,'wd':df.wd, \
            'freq_val':df.freq_val,'power_no_wake':df_base.power_no_wake, \
            'power_baseline':df_base.power_baseline,'power_opt':df_opt.power_opt})
        # initialize power rose
        df_yaw = pd.DataFrame([list(row) for row in df_opt['yaw_angles']],columns=[str(i) for i in range(1,N_turb+1)])
        df_yaw['ws'] = df.ws
        df_yaw['wd'] = df.wd
        df_turbine_power_no_wake = pd.DataFrame([list(row) for row in df_base['turbine_power_no_wake']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_no_wake['ws'] = df.ws
        df_turbine_power_no_wake['wd'] = df.wd
        df_turbine_power_baseline = pd.DataFrame([list(row) for row in df_base['turbine_power_baseline']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_baseline['ws'] = df.ws
        df_turbine_power_baseline['wd'] = df.wd
        df_turbine_power_opt = pd.DataFrame([list(row) for row in df_opt['turbine_power_opt']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_opt['ws'] = df.ws
        df_turbine_power_opt['wd'] = df.wd
        
        # Summarize using the power rose module
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        plt.show()
        
        # Summarize AEP gains WITH uncertainty
        # combine wind farm-level power into one dataframe
        df_power = pd.DataFrame({'ws':df.ws,'wd':df.wd, \
            'freq_val':df.freq_val,'power_no_wake':df_base_unc.power_no_wake, \
            'power_baseline':df_base_unc.power_baseline,'power_opt':df_opt_unc.power_opt})
        
        # initialize power rose
        df_yaw = pd.DataFrame([list(row) for row in df_opt_unc['yaw_angles']],columns=[str(i) for i in range(1,N_turb+1)])
        df_yaw['ws'] = df.ws
        df_yaw['wd'] = df.wd
        df_turbine_power_no_wake = pd.DataFrame([list(row) for row in df_base_unc['turbine_power_no_wake']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_no_wake['ws'] = df.ws
        df_turbine_power_no_wake['wd'] = df.wd
        df_turbine_power_baseline = pd.DataFrame([list(row) for row in df_base_unc['turbine_power_baseline']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_baseline['ws'] = df.ws
        df_turbine_power_baseline['wd'] = df.wd
        df_turbine_power_opt = pd.DataFrame([list(row) for row in df_opt_unc['turbine_power_opt']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_opt['ws'] = df.ws
        df_turbine_power_opt['wd'] = df.wd
        # Summarize using the power rose module
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        plt.show()
    
    elif Optimization_case == "Just_Unc":
        
        yaw_opt = YawOptimizationWindRose(fi, df.wd, df.ws,df.ti,
                                       minimum_yaw_angle=min_yaw,
                                       maximum_yaw_angle=max_yaw,
                                       minimum_ws=minimum_ws,
                                       maximum_ws=maximum_ws,
                                       include_unc=True,
                                       unc_options=unc_options)
        
        # Determine baseline power with and without wakes
        df_base_unc = yaw_opt.calc_baseline_power()
        # Perform optimization
        df_opt_unc = yaw_opt.optimize()
        
        # Summarize AEP gains WITH uncertainty
        # combine wind farm-level power into one dataframe
        df_power = pd.DataFrame({'ws':df.ws,'wd':df.wd, \
            'freq_val':df.freq_val,'power_no_wake':df_base_unc.power_no_wake, \
            'power_baseline':df_base_unc.power_baseline,'power_opt':df_opt_unc.power_opt})
        
        # initialize power rose
        df_yaw = pd.DataFrame([list(row) for row in df_opt_unc['yaw_angles']],columns=[str(i) for i in range(1,N_turb+1)])
        df_yaw['ws'] = df.ws
        df_yaw['wd'] = df.wd
        df_turbine_power_no_wake = pd.DataFrame([list(row) for row in df_base_unc['turbine_power_no_wake']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_no_wake['ws'] = df.ws
        df_turbine_power_no_wake['wd'] = df.wd
        df_turbine_power_baseline = pd.DataFrame([list(row) for row in df_base_unc['turbine_power_baseline']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_baseline['ws'] = df.ws
        df_turbine_power_baseline['wd'] = df.wd
        df_turbine_power_opt = pd.DataFrame([list(row) for row in df_opt_unc['turbine_power_opt']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_opt['ws'] = df.ws
        df_turbine_power_opt['wd'] = df.wd
        # Summarize using the power rose module
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        plt.show()
    
    elif Optimization_case == "Just_Base":
        
        # Instantiate the Optimization object FOR NOW ASSUME TI WORKS
        yaw_opt = YawOptimizationWindRose(fi, df.wd, df.ws, df.ti,
                                       minimum_yaw_angle=min_yaw, 
                                       maximum_yaw_angle=max_yaw,
                                       minimum_ws=minimum_ws,
                                       maximum_ws=maximum_ws)
    
        # Determine baseline power with and without wakes
        df_base = yaw_opt.calc_baseline_power()
        # Perform optimization
        df_opt = yaw_opt.optimize()
    
        ##For cases WITHOUT UNC.
        # combine wind farm-level power into one dataframe
        df_power = pd.DataFrame({'ws':df.ws,'wd':df.wd, \
            'freq_val':df.freq_val,'power_no_wake':df_base.power_no_wake, \
            'power_baseline':df_base.power_baseline,'power_opt':df_opt.power_opt})
        # initialize power rose
        df_yaw = pd.DataFrame([list(row) for row in df_opt['yaw_angles']],columns=[str(i) for i in range(1,N_turb+1)])
        df_yaw['ws'] = df.ws
        df_yaw['wd'] = df.wd
        df_turbine_power_no_wake = pd.DataFrame([list(row) for row in df_base['turbine_power_no_wake']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_no_wake['ws'] = df.ws
        df_turbine_power_no_wake['wd'] = df.wd
        df_turbine_power_baseline = pd.DataFrame([list(row) for row in df_base['turbine_power_baseline']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_baseline['ws'] = df.ws
        df_turbine_power_baseline['wd'] = df.wd
        df_turbine_power_opt = pd.DataFrame([list(row) for row in df_opt['turbine_power_opt']],columns=[str(i) for i in range(1,N_turb+1)])
        df_turbine_power_opt['ws'] = df.ws
        df_turbine_power_opt['wd'] = df.wd
        
        # Summarize using the power rose module
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        plt.show()
    
    else: 
        raise SystemExit("None Valid Optimization Method Chosen")
    