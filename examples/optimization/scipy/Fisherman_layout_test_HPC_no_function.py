#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:44:08 2020

@author: dbensaso
"""

import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.optimization import YawOptimizationWindRoseParallel
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import os
import math 
import six
from mpi4py import MPI

#import copy
#import floris.tools.wind_rose as rose
#from mpi4py.futures import MPIPoolExecutor
#from itertools import repeat
#from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank > 0:
    print('DUNZO',rank)

else:
    print('Going forward with rank 0')
    # Instantiate the FLORIS object
    file_dir = os.path.dirname(os.path.abspath(__file__))
    fi = wfct.floris_interface.FlorisInterface(
        os.path.join(file_dir, '../../example_input.json')
    )
    
    # Function for plotting final tabular data as an image 
    def render_mpl_table(table_new, col_width=2.8, row_height=0.625, font_size=14,
                         header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                         bbox=[0, 0.1, 1, 0.8], header_columns=1,
                         ax=None, **kwargs):
        if ax is None:
            size = (np.array(table_new.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            ax.axis('off')
    
        mpl_table = ax.table(cellText=table_new.values, bbox=bbox, colLabels=table_new.columns, **kwargs)
    
        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)
    
        for k, cell in  six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
        return fig, ax
      
    # Define wind farm coordinates and layout
    wf_coordinate = [41.05, -70.2]
    
    # set min and max yaw offsets for optimization
    min_yaw = 0.0
    max_yaw = 25.0
    
    # Define minimum and maximum wind speed for optimizing power. 
    # Below minimum wind speed, assumes power is zero.
    # Above maximum_ws, assume optimal yaw offsets are 0 degrees
    minimum_ws = 3.0
    maximum_ws = 15.0
    
    # Instantiate the FLORIS object
    #file_dir = os.path.dirname(os.path.abspath(__file__))
    #fi = wfct.floris_interface.FlorisInterface(
    #    os.path.join(file_dir, '../../example_input.json')
    #)
    
    #Define Wake Model and input parameters 
    fi.floris.farm.wake.velocity_model.use_yaw_added_recovery = True
    fi.floris.farm.wake.deflection_model.use_secondary_steering = True
    
    fi.floris.farm.flow_field.wake.deflection_model.deflection_multiplier = 1.2
    fi.floris.farm.flow_field.wake.deflection_model.ka = 0.3
    
    
    # Set wind farm to N_row x t_row grid with constant spacing 
    kf= "Fishermans"
    
    D = 164
    N_row =12 
    T_row = 8
    Num_Turb = N_row*T_row
    spc_N = 11.2926829268 *(12/N_row)  #(1Nautical mile/164)
    spc_T= 11.2926829268 * (8/T_row) #(1Nautical mile/164)
    layout_x = []
    layout_y = []
    for i in range(N_row):
    	for k in range(T_row):
    		layout_x.append(i*spc_N*D)
    		layout_y.append(k*spc_T*D)
    
    
    new_layout_x= []
    new_layout_y=[]
    for i in layout_x:
        for j in layout_y:
            new_layout_x.append(i*math.cos(-45) -j *math.cos(-45))
            new_layout_y.append(i*math.cos(-45) +j *math.cos(-45))
    
    N_turb = len(layout_x)
    
    
    
    fi.reinitialize_flow_field(layout_array=(new_layout_x, new_layout_y),wind_direction=[270.0],wind_speed=[8.0])
    fi.calculate_wake()
    
    cp_8MW= [0,0,0,0,0.13,0.3,0.37,0.39,0.41,0.42,0.43,0.43,0.44,0.44,0.44,0.44,0.44,0.43,0.42,0.39,0.35,
             0.31,0.28,0.25,0.23,0.2,0.18,0.17,0.15,0.14,0.13,0.12,0.11,0.1,0.09,0.08,0.08,0.07,0.07,0.06,
             0.06,0.05,0.05,0.05,0.04,0.04,0.04,0.0]
    
    ct_8MW= [1.18,1.1,1.03,0.97,0.92,0.88,0.85,0.83,0.82,0.81,0.8,0.79,0.78,0.77,0.76,0.75,0.73,0.71,0.67,
             0.6,0.52,0.45,0.39,0.34,0.3,0.27,0.24,0.22,0.19,0.18,0.16,0.15,0.14,0.13,0.12,0.11,0.1,0.09,
             0.09,0.08,0.08,0.07,0.07,0.06,0.06,0.06,0.05,0.05]
    
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["power"] = cp_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["thrust"] = ct_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.rotor_diameter = 164
    #fi.floris.farm.flow_field.turbine_map.turbines.hub_height = 109
    
    #hub_h= 109
    #Diam= 164
    for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
            turbine.rotor_diameter = 109
            turbine.hub_height = 164
            cp_new = cp_8MW
            ct_new = ct_8MW
            turbine.power_thrust_table["power"] = cp_new
            turbine.power_thrust_table["thrust"] = ct_new
    
    unc_options={'std_wd': 4.95, 'std_yaw': 0.0,'pmf_res': 1.0, 'pdf_cutoff': 0.95}
    # ================================================================================
    print('Plotting the FLORIS flowfield...')
    # ================================================================================
    
    # Initialize the horizontal cut
    hor_plane = fi.get_hor_plane(
        height=fi.floris.farm.turbines[0].hub_height
    )
    
    # Plot and show
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
    ax.set_title('Baseline flow for U = 8 m/s, Wind Direction = 270$^\circ$')
    
    layout_name = str(kf) + "_layout.png"
    plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_layout/{}'.format(layout_name))
    # ================================================================================
    print('Importing wind rose data...')
    # ================================================================================
    
    # Create wind rose object and import wind rose dataframe using WIND Toolkit HSDS API.
    # Alternatively, load existing file with wind rose information.
    calculate_new_wind_rose = False
    
    wind_rose = rose.WindRose()
    
    if calculate_new_wind_rose:
    
    	wd_list = np.arange(0,360,5)
    	ws_list = np.arange(0,26,1)
    
    	df = wind_rose.import_from_wind_toolkit_hsds(wf_coordinate[0],
    	                                                    wf_coordinate[1],
    	                                                    ht = 100,
    	                                                    wd = wd_list,
    	                                                    ws = ws_list,
                                                            include_ti=True,
    	                                                    limit_month = None,
    	                                                    st_date = None,
    	                                                    en_date = None)
    
    else:
        file_name = str(kf) + "_Wind_Farm.p"
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle/{}'.format(file_name))
    
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["power"] = cp_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["thrust"] = ct_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.rotor_diameter = 164
    #f
    # plot wind rose
    wind_rose.plot_wind_rose()
    windrose_name = str(kf) + "_Wind_rose.png"
    plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose/{}'.format(windrose_name))
    #Plot ti rose
    wind_rose.plot_wind_rose_ti()
    ti_rose = str(kf) + "_ti_rose.png"
    plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_rose/{}'.format(ti_rose))
    # Plot ti vs. ws 
    wind_rose.plot_ti_ws()
    ti_ws = str(kf) + "_ti_ws.png"
    plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_ws/{}'.format(ti_ws))
    
    # =============================================================================
    print('Finding optimal yaw angles in FLORIS...')
    # =============================================================================
    
    #Three cases: 
    #1) Unc_and_base : Run optimization for both robust and non-robust
    #2) Just_Unc : Run for just robost 
    #3) Just_Base : Run for just non-robust
    
    Optimization_case= "Just_Base"
    
    if Optimization_case == "Unc_and_base":
    
    
        yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws,df.ti,
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
        yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws, df.ti,
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
        case_name = 'Example '+str(kf)+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_without_unc = str(kf) + "_report_without_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_without_unc))
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
        case_name_1 = 'Example '+str(kf)+ ' Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name_1, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_with_unc = str(kf) + "_report_with_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_with_unc))
        plt.show()
    
        #Save final data as a pickle (without_unc)
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(Num_Turb), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + str(kf) + "_without_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle))
        
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(data.iloc[0]['AEP_No_Wake']),3), round(float(data.iloc[0]['AEP_Baseline']),3), round(float(data.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data.iloc[0]['%_Baseline']),3), round(float(data.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(data.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        table_image = "Table_Image_" + str(kf)+ "_without_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
        #Save final data as a pickle (with unc)
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(Num_Turb), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + str(kf) + "_with_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle))
        
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(data.iloc[0]['AEP_No_Wake']),3), round(float(data.iloc[0]['AEP_Baseline']),3), round(float(data.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data.iloc[0]['%_Baseline']),3), round(float(data.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(data.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        table_image = "Table_Image_" + str(kf)+ "_with_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
        
    elif Optimization_case == "Just_Unc":
        
        yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws,df.ti,
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
        case_name = 'Example '+str(kf)+ ' Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_with_unc = str(kf) + "_report_with_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_with_unc))
        
        plt.show()
    
        #Save final data as a pickle 
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(Num_Turb), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + str(kf) + "_with_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle))
        
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(data.iloc[0]['AEP_No_Wake']),3), round(float(data.iloc[0]['AEP_Baseline']),3), round(float(data.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data.iloc[0]['%_Baseline']),3), round(float(data.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(data.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        table_image = "Table_Image_" + str(kf)+ "_with_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
    
    elif Optimization_case == "Just_Base":
        
        # Instantiate the Optimization object FOR NOW ASSUME TI WORKS
        yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws, df.ti,
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
        case_name = 'Example '+str(kf)+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_without_unc = str(kf) + "_report_without_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_without_unc))
        
        plt.show()
    
        #Save final data as a pickle 
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(Num_Turb), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + str(kf) + "_without_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle))
        
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(data.iloc[0]['AEP_No_Wake']),3), round(float(data.iloc[0]['AEP_Baseline']),3), round(float(data.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data.iloc[0]['%_Baseline']),3), round(float(data.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(data.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        table_image = "Table_Image_" + str(kf)+ "_without_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
    
    else: 
        raise SystemExit("None Valid Optimization Method Chosen")
    
    
    
