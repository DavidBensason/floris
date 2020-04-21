#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:34:44 2020

@author: dbensaso
"""
import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.yaw_wind_rose_parallel import YawOptimizationWindRoseParallel
import floris.tools.wind_rose as rose
import WakeSteering_US.cp_for_any_turb as cturb
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import os
import math
import six
from math import sqrt, floor
#import copy
#import floris.tools.wind_rose as rose
#from mpi4py.futures import MPIPoolExecutor
#from itertools import repeat
#from mpi4py import MPI

if __name__ == '__main__':

    # Instantiate the FLORIS object
    file_dir = os.path.dirname(os.path.abspath(__file__))
    fi = wfct.floris_interface.FlorisInterface(
        os.path.join(file_dir, '../../../example_input.json')
    )
    fi.reinitialize_flow_field(turbulence_intensity=[0.06])
    fi.set_gch(enable=False) ##Normal gauss model
    def isPerfect(N):
        """Function to check if a number is perfect square or not
    
        taken from:
        https://www.geeksforgeeks.org/closest-perfect-square-and-its-distance/
        by sahishelangia
        """
        if (sqrt(N) - floor(sqrt(N)) != 0):
            return False
        return True

    #Patrick layout code
    def getClosestPerfectSquare(N):
        """Function to find the closest perfect square taking minimum steps to
            reach from a number
    
        taken from:
        https://www.geeksforgeeks.org/closest-perfect-square-and-its-distance/
        by sahishelangia
        """
        if (isPerfect(N)):
            distance = 0
            return N, distance
    
        # Variables to store first perfect square number above and below N
        aboveN = -1
        belowN = -1
        n1 = 0
    
        # Finding first perfect square number greater than N
        n1 = N + 1
        while (True):
            if (isPerfect(n1)):
                aboveN = n1
                break
            else:
                n1 += 1
    
        # Finding first perfect square number less than N
        n1 = N - 1
        while (True):
            if (isPerfect(n1)):
                belowN = n1
                break
            else:
                n1 -= 1
    
        # Variables to store the differences
        diff1 = aboveN - N
        diff2 = N - belowN
    
        if (diff1 > diff2):
            return belowN, -diff2
        else:
            return aboveN, diff1
    
    
    def make_grid_layout(n_wt, D, grid_spc):
        """Make a grid layout (close as possible to a square grid)
    
        Inputs:
        -------
            n_wt : float
                Number of wind turbines in the plant
            D : float (or might want array_like if diff wt models are used)
                Wind turbine rotor diameter(s) in meters
            grid_spc : float
                Spacing between rows and columns in number of rotor diams D
            plant_cap_MW : float
                Total wind plant capacity in MW
    
        Returns:
        --------
            layout_x : array_like
                X positions of the wind turbines in the plant
            layout_y : array_like
                Y positions of the wind turbines in the plant
        """
    
        # Initialize layout variables
        layout_x = []
        layout_y = []
    
        # Find the closest square root
        close_square, dist = getClosestPerfectSquare(n_wt)
        close_square=int(close_square)
        dist=int(dist)
        side_length = int(sqrt(close_square))
    
        # Build a square grid
        for i in range(side_length):
            for k in range(side_length):
                layout_x.append(i*grid_spc*D)
                layout_y.append(k*grid_spc*D)
    
        # Check dist and determine what to do
        if dist == 0:
            # do nothing
            pass
        elif dist > 0:
            # square>n_wt : remove locations
            del(layout_x[close_square-dist:close_square])
            del(layout_y[close_square-dist:close_square])
        else:
            # square < n_w_t : add a partial row
            for i in range(abs(dist)):
                layout_x.append(sqrt(close_square)*grid_spc*D)
                layout_y.append(i*grid_spc*D)
    
        return layout_x, layout_y

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
      
    #Rated power of the turbine
    #turb_cap = [15]
    
    # Number of turbines that can fit in 1000MW farm
    #n_wt = [math.ceil(1000/turb_cap[0])]
    
    # D of turbine 
    #d_lst = [240]
    
    # Spacing we want to itterate through (5,1nm,10)
    #grid_spc = [5,(1852/d_lst[0]),10] #Will itterate through this 
    
    # Define wind farm coordinates and layout
    wf_coordinate = [40.99776, -124.671]
    #Rated power of the turbine
    turb_cap = [15]
    # Number of turbines that can fit in 1000MW farm
    n_wt = [math.ceil(1000/turb_cap[0])]
    # D of turbine 
    d_lst = [240]
    # Spacing we want to itterate through (5,1nm,10)
    grid_spc = [1852/d_lst[0]] #Will itterate through this 
    
    min_yaw = -25.0
    max_yaw = 25.0
    # Define minimum and maximum wind speed for optimizing power. 
    # Below minimum wind speed, assumes power is zero.
    # Above maximum_ws, assume optimal yaw offsets are 0 degrees
    minimum_ws = 3.0
    maximum_ws = 15.0
    data = pd.DataFrame([])
    data1 = pd.DataFrame([])
    #group = "Pat_duff_comp"
    N = n_wt[0]
    D= d_lst[0]
    spc = grid_spc[0]
    
    for ii in range(len(n_wt)):
        x, y = make_grid_layout(n_wt=N, D=D, grid_spc=spc)
        #xmax = max(x)
        #ymax = max(y)
        #plt.plot(x,y,'bo')
        #plt.show()

    # Set wind farm to N_row x t_row grid with constant spacing 
    kf= "Patrick_test_spc=" + str(round(spc,2))
    N_turb = len(x)
    
    fi.reinitialize_flow_field(layout_array=(x, y),wind_direction=[270.0],wind_speed=[8.0])
    fi.calculate_wake()
    
    cp_15MW = [0,0,0.100335552,0.272197208,0.359305118,0.409156685,0.46,0.441819614,0.475270298,
               0.484447627,0.488866488,0.489224161,0.489263048,0.48928802,0.489293116,0.489304304,
               0.489319143,0.489378059,0.435324299,0.380965216,0.335295525,0.3,0.263718683,0.235,
               0.211136579,0.19,0.171651956,0.15,0.14,0.12,0.108080744,0.1,0.092,0.085,0.078,0.072394937,
               0.065,0.06,0.058,0.051,0.050842443,0.045,0.043,0.041,0.04,0.037062292,0.037]

    ct_15MW= [0.81,0.81,0.819748943,0.801112031,0.808268424,0.821910918,0.8,0.823265981,0.834932456,
              0.829011103,0.811200233,0.805469658,0.804571567,0.803949121,0.803904895,0.803708734,
              0.80345211,0.801777393,0.607277698,0.501379105,0.425965654,0.37,0.32116631,0.28,0.2511023,
              0.22,0.201415182,0.18,0.161,0.134,0.125653944,0.11,0.101,0.1,0.09,0.08506697,0.075,
              0.072,0.068,0.065,0.061026446,0.056,0.05,0.049,0.045,0.045814967,0.045]

   
    for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
            turbine.rotor_diameter = d_lst[0]
            turbine.power_thrust_table["power"] = cp_15MW #Add new lists here
            turbine.power_thrust_table["thrust"] = ct_15MW
    for count, coord in enumerate(fi.floris.farm.flow_field.turbine_map.coords):
        coord.x3 = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height
    fi.floris.farm.flow_field.specified_wind_height = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height

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
    
    layout_name = str(kf) + "_layout_.png"
    plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Ti_Study/Saved_Fig/{}'.format(layout_name))
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
        #file_name = str(zf) + "_Wind_Farm.p"
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle/Patrick_Duffy_no_ti_Wind_Farm.p')
    
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["power"] = cp_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.power_thrust_table["thrust"] = ct_8MW
    #fi.floris.farm.flow_field.turbine_map.turbines.rotor_diameter = 164
    #f
    # plot wind rose
    #wind_rose.plot_wind_rose()
    #windrose_name = str(zf) + "_Wind_rose.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose/{}'.format(windrose_name))
    
    #Plot ti rose
    #wind_rose.plot_wind_rose_ti()
    #ti_rose = str(zf) + "_ti_rose.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_rose/{}'.format(ti_rose))
    # Plot ti vs. ws 
    #wind_rose.plot_ti_ws()
    #ti_ws = str(zf) + "_ti_ws.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_ws/{}'.format(ti_ws))
    
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
        
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(N_turb), 'Turbine_D':int(D),'Turb_spc_D': float(spc), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_D_Group_" + str(kf) + "_without_unc"
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
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
        #Save final data as a pickle (with unc)
        
        data1 = data1.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(N_turb),'Turbine_D':int(D),'Turb_spc_D': float(spc),'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle_1 = "Pickle_table_D_Group_" + str(kf) + "_with_unc"
        data1.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle_1))
        
        # Save final data as an image 
        farm_data_1 = [('AEP(GWh)',round(float(data1.iloc[0]['AEP_No_Wake']),3), round(float(data1.iloc[0]['AEP_Baseline']),3), round(float(data1.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data1.iloc[0]['%_Baseline']),3), round(float(data1.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data1.iloc[0]['Wk_Loss_Baseline']),3), round(float(data1.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data1.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data1.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new_1= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new_1)
        
        table_image_1 = "Table_Image_" + str(kf)+ "_with_unc"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image_1))
        
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
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(N_turb),'Turbine_D':int(D),'Turb_spc_D': float(spc), 'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_D_Group_" + str(kf) + "_with_unc"
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
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
    
    elif Optimization_case == "Just_Base":
        
        include_obv_ti = False
        if include_obv_ti:
            # Instantiate the Optimization object FOR NOW ASSUME TI WORKS
            yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws, df.ti,
                                           minimum_yaw_angle=min_yaw, 
                                           maximum_yaw_angle=max_yaw,
                                           minimum_ws=minimum_ws,
                                           maximum_ws=maximum_ws)
        else: 
            yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws,
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
        report_farm_without_unc = str(kf) +"_report_without_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Ti_Study/Saved_Fig/{}'.format(report_farm_without_unc))
        
        plt.show()
    
        #Save final data as a pickle 
        
        data = data.append(pd.DataFrame({'Farm Name': str(kf), '#Turbine': int(N_turb),'Turbine_D':int(D),'Turb_spc_D': round(spc,2),'Farm_lat':wf_coordinate[0], 'Farm_lon': wf_coordinate[1], 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + str(kf) + "_without_unc"
        data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Ti_Study/Saved_Fig/{}'.format(table_pickle))
        
        tabular = (data.loc[data['Turbine_spc_D'] == round(spc,2)])
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(tabular.iloc[0]['AEP_No_Wake']),3), round(float(tabular.iloc[0]['AEP_Baseline']),3), round(float(tabular.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(tabular.iloc[0]['%_Baseline']),3), round(float(tabular.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(tabular.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(tabular.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(tabular.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        table_image = "Table_Image_" + str(kf)+ "_without_unc"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Ti_Study/Saved_Fig/{}.png'.format(table_image))
    
    
    else: 
        raise SystemExit("None Valid Optimization Method Chosen")

# -*- coding: utf-8 -*-

