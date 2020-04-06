
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
from floris.tools.optimization.scipy.optimization import YawOptimizationWindRoseParallel
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import WakeSteering_US.namingfarm as nf
import WakeSteering_US.cp_for_any_turb as cturb
import pdb
import os
import six
#from mpi4py import MPI

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#if rank > 0:
#    print('DUNZO',rank)
#else:
#    print('Going forward with rank 0')
if __name__ == '__main__':
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
            
    # Instantiate the FLORIS object
    #file_dir = os.path.dirname(os.path.abspath(__file__))
    #fi = wfct.floris_interface.FlorisInterface(
    #    os.path.join(file_dir, '../../example_input.json')
    #)
    
    #Define Wake Model and input parameters 
    #fi.floris.farm.wake.velocity_model.use_yaw_added_recovery = True
    #fi.floris.farm.wake.deflection_model.use_secondary_steering = True
    
    #fi.floris.farm.flow_field.wake.deflection_model.deflection_multiplier = 1.2
    #fi.floris.farm.flow_field.wake.deflection_model.ka = 0.3
    
    # Define wind farm coordinates and layout
    #farm_name = "Jericho Mountain"
    mf =pd.read_pickle('/home/dbensaso/code/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
    kf = (mf.loc[mf['p_name'] == "Frontier I"])
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
    tf= pd.read_pickle('/home/dbensaso/code/WakeSteering_US/Working_dir_WS_US/Wind_Cp_look_up_table')
    
    ## Enumerate so for each turbine 
    for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
            turbine.rotor_diameter = D.iloc[count]
            turbine.hub_height = hub_h.iloc[count]
            T_Area = (np.pi* (D.iloc[count]**2)) /4
            U_turb_rated= (2* P_r.iloc[count]*(10**3)/ (C_p_rated * 1.225* T_Area))**(1/3)
            U_turb_norm =  tf.iloc[:,0] / U_turb_rated
            cp_new = cturb.cp_for_any_turb(U_turb_norm,tf)
            ct_new = cturb.ct_for_any_turb(U_turb_norm,tf)
            turbine.power_thrust_table["power"] = cp_new
            turbine.power_thrust_table["thrust"] = ct_new
    
    # set min and max yaw offsets for optimization 
    min_yaw = -25.0
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
    
    #layout_name = str(kf['p_name'].iloc[0]) + "_layout.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_layout/{}'.format(layout_name))
    #ax.set_title('Baseline flow for U = 8 m/s, Wind Direction = 270$^\circ$')
    #str(fi.reinitialize_flow_field.wind_speed)
    #str(fi.reinitialize_flow_field.wind_direction) #took these out from axis label 
    # ================================================================================
    print('Importing wind rose data...')
    # ================================================================================
    
    # Create wind rose object and import wind rose dataframe using WIND Toolkit HSDS API.
    # Alternatively, load existing .csv fi.le with wind rose information.
    calculate_wind_rose = False
    
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
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle/{}'.format(file_name))
        
        #    file_name = str(kf['p_name'].iloc[0]) + "_Wind Farm.p"
        #    df = wind_rose.load(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\pickle_files\{}'.format(file_name))
    
    # plot wind rose and save plots
    #file_name = str(kf['p_name'].iloc[0]) + "_Wind Farm.p"
    #wind_rose.save(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\pickle_files\{}'.format(file_name))
    
    #Plot Wind Rose
    #wind_rose.plot_wind_rose()
    #windrose_name = str(kf['p_name'].iloc[0]) + "_Wind_rose.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose/{}'.format(windrose_name))
    
    #Plot Ti rose 
    #wind_rose.plot_wind_rose_ti() 
    #ti_rose = kf['p_name'].iloc[0] + "_ti_rose.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_rose/{}'.format(ti_rose))
    
    #Plot Ti Winspeed dist.
    #wind_rose.plot_ti_ws() ## ALSO NOT WORKINg
    #ti_ws = kf['p_name'].iloc[0] + "_ti_ws.png"
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/ti_ws/{}'.format(ti_ws))
    #ti_ws_name = str(kf['p_name'].iloc[0]) + "_ti_ws.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\ti_ws_plots_farm\{}'.format(ti_ws_name))
    
    #ti_wd_name = str(kf['p_name'].iloc[0]) + "_ti_wd.jpg"
    #plt.savefig(r'C:\Users\dbensaso\Documents\Code\WakeSteering_US\Working_dir_WS_US\Saved_fig_data\ti_wd_plots_farm\{}'.format(ti_wd_name))
    #wind_rose.ti_plot_wd(kf)
    #print("Pause")
    #input("PRESS ENTER TO CONTINUE.")
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
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_without_unc = kf['p_name'].iloc[0] + "_report_without_unc.png"
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
        case_name_1 = 'Example '+kf['p_name'].iloc[0]+ 'Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name_1, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_with_unc = kf['p_name'].iloc[0] + "_report_with_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_with_unc))
        plt.show()
        
        #Save final data as a pickle (without_unc)
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': kf['p_name'].iloc[0], '#Turbine': len(kf), 'Farm_lat':kf["ylat"].mean(), 'Farm_lon': kf["xlong"].mean(), 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + kf['p_name'].iloc[0] + "_without_unc"
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
        
        table_image = "Table_Image_" + kf['p_name'].iloc[0]+ "_without_unc"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
        #Save final data as a pickle (with unc)
        data1 = pd.DataFrame([])
        data1 = data1.append(pd.DataFrame({'Farm Name': kf['p_name'].iloc[0], '#Turbine': len(kf), 'Farm_lat':kf["ylat"].mean(), 'Farm_lon': kf["xlong"].mean(), 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle_1 = "Pickle_table_" + kf['p_name'].iloc[0] + "_with_unc"
        data1.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle_1))
        
        # Save final data as an image 
        farm_data_1 = [('AEP(GWh)',round(float(data1.iloc[0]['AEP_No_Wake']),3), round(float(data1.iloc[0]['AEP_Baseline']),3), round(float(data1.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data1.iloc[0]['%_Baseline']),3), round(float(data1.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data1.iloc[0]['Wk_Loss_Baseline']),3), round(float(data1.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data1.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data1.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new_1= pd.DataFrame(farm_data_1, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new_1)
        
        table_image_1 = "Table_Image_" + kf['p_name'].iloc[0]+ "_with_unc"
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
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm with UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        report_farm_with_unc = kf['p_name'].iloc[0] + "_report_with_unc.png"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_with_unc))
        plt.show()
        
        #Save final data as a pickle 
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': kf['p_name'].iloc[0], '#Turbine': len(kf), 'Farm_lat':kf["ylat"].mean(), 'Farm_lon': kf["xlong"].mean(), 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        table_pickle = "Pickle_table_" + kf['p_name'].iloc[0] + "_with_unc"
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
        
        table_image = "Table_Image_" + kf['p_name'].iloc[0]+ "_with_unc"
        plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
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
        case_name = 'Example '+kf['p_name'].iloc[0]+ ' Wind Farm without UNC'
        power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
        
        # Save farm report with designated name and path (this case  HPC)
        #report_farm_without_unc = kf['p_name'].iloc[0] + "_report_without_unc.png"
        #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/farm_report/{}'.format(report_farm_without_unc))
        #plt.show()
        
        #Save final data as a pickle 
        data = pd.DataFrame([])
        data = data.append(pd.DataFrame({'Farm Name': kf['p_name'].iloc[0], '#Turbine': len(kf), 'Farm_lat':kf["ylat"].mean(), 'Farm_lon': kf["xlong"].mean(), 'AEP_No_Wake': power_rose.total_no_wake, 
                                         'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                         '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                         'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                         'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                         index=[0]), ignore_index=True)
        #table_pickle = "Pickle_table_" + kf['p_name'].iloc[0] + "_without_unc"
        #data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_pickle/{}'.format(table_pickle))
        df_base_pickle = "Pickle_base_" + kf['p_name'].iloc[0] + "_without_unc"
        df_base.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Pat_Model_Comp/{}'.format(df_base_pickle))
        
        df_opt_pickle = "Pickle_opt_" + kf['p_name'].iloc[0] + "_without_unc"
        df_opt.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Pat_Model_Comp/{}'.format(df_opt_pickle))
        
        # Save final data as an image 
        farm_data = [('AEP(GWh)',round(float(data.iloc[0]['AEP_No_Wake']),3), round(float(data.iloc[0]['AEP_Baseline']),3), round(float(data.iloc[0]['AEP_Opt']),3)), 
                ('%', '--', round(float(data.iloc[0]['%_Baseline']),3), round(float(data.iloc[0]['%_Opt']),3)), 
                ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(data.iloc[0]['Wk_Loss_Opt']),3)),
                ('%AEP_Gain', '--', '--', round(float(data.iloc[0]['AEP_Gain_Opt']),3)), 
                ('Loss_Reduced', '--', '--', round(float(data.iloc[0]['Loss_Red_Opt']),3))]
    
        table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        fig, ax = render_mpl_table(table_new)
        
        #table_image = "Table_Image_" + kf['p_name'].iloc[0]+ "_without_unc"
        #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/tabular_data_image/{}.png'.format(table_image))
    
    
    else: 
        raise SystemExit("None Valid Optimization Method Chosen")
    
    


