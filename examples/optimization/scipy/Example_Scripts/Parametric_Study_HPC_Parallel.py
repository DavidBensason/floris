#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 08:23:38 2020

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

#Example for running parametric study on the HPC. For this example, the farm layout is a 7x7 grid, 
## with relative spacing of 7D, turbine specific power of 150W/m^2 and average wind speed of 9m/s

if __name__ == '__main__':      
     
    # Instantiate the FLORIS object
    file_dir = os.path.dirname(os.path.abspath(__file__))
    fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, '../../../example_input.json')
    )
    
    # Define opt_options
    opt_options = {'maxiter': 20, 'disp': False,'iprint': 1, 'ftol': 1e-6, 'eps': 0.01}

    data = pd.DataFrame([])
    # Define wind farm coordinates and layout
    wf_coordinate = [40.995293, -84.565796]
    # set min and max yaw offsets for optimization
    min_yaw = -25.0
    max_yaw = 25.0
    
    
    minimum_ws = 3.0
    maximum_ws = 15.0
    
    #data = pd.DataFrame([])
    data1 = pd.DataFrame([])
    group= "spc=5_g2" ##Change for with next run 
    
    # Specify desired average wind speed, relative spacing and turbine specific power
    scale_ws_avg = True
    if scale_ws_avg:
        des_avg_ws = 9.0
    
    rel_spc = 7
    
    D = 160
    SP = 150
  
    # Set wind farm to N_row x t_row grid with constant spacing 
    if scale_ws_avg:
        kf= "Rel_spc="+ str(rel_spc) +"_D=" + str(D) + '_avg_ws=' +str(des_avg_ws)\
        + '_SP=' + str(SP)
    else: 
        kf= "Rel_spc="+ str(rel_spc) +"_D=" + str(D) \
        + '_SP=' + str(SP)

    #Rated Power
    P_r = SP * ((D**2)/4) * 10**-3 * math.pi
    N_row =7 
    T_row = 7
    layout_x = []
    layout_y = []
    Num_Turb = N_row*T_row

    constant_area_layout = False
    #relative_original_spacing = False    #y_n and x_t are multiples of (1852/D)=11.29
    if constant_area_layout: ##THIS WORKS ATM 
        spc_N = rel_spc*1852  #(1Nautical mile/164)
        #spc_T= abs_spacing #(1Nautical mile/164)
       
    else: 
        spc_N = rel_spc
                
    for j in range(N_row):
        for k in range(T_row):
            layout_x.append(j*spc_N*D)
            layout_y.append(k*spc_N*D)

    
    N_turb = len(layout_x)
    
    
    # Reinitialize flow field and specify ambient TI. In this case, it is 8%. 
    fi.reinitialize_flow_field(layout_array=(layout_x, layout_y),wind_direction=[270.0],wind_speed=[8.0])
    fi.calculate_wake()
    fi.reinitialize_flow_field(turbulence_intensity=[0.08])
    # Using the NREL 5MW reference turbine     
    C_p_rated = 0.43003137
    C_t_rated = 0.70701647
    
    #Normalized wind speed for any turbine
    tf= pd.read_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/NREL_5MW_reference')

    
    for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
            turbine.rotor_diameter = int(D)
            turbine.hub_height = 95
            T_Area = (np.pi* (int(D)**2)) /4
            U_turb_rated= (2* P_r*(10**3)/ (C_p_rated * 1.225* T_Area))**(1/3)
            U_turb_norm =  tf.iloc[:,0] / U_turb_rated
            cp_new = cturb.cp_for_any_turb(U_turb_norm,U_turb_rated,T_Area,P_r,tf)
            ct_new = cturb.ct_for_any_turb(U_turb_norm,tf)
            turbine.power_thrust_table["power"] = cp_new
            turbine.power_thrust_table["thrust"] = ct_new
            turbine.change_turbine_parameters({})
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
    
    #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Parametric_Study/farm_layout/{}'.format(layout_name))
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
    	                                                    ht = 95,
    	                                                    wd = wd_list,
    	                                                    ws = ws_list,
                                                            include_ti=False,
    	                                                    limit_month = None,
    	                                                    st_date = None,
    	                                                    en_date = None)
    
    else:
        ## Load presaved parametric wind rose
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle/Onshore_case_Ohio.p')
    
    # Adjust wind speed column in wind rose to match the desired average wind speed   
    if scale_ws_avg: 
        ws_list = np.arange(0,26,1)
       # ws_list = np.arange(0,26,1)
        data_overall = pd.DataFrame([])
        for iii in ws_list: 
            dw = df.loc[df['ws'] == iii]
            freq_sum = dw['freq_val'].sum()
            product = iii * freq_sum
            data_overall = data_overall.append({'ws': iii, 'freq': freq_sum, 'product': product}, ignore_index=True)
        
        # Average ws in the farm 
        avg_ws = data_overall['product'].sum() / data_overall['freq'].sum()
        
        # Specify desired avg_ws
        scale_fac = des_avg_ws / avg_ws
        
        # Change the df ws column to match the des average 
        df['ws'] = df['ws'].apply(lambda x: x*scale_fac)
        
    
    # =============================================================================
    print('Finding optimal yaw angles in FLORIS...')
    # =============================================================================
    
    #Three cases: 
    #1) Just_Unc : Run for just robost 
    #2) Just_Base : Run for just non-robust
    
    Optimization_case= "Just_Base"
    
    
    
    
    if Optimization_case == "Just_Base":
        
        # Instantiate the Optimization object FOR NOW ASSUME TI WORKS
        yaw_opt = YawOptimizationWindRoseParallel(fi, df.wd, df.ws,
                                       minimum_yaw_angle=min_yaw, 
                                       maximum_yaw_angle=max_yaw,
                                       minimum_ws=minimum_ws,
                                       maximum_ws=maximum_ws, 
                                       opt_options= opt_options)
    
        # Determine baseline power with and without wakes
        df_base = yaw_opt.calc_baseline_power()
        # Perform optimization
        df_opt = yaw_opt.optimize()
        
        ## Save df_base and df_opt to pickle file 
        df_base_pickle = "Df_base_" + str(kf) + "_without_unc"
        #df_base.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Parametric_Study/df_base_pickle/{}'.format(df_base_pickle))
        
        df_opt_pickle = "Df_opt_" + str(kf) + "_without_unc"
        #df_opt.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Parametric_Study/df_opt_pickle/{}'.format(df_opt_pickle))
        
        
        # Summarize using the power rose module
        case_name = 'Example '+str(kf)+ ' Wind Farm without UNC'
        #power_rose = pr.PowerRose(case_name, df_power, df_turbine_power_no_wake, df_turbine_power_baseline,df_yaw, df_turbine_power_opt)
        power_rose = pr.PowerRose()
        power_rose.make_power_rose_from_user_data(
        	case_name,
        	df,
        	df_base['power_no_wake'],
        	df_base['power_baseline'],
            df_opt['power_opt']
        )
        
        fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(6.4, 6.5))
        power_rose.plot_by_direction(axarr)
        power_rose.report()
    
        # Save farm report with designated name and path (this case  HPC)
        report_farm_without_unc = str(kf) +"_report_without_unc.png"
        #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Parametric_Study/farm_report/{}'.format(report_farm_without_unc))
        
        #plt.show()
    
        #Save final data as a pickle 
        if scale_ws_avg:
            data = data.append(pd.DataFrame({'Farm Name': str(kf), 'Rated_Power': P_r, 'Specific_Power':SP,'#Turbine': int(Num_Turb), 'Turbine_D':int(D),
                                                 'Turb_spc_rel': float(rel_spc),'Avg_ws': float(des_avg_ws),'AEP_No_Wake': power_rose.total_no_wake, 
                                                 'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                                 '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                                 'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                                 'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                                 index=[0]), ignore_index=True)
        else: 
                    
            data = data.append(pd.DataFrame({'Farm Name': str(kf), 'Rated_Power': P_r, 'Specific_Power':SP,'#Turbine': int(Num_Turb), 'Turbine_D':int(D),
                                                 'Turb_spc_rel': float(rel_spc),'AEP_No_Wake': power_rose.total_no_wake, 
                                                 'AEP_Baseline': power_rose.total_baseline, 'AEP_Opt':power_rose.total_opt, 
                                                 '%_Baseline': 100.* power_rose.baseline_percent, '%_Opt': 100.* power_rose.opt_percent, 
                                                 'Wk_Loss_Baseline':100.* power_rose.baseline_wake_loss, 'Wk_Loss_Opt': 100.* power_rose.opt_wake_loss, 
                                                 'AEP_Gain_Opt': 100.* power_rose.percent_gain , 'Loss_Red_Opt':100.* power_rose.reduction_in_wake_loss}, 
                                                 index=[0]), ignore_index=True)

        table_pickle = "Pickle_table_" + str(kf) + "_without_unc"
        #data.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Parametric_Study/tabular_data_pickle/{}'.format(table_pickle))
        
        # Save final data as an image 
        #tabular = (data.loc[data['Turbine_D'] == D])
        # Save final data as an image 
        #farm_data = [('AEP(GWh)',round(float(tabular.iloc[0]['AEP_No_Wake']),3), round(float(tabular.iloc[0]['AEP_Baseline']),3), round(float(tabular.iloc[0]['AEP_Opt']),3)), 
        #        ('%', '--', round(float(tabular.iloc[0]['%_Baseline']),3), round(float(tabular.iloc[0]['%_Opt']),3)), 
        #        ('%Wake_Loss', '--',round(float(data.iloc[0]['Wk_Loss_Baseline']),3), round(float(tabular.iloc[0]['Wk_Loss_Opt']),3)),
         #       ('%AEP_Gain', '--', '--', round(float(tabular.iloc[0]['AEP_Gain_Opt']),3)), 
         #       ('Loss_Reduced', '--', '--', round(float(tabular.iloc[0]['Loss_Red_Opt']),3))]
    
        #table_new= pd.DataFrame(farm_data, columns = [' ','No-Wake','Baseline','Optimized'], index= None)
            
        # Render Table using above function 
        #fig, ax = render_mpl_table(table_new)
        
        #table_image = "Table_Image_" + str(kf)+ "_without_unc"
        #plt.savefig(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Fisherman_Study/table_image/{}.png'.format(table_image))
    
    
    else: 
        raise SystemExit("None Valid Optimization Method Chosen")






