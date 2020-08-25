#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 08:57:53 2020

@author: dbensaso
"""
import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.yaw_wind_rose_parallel import YawOptimizationWindRoseParallel
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import WakeSteering_US.namingfarm as nf
import WakeSteering_US.cp_for_any_turb as cturb
import os
import six

if __name__ == '__main__':
    tf= pd.read_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/TI_8_farm_results')
    
    farm_list = tf['Farm Name']
    #### Average wind speed for each farm 
    wind_rose = rose.WindRose()
    farm_avg_ws = pd.DataFrame([])
    
    total_fraction_farm = pd.DataFrame([])
    wind_rose = rose.WindRose()
    tf= pd.read_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/TI_8_farm_results')
    farm_list = tf['Farm Name']
    for farm in farm_list:
        ### First for every farm s
        file_dir = os.path.dirname(os.path.abspath(__file__))
        fi = wfct.floris_interface.FlorisInterface(
            os.path.join(file_dir, '../../example_input.json')
        )
        mf =pd.read_pickle('/home/dbensaso/code/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
                
        kf = (mf.loc[mf['p_name'] == farm])
        wf_coordinate = [kf["ylat"].mean(),kf["xlong"].mean()]
        
        df_opt_pickle = "Df_opt_" + str(kf['p_name'].iloc[0]) + "_without_unc"
        df_opt = pd.read_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Subset_2020/df_Opt_pickle/{}'.format(df_opt_pickle))
        
        file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle_subset_no_ti/{}'.format(file_name))
        
        D = fi.floris.farm.turbines[0].rotor_diameter
        lat_y = kf['ylat'].values
        long_x = kf['xlong'].values
        
        layout_x, layout_y= nf.longlat_to_utm(lat_y, long_x)
        #layout_x=layout_x1.tolist()
        #layout_y=layout_y1.tolist()
        
        file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
        df = wind_rose.load(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle_subset_no_ti/{}'.format(file_name))
            
        
        N_turb = len(layout_x)
        
        fi.reinitialize_flow_field(layout_array=(layout_x, layout_y), wind_direction=[270.0],wind_speed=[8.0])
        fi.reinitialize_flow_field(turbulence_intensity=[0.08]) ### Set turbuelence intensity here 
        fi.calculate_wake()
        #Diameter and Rated power based on the wind farm
        D = kf["t_rd"]
        P_r = kf["t_cap"]
        hub_h = kf["t_hh"]
        C_p_rated = 0.43003137
        C_t_rated = 0.70701647
        #Normalized wind speed for any turbine
        tf= pd.read_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/NREL_5MW_reference')
        
        ## Enumerate so for each turbine 
        for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
                turbine.rotor_diameter = D.iloc[count]
                turbine.hub_height = hub_h.iloc[count]
                T_Area = (np.pi* (D.iloc[count]**2)) /4
                U_turb_rated= (2* P_r.iloc[count]*(10**3)/ (C_p_rated * 1.225* T_Area))**(1/3)
                U_turb_norm =  tf.iloc[:,0] / U_turb_rated
                cp_new = cturb.cp_for_any_turb(U_turb_norm,U_turb_rated,T_Area,P_r.iloc[count],tf)
                ct_new = cturb.ct_for_any_turb(U_turb_norm,tf)
                turbine.power_thrust_table["power"] = cp_new
                turbine.power_thrust_table["thrust"] = ct_new
                turbine.change_turbine_parameters({})
        for count, coord in enumerate(fi.floris.farm.flow_field.turbine_map.coords):
            coord.x3 = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height
        fi.floris.farm.flow_field.specified_wind_height = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height
        
        fraction = pd.DataFrame([])
        for turb in range(len(fi.floris.farm.flow_field.turbine_map.turbines)):
            turbine_id = fi.floris.farm.flow_field.turbine_map.turbines[turb]
            freq_yaw = []
            for i in range(len(df_opt['wd'])):
                turbine_yaw = df_opt['yaw_angles'][i][turb]
                freq = df['freq_val'][i]
                #if abs(turbine_yaw) >=1:
                #freq_yaw.append(freq)
                total_fraction_farm = total_fraction_farm.append(pd.DataFrame({'Farm Name': farm,
                                                    'yaw': turbine_yaw,
                                                    'freq':freq}, 
                                                    index=[0]), ignore_index=True)
                
    yaw_pickle = "Yaw_dist_8_without_unc"
    total_fraction_farm.to_pickle(r'/home/dbensaso/code/floris/examples/optimization/scipy/Saved_Fig/Subset_2020/yaw_dist_test/{}'.format(yaw_pickle))
            
                #total_turbine_yaw_freq = sum(freq_yaw) * 100
                #fraction = fraction.append(pd.DataFrame({'Farm Name': farm,
                 #                                        'yaw_freq': total_turbine_yaw_freq}, 
                  #                                      index=[0]), ignore_index=True)
            
            #total_fraction_farm = total_fraction_farm.append(fraction,ignore_index=True)
