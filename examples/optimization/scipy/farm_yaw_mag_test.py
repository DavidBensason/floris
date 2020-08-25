#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 11:18:32 2020

@author: dbensaso
"""
import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.yaw_wind_rose import YawOptimizationWindRose
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import WakeSteering_US.WakeSteering_US.namingfarm as nf
import WakeSteering_US.WakeSteering_US.cp_for_any_turb as cturb
import pdb
import os
import time
import six
### Plotting magnitude of yaw angle in the farm

file_dir = os.path.dirname(os.path.abspath(__file__))
fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, '../../example_input.json')
)

#df_base =pd.read_pickle('/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Subset_2020/df_Base_pickle/Df_base_Bear Creek_without_unc')

mf =pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')
data = pd.DataFrame([])
data1 = pd.DataFrame([]) #Only filled for Unc case
#for i in dw:
wind_rose = rose.WindRose()

kf = (mf.loc[mf['p_name'] == "Border Winds Project"])
wf_coordinate = [kf["ylat"].mean(),kf["xlong"].mean()]

df_opt_pickle = "Df_opt_" + str(kf['p_name'].iloc[0]) + "_without_unc"
df_opt = pd.read_pickle(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Subset_2020/df_Opt_pickle/{}'.format(df_opt_pickle))
df_base_pickle = "Df_base_" + str(kf['p_name'].iloc[0]) + "_without_unc"
df_base = pd.read_pickle(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Subset_2020/df_Base_pickle/{}'.format(df_base_pickle))
file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
df = wind_rose.load(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle_subset_no_ti/{}'.format(file_name))
    
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
fi.reinitialize_flow_field(turbulence_intensity=[0.08]) ### Set turbuelence intensity here 
fi.calculate_wake()
opt_options = {'maxiter': 20, 'disp': False,'iprint': 1, 'ftol': 1e-6, 'eps': 0.01}

#Diameter and Rated power based on the wind farm
D = kf["t_rd"]
P_r = kf["t_cap"]
hub_h = kf["t_hh"]

C_p_rated = 0.43003137
C_t_rated = 0.70701647

#Normalized wind speed for any turbine
tf= pd.read_pickle(r'/home/dbensaso/floris/examples/optimization/scipy/NREL_5MW_reference')

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

# set min and max yaw offsets for optimization 
min_yaw = -25.0
max_yaw = 25.0

# Define minimum and maximum wind speed for optimizing power. 
# Below minimum wind speed, assumes power is zero.
minimum_ws = 3.0
maximum_ws = 15.0
asd

## For a single wd and ws combination case: Lets do 240 degress at 10m/ second 
overall_second_case = pd.DataFrame([])
c1 = df_opt[(df_opt['wd'] == 240) & \
                  (df_opt['ws'] == 10)]
c2 = df_base[(df_base['wd'] == 240) & \
                  (df_base['ws'] == 10)]
for turb in range(len(fi.floris.farm.flow_field.turbine_map.turbines)):
    turbine_id = fi.floris.farm.flow_field.turbine_map.turbines[turb]
    x_coord = fi.floris.farm.flow_field.turbine_map.coords[turb].x1
    y_coord = fi.floris.farm.flow_field.turbine_map.coords[turb].x2
    turbine_val = []
    turbine_opt_p = c1['turbine_power_opt'].iloc[0][turb]
    turbine_baseline_p = c2['turbine_power_baseline'].iloc[0][turb]
    turbine_yaw = abs(c1['yaw_angles'].iloc[0][turb])
    gain = ((turbine_opt_p - turbine_baseline_p)/ turbine_baseline_p) * 100 
    #freq = df['freq_val'][i]
    #product = abs(turbine_yaw * freq)
    overall_second_case = overall_second_case.append(pd.DataFrame({'turbine_ID': turbine_id,
                                                                   'xcoord': x_coord,
                                                                   'ycoord': y_coord,
                                                                   'yaw': turbine_yaw,
                                                                   'opt_power': turbine_opt_p,
                                                                   'percent_gain': gain}, 
                                                     index=[0]), ignore_index=True)

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.scatter(overall_second_case['xcoord'], overall_second_case['ycoord'], c = overall_second_case['yaw'], s = 60, cmap='bwr')
ax.set_xlabel('Easting (m)')
ax.set_ylabel('Norming (m)')
#ax.set_facecolor('xkcd:silver')
plt.colorbar(label='Yaw Misalaignment Magnitude')

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.scatter(overall_second_case['xcoord'], overall_second_case['ycoord'], c = overall_second_case['percent_gain'], s = 60, cmap='RdYlGn')
#ax.set_facecolor('xkcd:red')
ax.set_xlabel('Easting (m)')
ax.set_ylabel('Norming (m)')
plt.colorbar(label='% Gain')



####################################################################################
## Contour plot idea from papa
Z = df.pivot_table(index='wd', columns='ws', values='freq_val').T.values

X_unique = np.sort(df.wd.unique())
Y_unique = np.sort(df.ws.unique())
X, Y = np.meshgrid(X_unique, Y_unique)

fig,ax=plt.subplots(1,1)
cp1 = ax.contourf(X, Y, Z, 50)

v = fig.colorbar(cp1) # Add a colorbar to a plot
v.set_label('frequency')

ax.set_title('Case Frequency')
ax.set_xlabel('wd')
ax.set_ylabel('ws')
plt.show()

###################################################################################3
#### looking at sum of correction magnitudes at each direction 
part_b = pd.DataFrame([])
for case in range(len(df_opt)):
    freq_val = df['freq_val'][case]
    ws = df['ws'][case]
    wd = df['wd'][case]
    yaw_angles = df_opt['yaw_angles'].iloc[case]
    res =  [abs(ele) for ele in yaw_angles] 
    yaw_angles_sum = sum(res)
    yaw_angles_mean = sum(res)/ len(res)
    cases = []
    for i in res:
        if i >=5:
            cases.append(i)
        num_turb = len(cases)
        #if len(cases) > 0:
            #average = sum(cases) / len(cases)
    part_b = part_b.append(pd.DataFrame({'ws': ws,
                                                      'wd': wd,
                                                      'freq':freq_val,
                                                      'yaw_sum':yaw_angles_sum, 
                                                      'yaw_mean':yaw_angles_mean,
                                                      'case_above':num_turb}, 
                                                     index=[0]), ignore_index=True)

## Contour plot idea from papa
Z = part_b.pivot_table(index='wd', columns='ws', values='yaw_sum').T.values

X_unique = np.sort(df.wd.unique())
Y_unique = np.sort(df.ws.unique())
X, Y = np.meshgrid(X_unique, Y_unique)

fig,ax=plt.subplots(1,1)
cp1 = ax.contourf(X, Y, Z, 50)

v = fig.colorbar(cp1) # Add a colorbar to a plot
v.set_label('yaw sum')

ax.set_title('Sum of farm yaw angles')
ax.set_xlabel('wd')
ax.set_ylabel('ws')
plt.show()

Z = part_b.pivot_table(index='wd', columns='ws', values='yaw_mean').T.values

X_unique = np.sort(df.wd.unique())
Y_unique = np.sort(df.ws.unique())
X, Y = np.meshgrid(X_unique, Y_unique)

fig,ax=plt.subplots(1,1)
cp1 = ax.contourf(X, Y, Z, 50)

v = fig.colorbar(cp1) # Add a colorbar to a plot
v.set_label('yaw mean')

ax.set_title('Mean of farm yaw angles')
ax.set_xlabel('wd')
ax.set_ylabel('ws')
plt.show()

Z = part_b.pivot_table(index='wd', columns='ws', values='case_above').T.values

X_unique = np.sort(df.wd.unique())
Y_unique = np.sort(df.ws.unique())
X, Y = np.meshgrid(X_unique, Y_unique)

fig,ax=plt.subplots(1,1)
cp1 = ax.contourf(X, Y, Z, 50)

v = fig.colorbar(cp1) # Add a colorbar to a plot
v.set_label('Number of turbines with more than 5 degree')

ax.set_title('Number of turbines above 5')
ax.set_xlabel('wd')
ax.set_ylabel('ws')
plt.show()
##########################################################################################


wind_rose = rose.WindRose()
file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
df = wind_rose.load(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle_subset_no_ti/{}'.format(file_name))
wind_rose.plot_wind_rose()

wd_list = np.arange(0,360,5)

overall_second_case = pd.DataFrame([])
for turb in range(len(fi.floris.farm.flow_field.turbine_map.turbines)):
    turbine_id = fi.floris.farm.flow_field.turbine_map.turbines[turb]
    turbine_val = []
    for i in range(len(df_opt['wd'])):
        turbine_yaw = df_opt['yaw_angles'][i][turb]
        freq = df['freq_val'][i]
        product = abs(turbine_yaw * freq)
        turbine_val.append(product)
    
    total_for_turbine = sum(turbine_val)
    overall_second_case = overall_second_case.append(pd.DataFrame({'turbine_ID': turbine_id,
                                                                   'weighted_yaw': total_for_turbine}, 
                                                     index=[0]), ignore_index=True)     

### Get the x and y coordinates for each turbine id in the saem order 
coords_of_turbines =pd.DataFrame([])

for i in range(len(fi.floris.farm.flow_field.turbine_map.turbines)):
    x_coord = fi.floris.farm.flow_field.turbine_map.coords[i].x1
    y_coord = fi.floris.farm.flow_field.turbine_map.coords[i].x2
    
    coords_of_turbines = coords_of_turbines.append(pd.DataFrame({'xcoord': x_coord,'ycoord': y_coord}, 
                                                 index=[0]), ignore_index=True)


final_yaw_data = pd.DataFrame([])

final_yaw_data['Turbine_ID'] = overall_second_case['turbine_ID']
final_yaw_data['weighted_yaw'] = overall_second_case['weighted_yaw']
final_yaw_data['xcoord'] = coords_of_turbines['xcoord']
final_yaw_data['ycoord'] = coords_of_turbines['ycoord']



fig, ax = plt.subplots(nrows=1, ncols=1)

plt.scatter(final_yaw_data['xcoord'], final_yaw_data['ycoord'], c = final_yaw_data['weighted_yaw'], s = 60, cmap='Blues')
ax.set_facecolor('xkcd:red')
plt.colorbar(label='Average Yaw')


### Look at power defecit in turbines 
power_second_case = pd.DataFrame([])
for turb in range(len(fi.floris.farm.flow_field.turbine_map.turbines)):
    turbine_id = fi.floris.farm.flow_field.turbine_map.turbines[turb]
    tot_p_opt = []
    tot_b_base = []
    for i in range(len(df_opt['wd'])): 
        freq = df['freq_val'][i]
        turb_power_opt = df_opt['turbine_power_opt'][i][turb] * freq
        turb_power_baseline = df_base['turbine_power_baseline'][i][turb] * freq
        tot_p_opt.append(turb_power_opt)
        tot_b_base.append(turb_power_baseline)
    total_opt = sum(tot_p_opt)
    total_base = sum(tot_b_base)
    percentage_inc = ((total_opt - total_base) / total_base) *100
    power_second_case = power_second_case.append(pd.DataFrame({'turbine_ID': turbine_id,
                                                                       'per_power_dif': percentage_inc}, 
                                                         index=[0]), ignore_index=True)

final_power_data = pd.DataFrame([])

final_power_data['Turbine_ID'] = power_second_case['turbine_ID']
final_power_data['per_power_dif'] = power_second_case['per_power_dif']
final_power_data['xcoord'] = coords_of_turbines['xcoord']
final_power_data['ycoord'] = coords_of_turbines['ycoord']

fig, ax = plt.subplots(nrows=1, ncols=1)

plt.scatter(final_power_data['xcoord'], final_power_data['ycoord'], c = final_power_data['per_power_dif'], s = 60, cmap='Blues')
ax.set_facecolor('xkcd:red')
plt.colorbar(label='Perent Diff power')