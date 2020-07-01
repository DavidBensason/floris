#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:05:30 2020

@author: dbensaso
"""

### Making power curves on the the D changing data for a single turbine in a farm to see baseline power 

import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
#import floris.tools.cut_plane as cp
from floris.tools.optimization.scipy.yaw_wind_rose import YawOptimizationWindRose
import floris.tools.wind_rose as rose
import WakeSteering_US.WakeSteering_US.cp_for_any_turb as cturb
#import floris.tools.power_rose as pr
import numpy as np
import pandas as pd
import os
import math
#import six

# Instantiate the FLORIS object
file_dir = os.path.dirname(os.path.abspath(__file__))
fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, '../../../example_input.json')
)
"""
cp_15MW =[0.100335552,0.272197208,0.359305118,0.409156685,0.42729002,0.441819614,0.452708049,0.475270298,
0.479478644,0.482871263,0.484447627,0.485186337,0.485905846,0.487114154,0.488094793,0.488866488,
0.488988058,0.489050152,0.489088858,0.489114531,0.489127764,0.48915576,0.489165047,0.489154466,
0.489163867,0.489224161,0.489263048,0.48928802,0.489293116,0.489304304,0.489319143,0.48938292,
0.489378059,0.486507177,0.472991558,0.470349316,0.467724941,0.465119143,0.462537044,0.462021882,
0.461765284,0.461636557,0.46150846,0.461380941,0.461317766,0.459972382,0.447424184,0.435324299,
0.406934831,0.380965216,0.357157168,0.335295525,0.263718683,0.211136579,0.171651956,0.108080744,
0.072394937,0.050842443,0.037062292]

ct_15MW = [0.819748943,0.801112031,0.808268424,0.821910918,0.822876237,0.823265981,0.830989358,0.834932456,
0.833618598,0.83180478,0.829011103,0.826909201,0.824740997,0.820429675,0.816176257,0.811200233,
0.809740903,0.808780765,0.808102306,0.807566626,0.807251977,0.80662442,0.806495512,0.806806173,
0.806651158,0.805469658,0.804571567,0.803949121,0.803904895,0.803708734,0.80345211,0.801706154,
0.801777393,0.768657554,0.70731525,0.698507743,0.690211963,0.682335591,0.674835939,0.673371183,
0.672646111,0.672283185,0.671921569,0.671564033,0.671386994,0.667639697,0.635292304,0.607277698,
0.548965866,0.501379105,0.460982977,0.425965654,0.32116631,0.2511023,0.201415182,0.125653944,
0.08506697,0.061026446,0.045814967]

ws_15 = [2.999999831,3.499999916,4,4.500000084,4.750000126,5.000000169,5.249999874,5.999999663,
6.199999966,6.40000027,6.499999747,6.55000016,6.599999899,6.700000051,6.800000202,6.900000354,
6.919999845,6.929999928,6.94000001,6.950000093,6.960000175,6.969999584,6.980000341,6.989999749,
6.999999831,7.499999916,8,8.500000084,9.000000169,9.500000253,10.00000034,10.2499997,10.49999975,
10.60000057,10.70000005,10.72000022,10.73999971,10.76000055,10.78000004,10.78400034,10.78600049,10.78699989,
10.78799997,10.78899937,10.78950042,10.8000002,10.89999968,10.99999983,
11.24999987,11.50000059,11.75000063,11.99999933,12.99999949,13.99999966,14.99999983,17.50000025,20.00000067,
22.49999975,24.99999882]
"""
D_ls = [100]
ws_list = np.arange(0,26,1)
#overall_1 = pd.DataFrame([])
overall_2 = pd.DataFrame([])
P_r = 2000
hub_h = 95
y_n = 3
#C_p_rated = 0.31
#C_t_rated = 0.45
C_p_rated = 0.43003137
C_t_rated = 0.70701647

cP_test = pd.DataFrame([])
#tf= pd.read_pickle(r'/home/dbensaso/floris/examples/optimization/scipy/Lookup_table_8MW')
tf= pd.read_pickle('/home/dbensaso/floris/examples/optimization/scipy/NREL_5MW_reference')
for ii in D_ls:
    power_tot = pd.DataFrame([])
    for i in ws_list:
        
        for count, turbine in enumerate(fi.floris.farm.flow_field.turbine_map.turbines):
            turbine.rotor_diameter = int(ii)
            turbine.hub_height = hub_h
            T_Area = (np.pi* (int(ii)**2)) /4
            U_turb_rated= (2* P_r*(10**3)/ (C_p_rated * 1.225* T_Area))**(1/3)
            U_turb_norm =  tf.iloc[:,0] / U_turb_rated
            cp_new_2 = cturb.cp_for_any_turb(U_turb_norm,U_turb_rated,T_Area,P_r,tf)
            ct_new_2 = cturb.ct_for_any_turb(U_turb_norm,tf)
            turbine.power_thrust_table["power"] = cp_new_2
            turbine.power_thrust_table["thrust"] = ct_new_2
            #turbine.power_thrust_table["wind_speed"] = ws_15
            #turbine.generator_efficiency = 0.93
            turbine.change_turbine_parameters({})
        for count, coord in enumerate(fi.floris.farm.flow_field.turbine_map.coords):
            coord.x3 = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height
            #turbine.change_turbine_parameters({})
        fi.floris.farm.flow_field.specified_wind_height = fi.floris.farm.flow_field.turbine_map.turbines[0].hub_height
        
        D = ii
        N_row =1 
        T_row = 1
        layout_x = []
        layout_y = []
        Num_Turb = N_row*T_row
        #y_n = 1    ## For non-constant area option, spacing between N (vertical)
        #x_t = 1   ## For non-constant area option, spacing between T (horizontal)
        constant_area_layout = False
        #relative_original_spacing = False    #y_n and x_t are multiples of (1852/D)=11.29
        if constant_area_layout: ##THIS WORKS ATM 
            spc_N = (1852/D) *(13/N_row)  #(1Nautical mile/164)
            spc_T= (1852/D) * (8/T_row) #(1Nautical mile/164)
           
        else: 
            spc_N = y_n
                    
        for j in range(N_row):
            for k in range(T_row):
                layout_x.append(j*y_n*D*math.cos(-45) - k*y_n*D*math.cos(-45))
                layout_y.append(j*y_n*D*math.cos(-45) + k*y_n*D*math.cos(-45))
        
        fi.reinitialize_flow_field(layout_array=(layout_x, layout_y),wind_direction=[270.0],wind_speed=[i])
        fi.calculate_wake()
        
        power = fi.get_turbine_power()
        
        #if power >= P_r:
        #cP_test = cP_test.append(pd.DataFrame({'D': int(ii), 'cp': [cp_new], 'ct':[ct_new]}, 
        #                                 index=[0]), ignore_index=True)    
        power_tot['power_col'] = power
        power_tot['power_col'] = power_tot['power_col'].apply(lambda x: x/10**6)
        power_tot['ws'] = float(i)
        power_tot['D'] = float(ii)
        #power_tot['Rated'] = U_turb_rated
        #power_tot['ct'] = [ct_new]
        overall_2 = overall_2.append(power_tot, ignore_index = True)
        
        #self.air_density
        #self.rotor_radius
        #self.generator_efficiency
        #cptmp

#d_120 = overall.loc[overall['D'] == 120]
#d_140 = overall.loc[overall['D'] == 140]
#d_160 = overall.loc[overall['D'] == 160]
#r_5000 = overall_1.loc[overall_1['D'] == 126]
r_2000 = overall_2.loc[overall_2['D'] == 100]

#plt.plot(d_120['ws'], d_120['power_col'], '-', label= '120')
#plt.plot(d_140['ws'], d_140['power_col'], '-o', label= '140')
#plt.plot(d_160['ws'], d_160['power_col'], '-o', label= '160')
"""

### THIRD 
fig=plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.2)
  
plt.plot(r_5000['ws'], r_5000['power_col'], '-b',linewidth=3, label= 'NREL 5MW Reference')
plt.plot(r_2000['ws'], r_2000['power_col'], '--k', label= 'Border Winds Project')

plt.title('c) Turbine Power Curves', fontsize=15)
plt.xlabel('Wind Speed (m/s$^2$)', fontsize=10)
plt.ylabel('Power (MW)',fontsize=10)
#plt.legend(loc= 'upper left')
plt.grid(True)
table_power = "Power_for_report"
plt.savefig(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Journal_images/{}.png'.format(table_power))
 
plt.show()

##FIRST
fig=plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.2)

plt.plot(tf['Set_Wind_Speed'], cp_new_1, '-b',linewidth=3, label = 'NREL 5MW Reference')
plt.plot(tf['Set_Wind_Speed'], cp_new_2, '--k', label = 'Border Winds Project')
plt.title('a) Turbine cp Curves', fontsize=15)

plt.xlabel('Wind Speed (m/s$^2$)', fontsize=10)
plt.ylabel('cp',fontsize=10)
plt.legend()
plt.grid(True)
table_cp = "Cp_for_report"
plt.savefig(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Journal_images/{}.png'.format(table_cp))

plt.show()
   

## SECOND 
fig=plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.2)
  
plt.plot(tf['Set_Wind_Speed'], ct_new_1, '-b',linewidth=3, label = 'NREL 5MW Reference')
plt.plot(tf['Set_Wind_Speed'], ct_new_2, '--k', label = 'Border Winds Project')
plt.title('b) Turbine ct Curves', fontsize=15)
plt.xlabel('Wind Speed (m/s$^2$)', fontsize=10)
plt.ylabel('ct',fontsize=10)
plt.grid(True)
table_ct = "Ct_for_report"
plt.savefig(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/Journal_images/{}.png'.format(table_ct))
 
plt.show()









cp_160 = cP_test.loc[cP_test['D'] == 160]
cp_180 = cP_test.loc[cP_test['D'] == 180]
cp_200 = cP_test.loc[cP_test['D'] == 200]

#plt.plot(d_120['ws'], d_120['power_col'], '-', label= '120')
#plt.plot(d_140['ws'], d_140['power_col'], '-o', label= '140')
plt.plot(tf['Set_Wind_Speed'], cp_new_160, '-o', label= '160')
plt.plot(tf['Set_Wind_Speed'], cp_new_180, '-o', label= '180')
plt.plot(tf['Set_Wind_Speed'], cp_new_200, '-o', label= '200')
plt.legend(loc= 'upper right')
plt.show()

"""
"""
plt.plot(d_120['D'], d_120['Rated'], 'o', label= '120')
plt.plot(d_140['D'], d_140['Rated'], 'o', label= '140')
plt.plot(d_160['D'], d_160['Rated'], 'o', label= '160')
plt.plot(d_180['D'], d_180['Rated'], 'o', label= '180')
plt.plot(d_200['D'], d_200['Rated'], 'o', label= '200')
plt.legend(loc= 'upper right')
plt.show()


plt.plot(tf['Set_Wind_Speed'], fi.floris.farm.flow_field.turbine_map.turbines[0].power_thrust_table["power"], '-o', label= 'Adjusted for 200')
plt.plot(tf['Set_Wind_Speed'], tf['C_p'], '-o', label = 'Standard 8MW')
plt.legend(loc= 'best')
plt.show()

plt.plot(tf['Set_Wind_Speed'], ct_new, '-o', label= 'Adjusted for 200')
plt.plot(tf['Set_Wind_Speed'], tf['C_t'], '-o', label = 'Standard 8MW')
plt.legend(loc= 'best')
plt.show()
"""