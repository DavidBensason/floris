#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 12:47:43 2020

@author: dbensaso
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:11:35 2020

@author: dbensaso
"""
from itertools import combinations 
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
import math
# Instantiate the FLORIS object
file_dir = os.path.dirname(os.path.abspath(__file__))
fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, '../../example_input.json')
)


mf =pd.read_pickle('/home/dbensaso/WakeSteering_US/Working_dir_WS_US/Wind_US_Database')

#kf1 = (mf.loc[mf['p_name'] == "Grant Wind"])
#kf2 = (mf.loc[mf['p_name'] == "Grant Plains"])
#kf = pd.concat([kf1,kf2], ignore_index=True)
#mf =pd.read_pickle('/home/dbensaso/floris/examples/optimization/scipy/Peetz_Logan_Data')
kf = (mf.loc[mf['p_name'] == "Burley Butte"])
#kf1 = (mf.loc[mf['p_name'] == "Chisholm View"])
#kf2 = (mf.loc[mf['p_name'] == "Chisholm View II"])
#kf = pd.concat([kf1,kf2], ignore_index=True)

opt_options = {
    "maxiter": 5,
    "disp": True,
    "iprint": 2,
    "ftol": 1e-7, #     "eps": 0.01,
}

#kf = (mf.loc[mf['p_name'] == "Wildcat Ranch"])
wf_coordinate = [kf["ylat"].mean(),kf["xlong"].mean()]
#wf_coordinate = [40.995293, -84.565796]
# Set wind farm to N_row x N_row grid with constant spacing 
# (2 x 2 grid, 5 D spacing)
D = fi.floris.farm.turbines[0].rotor_diameter
lat_y = kf['ylat'].values
long_x = kf['xlong'].values

layout_x, layout_y= nf.longlat_to_utm(lat_y, long_x)

#layout_x=layout_x1.tolist()
#layout_y=layout_y1.tolist()

N_turb = len(layout_x)
fi.reinitialize_flow_field(layout_array=(layout_x, layout_y), wind_direction=[270],wind_speed=[10])
fi.calculate_wake()

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
asas

# Identify which turbines wake eachother 
def sub_farms_1(wd):
    waked_list=fi.floris.farm.turbine_map.number_of_turbines_iec(wd)
    
    # Convert waked_list to non-tuple list 
    new_list = [i[0] for i in waked_list]
    
    # Generate combinations for each turbine waked pair 
    combination_list = []
    for i in range(len(new_list)):
        comb = list(combinations(new_list[i], 2))
        for j in range(len(comb)):
            if comb[j][0] ==i or comb[j][1] ==i:
                combination_list.append(comb[j])
                
    # Generate x and y corrdinate that match the index values of the items 
    x_cord = []
    y_cord = []
    for i, (coord, turbine) in enumerate (fi.floris.farm.flow_field.turbine_map.items):
        x_cord.append(coord.x1)
        y_cord.append(coord.x2)
    
    # Generate edges for each combination and return list of farm groups 
    #if __name__=="__main__": #Not sure if this is required 
    g = Graph(len(x_cord))
    for k in combination_list:
        g.addEdge(k[0],k[1])
    cc = g.connectedComponents()
    
    return cc,x_cord, y_cord, combination_list


#Python program to print connected  
# components in an undirected graph 
class Graph: 
      
    # init function to declare class variables 
    def __init__(self,V): 
        self.V = V
        self.adj = [[] for i in range(V)] 
  
    def DFSUtil(self, temp, v, visited): 
  
        # Mark the current vertex as visited 
        visited[v] = True
  
        # Store the vertex to list 
        temp.append(v) 
  
        # Repeat for all vertices adjacent 
        # to this vertex v 
        for i in self.adj[v]: 
            if visited[i] == False: 
                  
                # Update the list 
                temp = self.DFSUtil(temp, i, visited) 
        return temp 
  
    # method to add an undirected edge 
    def addEdge(self, v, w): 
        self.adj[v].append(w) 
        self.adj[w].append(v) 
  
    # Method to retrieve connected components 
    # in an undirected graph 
    def connectedComponents(self): 
        visited = [] 
        cc = [] 
        for i in range(self.V): 
            visited.append(False) 
        for v in range(self.V): 
            if visited[v] == False: 
                temp = [] 
                cc.append(self.DFSUtil(temp, v, visited)) 
        return cc 

wind_rose = rose.WindRose()
file_name = str(kf['p_name'].iloc[0]) + "_Wind_Farm.p"
df = wind_rose.load(r'/home/dbensaso/floris/examples/optimization/scipy/Saved_Fig/wind_rose_pickle_subset_no_ti/{}'.format(file_name))
               
overall = pd.DataFrame([])
wd_list = np.arange(0,360,5)
for wind_dir in wd_list:
    wd_case = df[df['wd'] == wind_dir]
    freq_sum = wd_case['freq_val'].sum()
    cc,x_cord,y_cord, combination_list =sub_farms_1(wind_dir)
    
    # Drop duplicates
    res = [] 
    for icomb in combination_list: 
        if icomb not in res: 
            res.append(icomb)
    
    for j in range(len(cc)): #through each sub farm 
        sub_farm_xcoord = []
        sub_farm_ycoord = []		
        for k in cc[j]:
            sub_farm_xcoord.append(x_cord[k])
            sub_farm_ycoord.append(y_cord[k])
    
    rel_distances = []
    for j in res: #through each sub farm 
        sub_farm_xcoord = []
        sub_farm_ycoord = []		
        for k in j:
            sub_farm_xcoord.append(x_cord[k])
            sub_farm_ycoord.append(y_cord[k])
            
        distance = math.sqrt( ((sub_farm_xcoord[0]-sub_farm_xcoord[1])**2)+((sub_farm_ycoord[0]-sub_farm_ycoord[1])**2) ) / D
        
        rel_distances.append(distance)
    
    avg_rel_dist = sum(rel_distances) / len(rel_distances)
    
    product = avg_rel_dist * freq_sum
    overall = overall.append(pd.DataFrame({'wd': wind_dir,'avg_rel_spc': avg_rel_dist, 'freq_wd':freq_sum,
                                           'product': product}, 
                                                 index=[0]), ignore_index=True)

overall['product'].sum()


## Method, just find relative shortest distance to nearest turbine for each turb in farm and find average
node_list = []
from scipy.spatial import distance

def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]

for i in range(len(layout_x)):
    x = layout_x[i]
    y = layout_y[i]
    point = (x,y)
    node_list.append(point)

rel_distances = []
for coord in range(len(node_list)):
    to_exclude = {coord}
    nodes2 = [element for ii, element in enumerate(node_list) if ii not in to_exclude]
    close_p = closest_node(node_list[coord], nodes2)
    distances = math.sqrt( ((node_list[coord][0]-close_p[0])**2)+((node_list[coord][1]-close_p[1])**2) )
    rel_distance = distances / D
    rel_distances.append(rel_distance)

avg_rel_spc = sum(rel_distances) / len(rel_distances)