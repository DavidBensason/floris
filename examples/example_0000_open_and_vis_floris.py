# Copyright 2019 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See read the https://floris.readthedocs.io for documentation

import matplotlib.pyplot as plt
import floris.tools as wfct
from floris.simulation.farm import Farm
import copy
# Initialize the FLORIS interface fi
fi = wfct.floris_interface.FlorisInterface("example_input.json")
farm_instance_dict = fi.floris.farm.instance_dictionary

upper = copy.deepcopy(farm_instance_dict)[0]
upper['properties']['layout_x'] = [0, 750]
upper['properties']['layout_y'] = [600,600]
lower = copy.deepcopy(farm_instance_dict)[0]
lower['properties']['layout_x'] = [0, 750]
lower['properties']['layout_y'] = [0, 0]
farm_up = Farm(upper, fi.floris.farm.flow_field.turbine_map.turbines[0], fi.floris.farm.flow_field.wake)
farm_low = Farm(lower, fi.floris.farm.flow_field.turbine_map.turbines[0], fi.floris.farm.flow_field.wake)
farm_up.flow_field.calculate_wake()
farm_low.flow_field.calculate_wake()

# Calculate wake
fi.calculate_wake()
print(fi.floris.farm.wind_direction)

# Get horizontal plane at default height (hub-height)
hor_plane = fi.get_hor_plane()

# Plot and show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
plt.show()
