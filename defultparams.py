# -*- coding: utf-8 -*-
from math import sqrt

################## Model params############################
X0 = 1.e4  # Length of bay
Y0 = X0 / 20.0  # Half-width of bay
sigmaG = X0 / sqrt(2.)  # Gauss sigma value (at t=0)
eta0 = 2.0  # Gauss amplitude at t=0
Z0 = -2.e2  # Open sea depth

################## Step cotrol params ########
p_stepCheckSW = True  # Check S-W conditions
p_stepMesh = True  # Create mesh
p_stepSWW = True  # Execute simulation
p_stepReport = True  # Extract data from SWW
p_stepGauge = True  # Generate gauge file
p_stepGaugeChart = True  # Extract data for Gauge points ant generate charts

################## Time params ###################
timestep1 = 5.  # 1st stage time step
timestep2 = 1.  # 2nd stage time step
duration1 = 500  # 1st stage duration
duration2 = 1400  # 2nd stage duration
################## Resolution #####################################
max_tri_area_base = 5.e3  # Default max. tirangle area (m^2)
max_tri_area_focus = 5.e2  # Max triangle area (bay + reflective wall)
max_tri_area_inu = 5.e1  # Max triangle area (Run-up region)
###################################################################
flow_alg = 'DE2'  # Algorithm
p_friction = 0.  # Friction
infgcx = 4.0  # Gauss initial central point( acc. to bay-toe )
p_minimum_storable_height = 0.001  # Don't record water less than
ww = 60  # wall width ( if wall exist)
###################################################################
p_cache = False
p_verbose = True
###################################################################
elv_max=30. # Expected maximum elevation. Performance issue