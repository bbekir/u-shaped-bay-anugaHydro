# -*- coding: utf-8 -*-
"""Created on Fri Apr 05 21:55:47 2019
@author: Ebubekir Çelik
Analysis on u shaped bay
Bay equation: z=-a*( (x-xk)-(X0/Y0**2)*(y-yk)**2 )
Bay eq. at x=x0 --> z=-a*X0(1-y**2/Y0**2)  --> y= sqrt(z/Z0 +1)*Y0

To enable/disable steps change the values of p_step.. variables under "Runtime Behaviour" section

Steps:
- Preparation  : Calculates/assigns variables, coordinates etc. Prepare directory structure for output
- SW Control   : Check Shallow Water conditions.Can be disabled if not needed
- Create Mesh  : Creates mesh with .msh extension, Creates domain object. Records configuration to log file  log1
- Simülation   : Runs the simulation. Generates sww file and log
- Report       : Exracts  max-runup time and coordinate from output. Records calculation in log3 file. Extract snapshot of max-runup time
- Gauge        : Prepares appropriate point list to track time-dependent measurables for comparison  purpose and ceates gauge file
- Gauge Charts : Extracts time-varying values for coordinates at gauge file that  prepared at Gauge step. Creates seperate csv files for every coordinates
                  Creates time-series plots from csv files
"""

from math import sqrt, ceil
from anuga_util.spcutil import dirStruct
import time
import anuga
import anuga_util.utility as util
import os
# import sys
import defultparams as default
import numpy as np
import matplotlib.pyplot as plt
import anuga_util.animate as animate

# run_note = raw_input("Enter definitive comment.  ")
# prj_suffix =raw_input( "Enter suffix for run")
run_note = " default çözünülük*10."
prj_suffix = '02'
prj_prefix = r'modeldefult'
# ------------------------------------------------------------------------------
# Runtime behaviour
# ------------------------------------------------------------------------------
p_cache = default.p_cache
p_verbose = default.p_verbose
p_stepCheckSW = True  # default.p_stepCheckSW
p_stepMesh = True  # default.p_stepMesh
p_stepSWW = True  # default.p_stepSWW
p_stepReport = True  # default.p_stepReport
p_stepGauge = True  # default.p_stepGauge
p_stepGaugeChart = True  # default.p_stepGaugeChart

# ------------------------------------------------------------------------------
# Bay geometry
# ------------------------------------------------------------------------------

# Min depth of bay -at X0. !!Negative . type: float
Z0 = default.Z0

# Length of bay. type: float
X0 = default.X0

# Half-width of bay .type: float
Y0 = default.Y0

# Initial gauss shaped stage sigma value
sigmaG = default.sigmaG

# Initial gauss stage max height ( at the center of gauss) . type: float
eta0 = default.eta0

# Friction
p_friction = default.p_friction

# Minimum storable stage
p_minimum_storable_height = default.p_minimum_storable_height

# ------------------------------------------------------------------------------
# Execution parameters
# ------------------------------------------------------------------------------

# Flow algorithm
flow_alg = default.flow_alg

# Defult maximum triangle size (m2). type: float
max_tri_area_base = default.max_tri_area_base

# Focused area maximum triangle size (m2). type: float
max_tri_area_focus = default.max_tri_area_focus

# Inundation area maximum triangle size (m2). type: float
max_tri_area_inu = default.max_tri_area_inu

# Simulation timestep (sec)
timestep1 = default.timestep1
timestep2 = default.timestep2

# Simulation duration (sec)
duration1 = default.duration1
duration2 = default.duration2

# Approaching..
# To be used for calculating gauss stage center position gauss_x : X0 +xk +sigmaG*infgcx
infgcx = default.infgcx
elv_max = default.elv_max

###END OF PARAMETER SECTION

# ------------------------------------------------------------------------------
# Define file-directory structure
# ------------------------------------------------------------------------------

main_output_dir = os.getcwd() + os.sep + 'output'
domain_name = prj_prefix + '_' + prj_suffix
prj_output_dir = main_output_dir + os.sep + domain_name
pts_file = domain_name + '.pts'
mesh_file = domain_name + '.msh'
sww_file = domain_name + '.sww'
plot_file = domain_name + '.png'
gauge_file = domain_name + '_gauges.csv'
log1_file = domain_name + '_log01.log'
log2_file = domain_name + '_log02.log'
log3_file = domain_name + '_log03.log'

# ------------------------------------------------------------------------------
# Calculations
# ------------------------------------------------------------------------------

xk = ceil(-(X0 / Z0) * eta0 * elv_max)  # extra space (X direction)
yk = Y0 * 1.2  # ceil (sqrt(1 + elv_max/Z0)*Y0 )
yw = 2. * Y0  # wings width--2.---1.9
ww = default.ww  # Wall width
gauss_x = xk + X0 + sigmaG * infgcx  # initial gauss center coordinate
s_x = gauss_x + sigmaG * infgcx  # mesh rightside  x coord
s_y = 2. * (yw + Y0)  # mesh rightside +y coord

bay_slope = -Z0 / X0  # Slope of bay at x=0 line. Acc. to eq.

xi = xk
yi = s_y / 2.

# ------------------------------------------------------------------------------
# Check
# ------------------------------------------------------------------------------

if p_stepCheckSW:
    if X0 / Y0 < 10:
        raise Exception("X0/Y0 must bigger than 10")
    if X0 / sigmaG < sqrt(2):
        raise Exception("X0/sigma must equal or bigger than sqrt(2)")
    if sigmaG / -Z0 < 20:
        raise Exception("Wavelength(sigma) must bigger than 20*Z0 ( S-W conditions)")
    if -Z0 / eta0 < 20:
        raise Exception("Amplitude must equal or smaller than  depth/20 (S-W conditions)")


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def gaussian(x, xmean, sigma, amplitude):
    """
    Gaussian function. Generates y value for x
    Parameters: sigma     : Half width of Gauss
                xmean     : center point
                amplitude : amplitude. height at xmean
    """
    return amplitude * util.np.exp(-(x - xmean) ** 2. / (2. * (sigma ** 2.)))


def initial_stage(x, y):
    z = 0 * x
    n = len(x)  # number of points at x
    for i in range(n):
        if gauss_x - sigmaG * 2.5 + ww <= x[i] <= gauss_x + sigmaG * 2.5 - ww:
            z[i] = gaussian(x[i], gauss_x, sigmaG, eta0)
        else:
            z[i] = 0
    return z


def topography(x, y):
    z = 0 * x
    n = len(x)
    for i in range(n):
        if x[i] < X0 + xk:
            z[i] = -bay_slope * ((x[i] - xi) - (X0 / Y0 ** 2.) * (y[i] - yi) ** 2.)
        else:
            z[i] = Z0
    #  if  X0 + xk <= x[i] <= X0 + xk +ww*3. and  (y[i]>=s_y / 2 + yk or y[i]<=s_y / 2 - yk) :
    #      z[i]= 1500.
    return z


# ------------------------------------------------------------------------------
# Setup domain
# ------------------------------------------------------------------------------

p0 = [0, s_y / 2 - yk]
# p1 = [X0 + xk - ww, s_y / 2 - yk]
p1 = [X0 + xk, s_y / 2 - yk]
p4 = [X0 + xk, 0]
p5 = [s_x, 0]
p6 = [s_x, s_y]
p7 = [X0 + xk, s_y]
# p10 = [X0 + xk - ww, s_y / 2 + yk]
p10 = [X0 + xk, s_y / 2 + yk]
p11 = [0, s_y / 2 + yk]
p12 = [xk, s_y / 2 - yk]
p13 = [xk, s_y / 2 + yk]
p14 = [X0 + xk + 10 * ww, 0]
p15 = [X0 + xk + 10 * ww, s_y]

bounding_polygon = [p0, p1, p4, p5, p6, p7, p10, p11]
boundary_tags = {'bottom': [0, 2],
                 'right': [3],
                 'top': [4, 6],
                 'left': [7],
                 'refWall': [1, 5]}

poly1 = [p0, p12, p13, p11]
poly2 = [p12, p1, p4, p14, p15, p7, p10, p13]
poly3 = [p14, p5, p6, p15]

interiors = [[poly1, max_tri_area_inu], [poly2, max_tri_area_focus], [poly3, max_tri_area_base]]

# ------------------------------------------------------------------------------
###STEP-1 Mesh Creation
# ------------------------------------------------------------------------------

if p_stepMesh:
    dirStruct(domain_name)
    os.chdir(prj_output_dir)
    domain = anuga.create_domain_from_regions(bounding_polygon,
                                              boundary_tags,
                                              maximum_triangle_area=max_tri_area_base,
                                              mesh_filename=prj_output_dir + os.sep + mesh_file,
                                              interior_regions=interiors,
                                              use_cache=p_cache,
                                              verbose=p_verbose)
    domain.set_name(domain_name)
    domain.set_datadir(prj_output_dir)
    domain.set_flow_algorithm(flow_alg)
    domain.set_minimum_storable_height(p_minimum_storable_height)
    domain.set_quantity('elevation', topography, location='vertices')
    domain.set_quantity('friction', p_friction)
    domain.set_quantity('stage', initial_stage)

    # ------------------------------------------------------------------------------
    Br = anuga.Reflective_boundary(domain)  # Solid reflective wall
    Bt = anuga.Transmissive_boundary(domain)  # Continue all values on boundary
    Bs = anuga.Transmissive_stage_zero_momentum_boundary(domain)
    Bd = anuga.Dirichlet_boundary([-0.2, 0., 0.])  # Constant boundary values
    # Bw = anuga.Time_boundary(domain=domain,  # Time dependent boundary
    #                         function=lambda t: [(sin(t * pi / 20) - 0.2), 0.0, 0.0])
    # Bw = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(  domain=domain,
    #    function=lambda t: [ 2, 0, 0])
    Bw = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
        domain=domain,
        function=lambda t: [(0 < t < 3660) * 10, 0, 0])

    domain.set_boundary({'left': Bt, 'right': Bt, 'top': Bt, 'bottom': Bt, 'refWall': Br})
    # domain.set_quantities_to_be_monitored('stage-elevation'', polygon=None, time_interval = None)

    # Log section. Record important parameters

    log1 = open(prj_output_dir + os.sep + log1_file, "w+")

    util.printlog("Notes :" + run_note, log1)
    util.printlog("Main parameters ;", log1)
    util.printlog("     Bay length (X0)        :" + str(X0), log1)
    util.printlog("     Bay width (YO)         :" + str(Y0), log1)
    util.printlog("     Open sea depth(Z0)     :" + str(Z0), log1)
    util.printlog("     Gauss wave sigma       :" + str(sigmaG), log1)
    util.printlog("     Gauss amplitude(A/2)   : " + str(eta0 / 2.), log1)
    util.printlog("     Bay toe (x,y_up,y_down): " + str(xi + X0) + "," + str(yi + Y0) + "," + str(yi - Y0), log1)
    util.printlog("Boundaries ;", log1)
    for b in domain.boundary_map:
        util.printlog('   ' + b + '     :   ' + domain.boundary_map[b].__repr__(), log1)
    util.printlog("Timestep           :" + str(timestep1) + "-->" + str(timestep2), log1)
    util.printlog("Simulation time    :" + str(duration1) + "-->" + str(duration2), log1)
    util.printlog("Algorithm          :" + flow_alg, log1)
    util.printlog("Mesh coords;", log1)
    util.printlog(bounding_polygon.__str__(), log1)
    util.printlog("Bay start coords      :" + " (" + str(xi) + " , " + str(yi) + ")", log1)
    util.printlog("Gauss initial x coord :" + str(gauss_x), log1)
    util.printlog("X0/Z0  :" + str(abs(Z0) / X0), log1)
    util.printlog(" Mesh ;", log1)
    util.printlog("   length x width     :" + str(s_x) + " x " + str(s_y), log1)
    util.printlog("   area (km2)         :" + str(anuga.polygon_area(bounding_polygon) / 1000000.0), log1)
    util.printlog("   max triangle base  :" + str(max_tri_area_base), log1)
    util.printlog("   max triangle focus :" + str(max_tri_area_focus), log1)
    util.printlog("   max triangle inund :" + str(max_tri_area_inu), log1)
    util.printlog("   triang. cnt        :" + str(len(domain)), log1)
    util.printlog("   extends            :" + str(domain.get_extent()), log1)
    util.printlog(domain.statistics(), log1)
    log1.close()

    log2 = open(prj_output_dir + os.sep + log2_file, "w+")
    util.printlog(
        "X0;YO;Z0;sigma;eta0;max_tri_area_base;max_tri_area_focus;max_tri_area_inu;timestep2;flow_alg;Area;s_x;s_y;TriangCnt",
        log2)
    util.printlog(str(X0) + ";" + str(Y0) + ";" + str(-Z0) + ";" + str(sigmaG) + ";" + str(eta0 / 2.) + ";" +
                  str(max_tri_area_base) + ";" + str(max_tri_area_focus) + ";" + str(max_tri_area_inu) + ";" + str(
        timestep2) + ";" +
                  str(flow_alg) + ";" + str(anuga.polygon_area(bounding_polygon) / 1000000.0) + ";" +
                  str(s_x) + ";" + str(s_y) + ";" + str(str(len(domain))), log2)

    log2.close()

# ------------------------------------------------------------------------------
###STEP-2  Simulation
# ------------------------------------------------------------------------------

if p_stepSWW:
    os.chdir(prj_output_dir)
    tIndice, tCentroide, tVertex = util.triangCoversPoint(domain.edgelengths, domain.centroid_coordinates,
                                                          domain.vertex_coordinates, xi, yi)
    log1 = open(prj_output_dir + os.sep + log1_file, "a")
    t1 = time.time()

    for t in domain.evolve(yieldstep=timestep1, finaltime=duration1):
        util.printlog(domain.timestepping_statistics(), log1)
    for t in domain.evolve(yieldstep=timestep2, finaltime=duration2):
        # print  str(time.time() - t1)  # her bir adımın çalışma süresi
        t1 = time.time()
        util.printlog(domain.timestepping_statistics(), log1)
    log1.close()
# ------------------------------------------------------------------------------
###STEP-3 Extract Summary
# ------------------------------------------------------------------------------

if p_stepReport:
    os.chdir(prj_output_dir)
    splotter = animate.SWW_plotter(prj_output_dir + os.sep + sww_file)

    # max inund
    max_runup = np.where(splotter.depth > 0, splotter.stage, 0).max()
    idx_max_runup = np.unravel_index(np.where(splotter.depth > 0, splotter.stage, 0).argmax(), splotter.stage.shape)
    max_r_x = splotter.xc[idx_max_runup[1]]
    max_r_y = splotter.yc[idx_max_runup[1]]
    max_r_time = splotter.time[idx_max_runup[0]]

    log3 = open(prj_output_dir + os.sep + log3_file, "w+")
    util.printlog("max runup      :" + str(max_runup), log3)
    util.printlog("max runup loc  :" + str(max_r_x) + " , " + str(max_r_y), log3)
    util.printlog("max runup time :" + str(max_r_time), log3)

    log3.close()
    # Plot Depth and Speed at the max runup time
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(10, 6))

    splotter.triang.set_mask(None)
    c1 = ax1.tripcolor(splotter.triang,
                       facecolors=splotter.depth[idx_max_runup[0], :],
                       cmap='Blues')
    c1bar = plt.colorbar(c1, ax=ax1)
    c1bar.ax.tick_params(labelsize=5)
    ax1.set_title("Depth")

    splotter.triang.set_mask(None)
    c2 = ax2.tripcolor(splotter.triang,
                       facecolors=splotter.speed[idx_max_runup[0], :],
                       cmap='seismic')
    c2bar = plt.colorbar(c2, ax=ax2)
    c2bar.ax.tick_params(labelsize=5)
    ax2.set_title("Speed")

    splotter.triang.set_mask(None)
    c3 = ax3.tripcolor(splotter.triang,
                       facecolors=np.where(splotter.depth[idx_max_runup[0], :] > 0, splotter.stage[idx_max_runup[0], :],
                                           splotter.depth[idx_max_runup[0], :]),
                       cmap='Blues')
    c3bar = plt.colorbar(c3, ax=ax3)
    c3bar.ax.tick_params(labelsize=5)
    ax3.set_title("Runup")
    fig.tight_layout()
    fig.savefig(prj_output_dir + os.sep + plot_file)

    plt.close()

# ------------------------------------------------------------------------------
###STEP-4 Generate Master Gauge File
# ------------------------------------------------------------------------------

if p_stepGauge:
    os.chdir(prj_output_dir)
    import csv

    gaugefile = open(prj_output_dir + os.sep + gauge_file, 'wb+')
    gauge_writer = csv.writer(gaugefile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    maxr_elv = bay_slope * ((max_r_x - xi) - (X0 / Y0 ** 2.) * (max_r_y - yi) ** 2.)

    gauge_writer.writerow(['easting', 'northing', 'name', 'elevation'])
    gauge_writer.writerow([max_r_x, max_r_y, 'maxrunup', 0.])  # maxr_elv])
    gauge_writer.writerow([xi, yi, 'bay00', 0.])
    gauge_writer.writerow([xi + X0, yi, 'baymouth', 0.])
    gauge_writer.writerow([xi + X0 - Y0, yi, 'bay_i01', 0.])
    gauge_writer.writerow([xi + X0 + Y0, yi, 'bay_i02', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi, 'm1', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.2, 'm2', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.4, 'm3', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.6, 'm4', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.8, 'm5', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.9, 'm6', 0.])
    gauge_writer.writerow([xi + X0 + sigmaG, yi * 1.98, 'm7', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi, 'n1', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.2, 'n2', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.4, 'n3', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.6, 'n4', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.8, 'n5', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.9, 'n6', 0.])
    gauge_writer.writerow([xi + X0 + 2.6 * sigmaG, yi * 1.98, 'n7', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi, 'p1', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.2, 'p2', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.4, 'p3', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.6, 'p4', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.8, 'p5', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.9, 'p6', 0.])
    gauge_writer.writerow([xi + X0 + 0.4 * sigmaG, yi * 1.98, 'p7', 0.])

    gaugefile.close()
# ------------------------------------------------------------------------------
###STEP-5 Generate Gauge Files
# ------------------------------------------------------------------------------

if p_stepGaugeChart:
    os.chdir(prj_output_dir)
    anuga.sww2csv_gauges(prj_output_dir + os.sep + sww_file,
                         prj_output_dir + os.sep + gauge_file,
                         domain_name + '_gauge_',
                         quantities=['stage', 'speed', 'elevation'],
                         verbose=True)

    import pylab

    anuga.csv2timeseries_graphs(directories_dic={prj_output_dir + os.sep: [r'', 0, 0]},
                                output_dir=prj_output_dir + os.sep + 'gauges' + os.sep,
                                base_name=domain_name + '_gauge_',
                                plot_numbers='',
                                quantities=['stage', 'speed', 'elevation'],
                                extra_plot_name='',
                                assess_all_csv_files=True,
                                create_latex=False,
                                verbose=True)

print "Executions:"
print "S-W Conditions Check   :" + p_stepCheckSW.__str__()
print "Mesh generation        :" + p_stepMesh.__str__()
print "Simulation             :" + p_stepSWW.__str__()
print "Report                 :" + p_stepReport.__str__()
print "Generate master gauge  :" + p_stepGauge.__str__()
print "Generate gauges+charts :" + p_stepGaugeChart.__str__()
print "Finished!"
