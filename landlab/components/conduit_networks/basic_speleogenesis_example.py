from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.components import PresFlowNetwork, MeltCreep
import numpy as np
from matplotlib import colors,  cm
import matplotlib.animation as animation
from landlab.components.conduit_networks.calcite import calcite_diss_palmer_transfer
from pylab import *


from PIL import Image
from PIL import ImageDraw

def plot_links(grid, value_name, autoscale=True,
               vmin=0., vmax=0., cmap_name='viridis',
               magnitude= False, lw=5,
               pixels_per_node=31, x_name='X',
               y_name='Y', var_name='', use_PIL=True):#, logvalues=False ):
    if use_PIL:
        pixel_xscaling_factor = pixels_per_node/grid.dx
        pixel_yscaling_factor = pixels_per_node/grid.dy
        im = Image.new('RGBA', ( int(round(grid.number_of_cell_columns*pixels_per_node)), int(round(grid.number_of_cell_rows*pixels_per_node)) ), (255,255,255,255))
        draw = ImageDraw.Draw(im)

    link_head_x = grid.node_x[grid.node_at_link_head]
    link_head_y = grid.node_y[grid.node_at_link_head]
    link_tail_x = grid.node_x[grid.node_at_link_tail]
    link_tail_y = grid.node_y[grid.node_at_link_tail]
    if magnitude:
        values = abs(grid.at_link[value_name])
    else:
        values = grid.at_link[value_name]
    #if logvalues:
    #    values = np.log10(abs(grid.at_link[value_name]))


    #Normalize color values
    if autoscale:
        cnorm = colors.Normalize()
        cnorm.autoscale(values)
    else:
        cnorm = colors.Normalize(vmin,vmax)
    scalarMap = cm.ScalarMappable(norm=cnorm, cmap = get_cmap(cmap_name))
    scalarMap.set_array(values)
    #set_cmap()
    if use_PIL:
        for i, value in enumerate(values):
            draw.line( ((link_head_x[i]*pixel_xscaling_factor,link_head_y[i]*pixel_yscaling_factor),(link_tail_x[i]*pixel_xscaling_factor,link_tail_y[i]*pixel_yscaling_factor)),fill=scalarMap.to_rgba(value,bytes=True),width=lw)
        imshow(np.asarray(im), origin='lower', extent=(0,grid.number_of_cell_columns,0,grid.number_of_cell_rows))
    else:
        for i, value in enumerate(values):
            xs = [link_head_x[i],link_tail_x[i]]
            ys = [link_head_y[i],link_tail_y[i]]
            img = plot(xs,ys,lw=lw,color=scalarMap.to_rgba(value))
    cb = colorbar(scalarMap)
    cb.ax.set_ylabel(var_name)
    xlabel(x_name)
    ylabel(y_name)

mg = RasterModelGrid((25,25),500)
junc_elev = mg.add_zeros('node', 'junction__elevation')
R = mg.add_zeros('node', 'input__discharge')
h = mg.add_zeros('node', 'hydraulic__head')
Q = mg.add_zeros('link', 'conduit__discharge')

#set heads at edges
h[mg.nodes_at_left_edge] = 10.
h[mg.nodes_at_right_edge] = 0.
h[mg.nodes_at_top_edge] = 0.
h[mg.nodes_at_bottom_edge] = 0.

mg.set_closed_boundaries_at_grid_edges(False,True,False,True)
Q[mg.active_links] = 0.1*np.random.rand(mg.number_of_active_links)
n_core = mg.number_of_core_nodes
links = mg.links_at_node
print("Number of links = ", mg.number_of_links)
print("Number of nodes = ", mg.number_of_nodes)
print("Number of active links = ", mg.number_of_active_links)
print("Number of core nodes = ", mg.number_of_core_nodes)

d_h = mg.add_zeros('link','hydraulic__diameter')
mg.at_link['hydraulic__diameter'][mg.active_links]= 0.01*np.random.rand(mg.number_of_active_links)

rates = mg.add_zeros('link','diss__rates')

#code added for fsolve algorithm
#dhdx = mg.add_zeros('link', 'hydraulic_head_gradient')
#net_node_flux = mg.add_ones('node', 'net_node_flux')


pfn = PresFlowNetwork(mg, solutes=['Ca', 'PCO2'],
                      transfer_func=calcite_diss_palmer_transfer,
                     transfer_func_link_variables = ['hydraulic__diameter','conduit__discharge'],
                     transfer_kwd_args = {'T_C':10.})
#Here need to set boundary conditions for transport calculation
#add input to node 310
Ca = mg.at_node['concentration__Ca']
PCO2 = mg.at_node['concentration__PCO2']

input_Q=.005
input_C = 0.001
input_idx = 30+250#11
mg.at_node['input__discharge'][input_idx] = input_Q
mg.at_node['Ca__inflow_conc'][input_idx] = input_C
mg.at_node['PCO2__inflow_conc'][input_idx] = 0.05

Ca[mg.nodes_at_left_edge] = 0.0005
Ca[mg.nodes_at_right_edge] = 0.
Ca[mg.nodes_at_top_edge] = 0.
Ca[mg.nodes_at_bottom_edge] = 0.

PCO2[mg.nodes_at_left_edge] = 0.01
PCO2[mg.nodes_at_right_edge] = 0.
PCO2[mg.nodes_at_top_edge] = 0.
PCO2[mg.nodes_at_bottom_edge] = 0.

nsteps = 10000
step_start=0
dt = 1. #yr
every = 10
fig1, axs = subplots(2,3,figsize=(9,10))
for step in arange(nsteps)+step_start:
    print("Step: ",step)
    pfn.run_one_step(transport=True)
    #Calculate change of conduit diameters
    ddh_dt = mg.at_link['diss__rates']*0.001*dt
    mg.at_link['hydraulic__diameter'] += ddh_dt#[mg.active_links]
    print("max ddh_dt = ",str(max(ddh_dt)))
    print("min ddh_dt = ",str(min(ddh_dt)))
    print("max d_h=",str(max(d_h)))
    #Reset boundary node values
    Ca[mg.nodes_at_left_edge] = 0.0005
    Ca[mg.nodes_at_right_edge] = 0.
    Ca[mg.nodes_at_top_edge] = 0.
    Ca[mg.nodes_at_bottom_edge] = 0.

    PCO2[mg.nodes_at_left_edge] = 0.01
    PCO2[mg.nodes_at_right_edge] = 0.
    PCO2[mg.nodes_at_top_edge] = 0.
    PCO2[mg.nodes_at_bottom_edge] = 0.

    if (step % every)==0: #make an animation frame
        subplot(3,2,1)
        imshow_grid_at_node(mg, h, var_name='Head ($m$)', cmap='viridis')
        subplot(3,2,2)
        plot_links(mg, 'conduit__discharge', magnitude=True, var_name='Discharge ($m^3/s$)', lw=9,use_PIL=True)
        subplot(3,2,3)
        plot_links(mg, 'hydraulic__diameter', var_name='Hydraulic diameter ($m$)', lw=9,use_PIL=True)
        subplot(3,2,4)#, position=[0.2,0.2,0.9,0.9])
        plot_links(mg, 'diss__rates', var_name='Dissolution rate ($mm/yr$)', lw=9,use_PIL=True)
        subplot(3,2,5)
        imshow_grid_at_node(mg, Ca, var_name='Calcium ($mol/L$)', cmap='viridis')
        subplot(3,2,6)
        imshow_grid_at_node(mg, PCO2, var_name='PCO2 (atm)', cmap='viridis')
        image_name = 'speleo_example/speleogenesis_output'+str(step).zfill(6)+'.png'
        subplots_adjust(wspace=0.3)
        savefig(image_name)
        fig1.clf()
