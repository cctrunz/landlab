from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.components import PresFlowNetwork, MeltCreep
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


nsteps=1000
nx=50 #number of nodes
dt=500.
L=10000.
Zmax=500. # linear change in ice thickness between max and min
Zmin = 1
Zslope = 0. # or can set a slope
Qpeak = 2.
Qbase = 0.
Qsteady=False
period = 24.*60.*60. 
bedslope =0.
A_R = np.pi*1.5**2 # moulin cross-section area
D0 = 1 # initial hydraulic diameter of the conduit
hin=400. # initial upstream head
hout=0. # intial downstrem head
every=10 # how frequently to record data or make plots
hymin=0. # plot limits
hymax=1000 # plot limits
dymin=0 # plot limits
dymax=3 # plot limits

#figure out what the node spacing is
dx = L/(nx-1)
mg = RasterModelGrid((3,nx),dx) #3 cells wide because no flow boundary on the edges
junc_elev = mg.add_zeros('node', 'junction__elevation') #create elevation for the nodes
mg.at_node['junction__elevation'] = bedslope*(mg.node_x - L) # add values in the nodes for elevation based on baseslope and distance from the edge
R = mg.add_zeros('node', 'input__discharge') # possibility for recharge for all nodes
#d_h = mg.add_ones('link', 'hydraulic__diameter')
d_h = mg.add_zeros('link','hydraulic__diameter') # 4x.. look up formula semi-circular
mg.at_link['hydraulic__diameter'][mg.active_links]= D0 

h = mg.add_zeros('node', 'hydraulic__head') #initialize head at 0
#choose to use slope of thickness
if Zmin != None:
    Zslope = (Zmax - Zmin)/L
thickness = Zmax*np.ones(mg.number_of_nodes) - Zslope*mg.node_x #initialize ice thickness at 0
Z = mg.add_field('node', 'ice__thickness', thickness) #create variables and set them to values of thickness
Q = mg.add_zeros('link', 'conduit__discharge') #create variable and set to zero

#Boundary conditions
#set heads at edges
h[mg.nodes_at_left_edge] = hin
h[mg.nodes_at_right_edge] = hout
h[mg.nodes_at_top_edge] = 0.
h[mg.nodes_at_bottom_edge] = 0.
mg.set_closed_boundaries_at_grid_edges(False,True,False,True) # for the edges
Q[mg.active_links] = 1. # set all the discharge to 1 (initialization)
#
pfn = PresFlowNetwork(mg) # pressurized flow network conduit -- python object that contains the ... for flow
pfn.run_one_step() # solves for the flow -- 
mc = MeltCreep(mg, dt=dt) # solves for melt and creep -- output hydraulic diameter

#type of QIN
def moulin_recharge(t):
    if not Qsteady:
        return Qpeak*(1. + np.sin(2*np.pi*t/period))/2. + Qbase
    else:
        return Qpeak

#PLOTS
fig1, axs = plt.subplots(3,1,figsize=(6,10))
time = 0.
time_list = []
h_list = []
d_list = []
q_list = [] #discharge
r_list = []


for step in np.arange(nsteps):
    pfn.run_one_step()
    Qnow = Q[mg.active_links][0]
    dh_dt = (moulin_recharge(time) - Qnow)/A_R
    h[mg.nodes_at_left_edge] = h[mg.nodes_at_left_edge] + dh_dt*dt
    q_list.append(Qnow)
    r_list.append(moulin_recharge(time))
    time_list.append(time)
    h_list.append(mg.at_node['hydraulic__head'][mg.core_nodes])#[nx:2*nx])
    d_list.append(mg.at_link['hydraulic__diameter'][mg.active_links])
    mc.run_one_step() #calculate new hydraulic diameter
    time += dt #update time
    # if (step % every)==0: #make an animation frame -- if remainder is 0 create a plot
    #     print("step =",step, " avg d_h=",mg.at_link['hydraulic__diameter'].mean())
    #     plt.subplot(3,1,1)
    #     #imshow_grid_at_node(mg, h)
    #     plt.plot(mg.at_node['hydraulic__head'][mg.core_nodes])
    #     plt.ylabel('Hydraulic head')
    #     plt.ylim([hymin,hymax])
    #     plt.subplot(3,1,2)
    #     plt.plot(mg.at_link['conduit__discharge'][mg.active_links])
    #     plt.ylabel('Discharge')
    #     plt.ylim([0,Qpeak*1.25])
    #     plt.subplot(3,1,3)
    #     plt.plot(mg.at_link['hydraulic__diameter'][mg.active_links])
    #     plt.ylabel('Diameter')
    #     plt.ylim([dymin,dymax])
    #     image_name = 'heads_and_discharge'+str(step).zfill(6)+'.png'
    #     plt.tight_layout()
    #     plt.savefig(image_name)
    #     fig1.clf()
h_arr = np.array(h_list) # head[time,position]
q_arr = np.array(q_list) # discharge[time,position]
r_arr = np.array(r_list) # recharge[time,position]
d_arr = np.array(d_list) # subglacial channel diameter[time,position]
time_arr = np.array(time_list)

results = {'h':h_arr, 'q':q_arr, 'r':r_arr,'d':d_arr,'t':time_arr}


#%%to plot 1 position for all timesteps:
time_day = time_arr/3600/24

position = 0
hx = h_arr[:,position]
plt.figure()
plt.plot(time_list,hx)

#to plot 1 timestep for all position:
x = np.linspace(dx/2,nx*dx,nx-1)
timestep = 500
diameter = d_arr[timestep,:]
plt.figure()
plt.plot(x,diameter)

#%%to plot multiple timesteps for all position:
plt.figure()
x = np.linspace(dx/2,nx*dx,nx-1)
iterations = np.arange(0,nsteps,10)
color=cm.rainbow(np.linspace(0,1,len(iterations)))
for timestep,c in zip(iterations,color): 
    diameter = d_arr[timestep,:]
    plt.plot(x,diameter,c=c)
    
#%%
plt.figure()
x = np.linspace(dx,nx*dx,nx-2)
iterations = np.arange(0,nsteps,10)
color=cm.rainbow(np.linspace(0,1,len(iterations)))
plt.plot(x,mg.at_node['ice__thickness'][mg.core_nodes],color='black',linewidth=3)
for timestep,c in zip(iterations,color): 
    head = h_arr[timestep,:]
    plt.plot(x,head,c=c)
    
#%%to plot multiple timesteps for all position:
plt.figure()
plt.plot(time_day,r_arr,color='black',label='Qin')
plt.plot(time_day,q_arr,color='blue',label='Qout')
plt.legend()
























