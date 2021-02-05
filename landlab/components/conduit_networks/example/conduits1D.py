from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.components import PresFlowNetwork, MeltCreep
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def run1Dsim(nsteps=1000,
             nx=50, #number of nodes
             dt=500.,
             L=10000.,
             #choose type of ice thickness profile: 
             #'constant' = same value everywhere. require Zmax value;
             #'linear' = decrease linearly from Zmax to Zmin or from Zmax with Zslope
             #'sqrt' = decrease with a sqrt function. Only requires L
             #'custom' = doesn't exist yet, but will be used to input a custom profile
             Ztype = 'sqrt', 
             Zmax = None, # linear change in ice thickness between max and min
             Zmin = None,
             Zslope = None, # or can set a slope
             Qpeak = 2.,
             Qbase = 0.,
             Qsteady=False,
             period = 24.*60.*60., 
             bedslope =0.,
             A_R = 50., # moulin cross-section area
             D0 = 0.5, # initial hydraulic diameter of the conduit
             hin=400., # initial upstream head
             hout=0., # intial downstrem head
             # every=10, # how frequently to record data or make plots
             # hymin=0., # plot limits
             # hymax=1000., # plot limits
             # dymin=0., # plot limits
             # dymax=3.): # plot limits
             ):
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

    #Ice profile

    if Ztype == 'constant':
    	thickness = Zmax*np.ones(mg.number_of_nodes)
    if Ztype == 'linear':
        if Zmin != None:
        	Zslope = (Zmax - Zmin)/L
        else:
        	Zslope = 0
        thickness = Zmax*np.ones(mg.number_of_nodes) - Zslope*mg.node_x #initialize ice thickness at 0
    if Ztype == 'sqrt':
        H0= 1500./np.sqrt(60000.) #assures ice thickness of 1500 m at 60 km from edge
        thickness = np.flip(H0*np.sqrt(mg.node_x))

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
    #fig1, axs = plt.subplots(3,1,figsize=(6,10))
    time = 0.
    time_list = []
    h_list = []
    d_list = []
    q_list = [] #discharge
    r_list = []
    x_node = np.linspace(dx,nx*dx,nx-2)/1000 #distance in km
    x_link = np.linspace(dx,nx*dx,nx-1)/1000 #distance in km


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
    h_arr = np.array(h_list)
    q_arr = np.array(q_list)
    r_arr = np.array(r_list)
    d_arr = np.array(d_list)
    time_arr = np.array(time_list)
    results = {'h':h_arr, 'q':q_arr, 'r':r_arr,'d':d_arr,'t':time_arr,\
               'ice_thickness':mg.at_node['ice__thickness'][mg.core_nodes],\
                'x_node':x_node, 'x_link':x_link}
    return results

def plot_3panels_singlestep(results,timestep):
    '''
    Plot figure with 3 panels. 
   - upper panel is Qin and Qout over time; 
   - middle panel is head and ice thickness profile; 
   - lower panel is the subglacial channel hydraulic diameter.

    Parameters
    ----------
    results : dictionnary
        output of run1Dsim.
    timestep : TYPE
        DESCRIPTION.

    Returns
    -------
    figure

    '''
    
    Qin = results['r']
    Qout = results['q']
    head = results['h']
    diameter = results['d']
    time_day = results['t']/3600/24
    ice_thickness = results['ice_thickness']  
    x_node = results['x_node']
    x_link = results['x_link']

    fig = plt.figure(figsize=(5,6))
    
    plt.subplot(3,1,1)
    plt.plot(time_day[0:timestep],Qin[0:timestep],color='blue',label='Qin')
    plt.plot(time_day[0:timestep],Qout[0:timestep],color='red',label='Qout')
    plt.xlim([time_day[0],time_day[-1]])
    plt.ylim([np.min(Qout),np.max(Qout)])
    plt.ylabel('($m^3/s$)')
    plt.xlabel('(days)')
    plt.legend()

    plt.subplot(3,1,2)
    plt.plot(x_node,ice_thickness,color='black',linewidth=3)
    plt.plot(x_node,head[timestep,:])
    plt.ylim([0,ice_thickness[0]+ice_thickness[0]/2])
    plt.ylabel('head (m)')
    plt.xlabel('distance from moulin (km)')
    
    plt.subplot(3,1,3)
    plt.plot(x_link,diameter[timestep,:],color='black',linewidth=1.5)
    plt.ylim([np.min(diameter),np.max(diameter)])
    plt.ylabel('sub channel HD (m)')
    plt.xlabel('distance from moulin (km)')
    
    fig.tight_layout()
    


def plot_3panels_movie(results,path,every=10): # how frequently to record data or make plots
    '''
    
    Parameters
    ----------
    results : dictionnary
        output of run1Dsim
    
    path : string
        path where to save the figure. end widthout a slash!

    Returns
    -------
    save a png of the plot for each timestep

    '''
    Qin = results['r']
    Qout = results['q']
    head = results['h']
    diameter = results['d']
    time_day = results['t']/3600/24
    ice_thickness = results['ice_thickness']  
    x_node = results['x_node']
    x_link = results['x_link']
    idx = 0
    
    for timestep in np.arange(len(time_day)): 
        
        if (timestep % every)==0: #make an animation frame -- if remainder is 0 create a plot
            fig = plt.figure(figsize=(5,6))
            
            plt.subplot(3,1,1)
            plt.plot(time_day[0:timestep],Qin[0:timestep],color='blue',label='Qin')
            plt.plot(time_day[0:timestep],Qout[0:timestep],color='red',label='Qout')
            plt.xlim([time_day[0],time_day[-1]])
            plt.ylim([np.min(Qout),np.max(Qout)])
            plt.ylabel('($m^3/s$)')
            plt.xlabel('(days)')
            plt.legend()
        
            plt.subplot(3,1,2)
            plt.plot(x_node,ice_thickness,color='black',linewidth=3)
            plt.plot(x_node,head[timestep,:])
            plt.ylim([0,ice_thickness[0]+ice_thickness[0]/2])
            plt.ylabel('head (m)')
            plt.xlabel('distance from moulin (km)')
            
            plt.subplot(3,1,3)
            plt.plot(x_link,diameter[timestep,:],color='black',linewidth=1.5)
            plt.ylim([np.min(diameter),np.max(diameter)])
            plt.ylabel('sub channel HD (m)')
            plt.xlabel('distance from moulin (km)')
            
            fig.tight_layout()
            
            plt.savefig('%s/figure%s.png'%(path,idx))
            plt.clf()
            plt.close(fig)
            idx = idx+1
            
def plot_2panel_overtime_multiposition(results):

    head = results['h']
    diameter = results['d']
    time_day = results['t']/3600/24
    ice_thickness = results['ice_thickness']  
    x_link = results['x_link']
     
    iterations = np.arange(0,len(x_link),10)
    color=cm.rainbow(np.linspace(0,1,len(iterations)))
    
    plt.figure(figsize=(8,6))
    plt.subplot(2,1,1)
    for position,c in zip(iterations,color): 
        hx = head[:,position]
        plt.plot(time_day,hx,c=c,label='H=%s'%np.int(ice_thickness[position]))
    plt.ylim([0,ice_thickness[0]])
    plt.ylabel('head (m)')
    plt.xlabel('time (day)')
    plt.legend()
    
    plt.subplot(2,1,2)
    for position,c in zip(iterations,color): 
        dx = diameter[:,position]
        plt.plot(time_day,dx,c=c)
    #plt.ylim([0,ice_thickness[0]])
    plt.ylabel('sub. channel HD (m)')
    plt.xlabel('time (day)')