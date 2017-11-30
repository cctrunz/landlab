from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.components import PresFlowNetwork, MeltCreep
import numpy as np
import matplotlib.pyplot as plt


def run1Dsim(nsteps=1000,
             nx=50,
             dt=500.,
             L=10000.,
             Zmax=500.,
             Zmin = None,
             Zslope =0.,
             Qpeak = 2.,
             Qbase = 0.,
             Qsteady=False,
             period = 24.*60.*60., 
             bedslope =0.,
             A_R = 50.,
             D0 = 0.5,
             hin=400.,
             hout=0.,
             every=10,
             hymin=0.,
             hymax=1000.,
             dymin=0.,
             dymax=3.):

    dx = L/(nx-1)
    mg = RasterModelGrid((3,nx),dx)
    junc_elev = mg.add_zeros('node', 'junction__elevation')
    mg.at_node['junction__elevation'] = bedslope*(mg.node_x - L)
    R = mg.add_zeros('node', 'input__discharge')
    #d_h = mg.add_ones('link', 'hydraulic__diameter')
    d_h = mg.add_zeros('link','hydraulic__diameter')
    mg.at_link['hydraulic__diameter'][mg.active_links]= D0 

    h = mg.add_zeros('node', 'hydraulic__head')
    if Zmin != None:
        Zslope = (Zmax - ZMin)/L
    thickness = Zmax*np.ones(mg.number_of_nodes) - Zslope*mg.node_x
    Z = mg.add_field('node', 'ice__thickness', thickness)
    Q = mg.add_zeros('link', 'conduit__discharge')
    
    #set heads at edges
    h[mg.nodes_at_left_edge] = hin
    h[mg.nodes_at_right_edge] = hout
    h[mg.nodes_at_top_edge] = 0.
    h[mg.nodes_at_bottom_edge] = 0.
    mg.set_closed_boundaries_at_grid_edges(False,True,False,True)
    Q[mg.active_links] = 1.
    pfn = PresFlowNetwork(mg)
    pfn.run_one_step()
    mc = MeltCreep(mg, dt=dt)
    def moulin_recharge(t):
        if not Qsteady:
            return Qpeak*(1. + np.sin(2*np.pi*t/period))/2. + Qbase
        else:
            return Qpeak

    fig1, axs = plt.subplots(3,1,figsize=(6,10))
    time = 0.
    time_list = []
    h_list = []
    d_list = []
    q_list = []
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
        mc.run_one_step()
        time += dt
        if (step % every)==0: #make an animation frame
            print "step =",step, " avg d_h=",mg.at_link['hydraulic__diameter'].mean()
            plt.subplot(3,1,1)
            #imshow_grid_at_node(mg, h)
            plt.plot(mg.at_node['hydraulic__head'][mg.core_nodes])
            plt.ylabel('Hydraulic head')
            plt.ylim([hymin,hymax])
            plt.subplot(3,1,2)
            plt.plot(mg.at_link['conduit__discharge'][mg.active_links])
            plt.ylabel('Discharge')
            plt.ylim([0,Qpeak*1.25])
            plt.subplot(3,1,3)
            plt.plot(mg.at_link['hydraulic__diameter'][mg.active_links])
            plt.ylabel('Diameter')
            plt.ylim([dymin,dymax])
            image_name = 'heads_and_discharge'+str(step).zfill(6)+'.png'
            plt.tight_layout()
            plt.savefig(image_name)
            fig1.clf()
    h_arr = np.array(h_list)
    q_arr = np.array(q_list)
    r_arr = np.array(r_list)
    d_arr = np.array(d_list)
    time_arr = np.array(time_list)
    results = {'h':h_arr, 'q':q_arr, 'r':r_arr,'d':d_arr,'t':time_arr}
    return results
