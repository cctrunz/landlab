from landlab import Component, FieldError
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.optimize import fsolve
import time
#from numba import jit

class PresFlowNetwork(Component):
    """
    Calculate water flow within a pressurized network.

    Landlab component that solves for flow in a 2D pressurized
    network (e.g. karst conduits or subglacial conduits).

    Construction:

        PressurizedFlowNetwork(grid, [stuff to add later])

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid object.

    [add later]

    Notes
    -----
    Solver uses Newton-Raphson Global Algorithm from:

    Todini and Rossmann (2013), Unified Framework for Deriving Siultaneous Equation Algorithms for Water Distribution Networks, Journal of Hydraulic Engineering, 139, 5, 511-526.



    """
    _name = 'PresFlowNetwork'

    _input_var_names = (
        'junction__elevation',
        'hydraulic__diameter',
        'input__discharge',
    )

    _output_var_names = (
        'resistance__factor',
        'hydraulic__head',
        'conduit__discharge',
        'head_loss__exponent',
    )

    _var_units = {
        'junction__elevation':'m',
        'hydraulic__diameter':'m',
        'resistance__factor':'(s/m2)^head_loss__exponent',
        'head_loss__exponent':'unitless',
        'hydraulic__head':'m',
        'conduit__discharge':'m3/s',
        'input__discharge':'m3/s',
    }

    _var_mapping = {
        'junction__elevation': 'node',
        'hydraulic__diameter':'link',
        'resistance__factor':'link',
        'hydraulic__head':'node',
        'conduit__discharge':'link',
        'input__discharge':'node',
    }

    _var_doc = {
        'junction__elevation':
            'elevation of the node junction relative to some datum',
        'hydraulic__diameter':
            'equivalent diameter of the conduit segment (D_H=4A_c/P_w)',
        'resistance__factor':
            'resistance factor (r) used in head loss equation dh = r Q^a',
        'head_loss__exponent':
            'exponent (a) in head loss equation dh = r Q^a',
        'hydraulic__head':
            'elevation head at node',
        'conduit__discharge':
            'flow discharge through conduit in direction of link',
        'input__discharge':
            'flow discharge that is added at nodes (positive is inflow)',
    }

    def __init__(self, grid, flow_eqn='Darcy-Weisbach',
                    shape_factor=np.pi, g=9.81, f=0.1, solutes=None,**kwds):
        """
        Initialize the PresFlowNetwork

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        flow_eqn : str, optional (defaults to 'Darcy-Weisbach')
            Flow equation to be used in conduits.
        """
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        available_flow_eqns = ['Darcy-Weisbach']
        if flow_eqn in available_flow_eqns:
            self.flow_eqn = flow_eqn
        else:
            raise FieldError(
                'Flow equation '+ str(flow_eqn)+' not available!'
            )
        self.g = g
        self.f = f
        self.shape_factor = shape_factor
        ##Create fields
        if 'junction__elevation' in grid.at_node:
            self.junc_elev = grid.at_node['junction__elevation']
        else:
            raise FieldError(
                'Junction elevations are required as a component input!')
        if 'hydraulic__diameter' in grid.at_link:
            self.d_h = grid.at_link['hydraulic__diameter']
        else:
            raise FieldError(
                'Hydraulic diameters of conduits are required as a component input!')
        if 'input__discharge' in grid.at_node:
            self.input__discharge = grid.at_node['input__discharge']
        else:
            raise FieldError(
                'Node input discharges are required as a component input!')
        if 'resistance__factor' in grid.at_link:
            self.r = grid.at_link['resistance__factor']
        else:
            self.r = grid.add_ones('link', 'resistance__discharge')
        if 'head_loss__exponent' in grid.at_link:
            self.a = grid.at_link['head_loss__exponent']
        else:
            self.a = grid.add_ones('link', 'head_loss__exponent')
        if 'hydraulic__head' in grid.at_node:
            self.h = grid.at_node['hydraulic__head']
        else:
            self.h = grid.add_zeros('node', 'hyraulic__head')
        if 'conduit__discharge' in grid.at_link:
            self.Q = grid.at_link['conduit__discharge']
        else:
            self.Q = grid.add_ones('link', 'conduit__discharge')
        self.solutes = solutes
        if self.solutes != None:
            #Loop through list of solutes and create associated link and node values
            for solute in solutes:
                if type(solute)!=str:
                    raise FieldError(
                        'Solutes keyword argument must be a list of strings of solute names!')
                #create empty node solute values
                self.grid.add_zeros('node','concentration__'+solute)
                #create empty link solute concentration values
                self.grid.add_zeros('link',solute+'__conc_in')
                self.grid.add_zeros('link',solute+'__conc_out')
                #Create node inflow concentration on grid if not already available
                if not solute+'__inflow_conc' in grid.at_node:
                    self.grid.add_zeros('node', solute+'__inflow_conc')

    def calc_r(self):
        """
        Calculate flow resistance factors for conduits (links)
        """
        L = self.grid.length_of_link
        f = self.f
        g = self.g
        d_h = self.d_h
        r_s = self.shape_factor
        a = self.a
        if self.flow_eqn == 'Darcy-Weisbach':
            self.r = 8*f*L/(g*r_s**2.*d_h**5.)

    def calc_a(self):
        """
        Calculate and set head loss exponent.
        """
        if self.flow_eqn == 'Darcy-Weisbach':
            self.a = 2.0*np.ones(self.grid.number_of_links)

######################
#### Begin attempt to incorporate direct residual solver (not working yet)
##############

    def network_residuals(self,heads):
        a = self.a
        r = self.r
        self.h[self.grid.core_nodes] = heads
        dh = self.grid.length_of_link * self.grid.calc_grad_at_link(self.h)
        self.Q = np.sign(dh)*(np.fabs(dh)/r)**(1./a)
        return self.grid.calc_net_flux_at_node(self.Q)[self.grid.core_nodes]/self.grid.dx - self.grid.at_node['input__discharge'][self.grid.core_nodes]


    def run_one_step_fsolve(self, **kwds):
        self.calc_a()
        self.calc_r()
        self.h[self.grid.core_nodes] = fsolve(self.network_residuals, self.h[self.grid.core_nodes])

######################
#### End attempt to incorporate direct residual solver (not working yet)
##############

#    @jit(nopython=True)#This doesn't work directly. Probably have to create a jitclass

    def run_one_step(self, transport=False):#, **kwds):
        self.calculate_flow()
        if transport:
            self.steady_state_transport()

    def calculate_flow(self):
        #Calculate flow in network
        max_tol = 1e-5#0.001
        tol = 1.
        niter = 1
        links = self.grid.links_at_node
        n_core = self.grid.number_of_core_nodes
        self.calc_a()
        a = self.a
        self.calc_r()
        r = self.r
        while tol>max_tol:
            start = time.time()
            Q_ij = self.Q[links][self.grid.core_nodes]*self.grid.active_link_dirs_at_node[self.grid.core_nodes]
            F = np.sum(Q_ij*(1. - 1./a[links][self.grid.core_nodes] ), axis=1)
            #Add recharge to nodes
            F += self.grid.at_node['input__discharge'][self.grid.core_nodes]
            ADA_ij = np.zeros([n_core, n_core])
            for i, this_node in enumerate(self.grid.core_nodes):
                #calculate diagonal terms in the matrix A_{21} (D_{11})^{-1} A_{12}
                node_Qs = Q_ij[i]
                node_as = a[links][self.grid.core_nodes][i]
                node_rs = r[links][self.grid.core_nodes][i]
                ADA_ii = 0.
                for link_idx, link_Q in enumerate(node_Qs):
                    if link_Q != 0:
                        ADA_ii += 1./(node_as[link_idx]*node_rs[link_idx]*np.fabs(link_Q)**(node_as[link_idx]-1.))
                ADA_ij[i][i] = ADA_ii
                #loop through node links to calculate off-diagonal terms in A
                for this_link in links[this_node]:
                    if this_link >= 0:#Ignore -1 cases, which indicate no link
                        neighbor_node = self.grid.node_at_link_head[this_link]
                        #Check whether we have the neighbor node
                        if neighbor_node == this_node:
                            #neighbor is at tail, not head
                            neighbor_node = self.grid.node_at_link_tail[this_link]
                        #Make sure this neighbor is not a boundary node (these are treated differently)
                        if not self.grid.node_is_boundary(neighbor_node):
                            #Get proper j index for Array (neighbor core node number)
                            j = np.where(self.grid.core_nodes==neighbor_node)[0][0]
                            ADA_ij[i][j] = -1./(a[this_link]*r[this_link]*np.fabs(self.Q[this_link])**(a[this_link]-1))
                        else:#if boundary add term to F
                            if self.Q[this_link]!=0:#Check whether this is a closed boundary
                                F[i] += 1./(a[this_link]*r[this_link]*np.fabs(self.Q[this_link])**(a[this_link]-1))*self.h[neighbor_node]
            #Turn A into sparse matrix
            A_csr = csr_matrix(ADA_ij)
            ##Solve linear system for new approximation of heads
            # ADA_ij H = F
            self.h[self.grid.core_nodes] = spsolve(A_csr,F)
            #Evaluate scalar equations to update discharge
            dQ= -(1./a[self.grid.active_links])*self.Q[self.grid.active_links] - 1./(a[self.grid.active_links]*r[self.grid.active_links]*np.fabs(self.Q[self.grid.active_links])**(a[self.grid.active_links]-1)) * (self.grid.calc_diff_at_link(self.h)[self.grid.active_links])
            self.Q[self.grid.active_links] +=  dQ
            #tol = sum(np.fabs(dQ))/sum(np.fabs(self.Q[self.grid.active_links]))
            #Trying a different tolerance criteria, previous may not work for big grids
            #tol = np.mean(np.fabs(dQ/self.Q[self.grid.active_links]))
            #max_change = np.max(np.fabs(dQ/self.Q[self.grid.active_links]))
            tol = np.max(dQ)
            end = time.time()
            print("Number of iterations =", niter, "tolerance =", tol, " iteration time=",end-start) #,"max_change=",max_change, " max_dQ=",max_dQ
            niter += 1

    def steady_state_transport(self):
        if self.solutes==None:
            #Not really a field error, need to find more appropriate error to raise
            raise FieldError(
                'Cannot run transport if no solutes are provided!')
        #Get index of nodes sorted by descending head
        sorted_idx = np.argsort(-self.h)#Sort on negative to get descending order sort
        #Loop through nodes from highest to lowest hydraulic head
        for node_idx in sorted_idx:
            #Determine link indicies and directions for this node
            this_node_links = self.grid.links_at_node[node_idx]
            this_node_dirs = self.grid.active_link_dirs_at_node[node_idx]

            #Determine which links have inflow and outflow to/from node
            node_Qs = self.grid.at_link['conduit__discharge'][this_node_links]*this_node_dirs
            node_inflow_links = this_node_links[node_Qs>0]
            node_outflow_links = this_node_links[node_Qs<0]

            #Calculate node concentration from inflows (including boundary inflows)
            total_inflow = np.abs(self.grid.at_link['conduit__discharge'][node_inflow_links]).sum() + \
                            self.grid.at_node['input__discharge'][node_idx]

            for solute in self.solutes:
                #avoid cases with no inflow (e.g. boundaries)
                if total_inflow !=0:
                    conc_mult = np.abs(self.grid.at_link['conduit__discharge'][node_inflow_links])*self.grid.at_link[solute+'__conc_out'][node_inflow_links]
                    input_conc_mult = self.grid.at_node['input__discharge'][node_idx]*self.grid.at_node[solute+'__inflow_conc'][node_idx]
                    conc_node = (conc_mult.sum() + input_conc_mult)/total_inflow

                    self.grid.at_node['concentration__'+solute][node_idx] = conc_node

                    #set outflow link concentrations
                    self.grid.at_link[solute+'__conc_in'][node_outflow_links] = conc_node

                else: #No inflow links (e.g. boundary node)
                    #Assume node concentration is fixed and set outflow conc to node conc
                    self.grid.at_link[solute+'__conc_in'][node_outflow_links] = self.grid.at_node['concentration__'+solute][node_idx]

                #Calculate output link concentrations (for now conservative)
                self.grid.at_link[solute+'__conc_out'][node_outflow_links] = self.grid.at_link[solute+'__conc_in'][node_outflow_links]
