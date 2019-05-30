from landlab import Component, FieldError
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.optimize import fsolve

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
                    shape_factor=np.pi, g=9.81, f=0.1, **kwds):
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


    def run_one_step(self, **kwds):
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
            Q_ij = self.Q[links][self.grid.core_nodes]*self.grid.active_link_dirs_at_node[self.grid.core_nodes]
            F = np.sum(Q_ij*(1. - 1./a[links][self.grid.core_nodes] ), axis=1)
            #Add recharge to nodes
            F += self.grid.at_node['input__discharge'][self.grid.core_nodes]
            A_ij = np.zeros([n_core, n_core])
            for i, this_node in enumerate(self.grid.core_nodes):
                #calculate diagonal term in A matrix
                node_Qs = Q_ij[i]
                node_as = a[links][self.grid.core_nodes][i]
                node_rs = r[links][self.grid.core_nodes][i]
                A_ii = 0.
                for link_idx, link_Q in enumerate(node_Qs):
                    if link_Q != 0:
                        A_ii += 1./(node_as[link_idx]*node_rs[link_idx]*np.fabs(link_Q)**(node_as[link_idx]-1.))
                A_ij[i][i] = A_ii
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
                            A_ij[i][j] = -1./(a[this_link]*r[this_link]*np.fabs(self.Q[this_link])**(a[this_link]-1))
                        else:#if boundary add term to F
                            if self.Q[this_link]!=0:#Check whether this is a closed boundary
                                F[i] += 1./(a[this_link]*r[this_link]*np.fabs(self.Q[this_link])**(a[this_link]-1))*self.h[neighbor_node]
            #Turn A into sparse matrix
            A_csr = csr_matrix(A_ij)
            ##Solve linear system for new approximation of heads
            self.h[self.grid.core_nodes] = spsolve(A_csr,F)
            dQ= -(1./a[self.grid.active_links])*self.Q[self.grid.active_links] - 1./(a[self.grid.active_links]*r[self.grid.active_links]*np.fabs(self.Q[self.grid.active_links])**(a[self.grid.active_links]-1)) * (self.grid.calc_diff_at_link(self.h)[self.grid.active_links])
            self.Q[self.grid.active_links] +=  dQ
            #tol = sum(np.fabs(dQ))/sum(np.fabs(self.Q[self.grid.active_links]))
            #Trying a different tolerance criteria, previous may not work for big grids
            #tol = np.mean(np.fabs(dQ/self.Q[self.grid.active_links]))
            #max_change = np.max(np.fabs(dQ/self.Q[self.grid.active_links]))
            tol = np.max(dQ)
            print "Number of iterations =", niter, "tolerance =", tol #,"max_change=",max_change, " max_dQ=",max_dQ
            niter += 1
