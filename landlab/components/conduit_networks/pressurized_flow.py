from landlab import Component, FieldError
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

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
        'friction__factor'
    )

    _output_var_names = (
        'hydraulic__head',
        'conduit__discharge'
    )

    _var_units = {
        'junction__elevation':'m',
        'hydraulic__diameter':'m',
        'friction__factor':'unitless',
        'hydraulic__head':'m',
        'conduit__discharge':'m3/s'        
    }

    _var_mapping = {
        'junction__elevation': 'node',
        'hydraulic__diameter':'link',
        'friction__factor':'link',
        'hydraulic__head':'node',
        'conduit__discharge':'link',
    }

    _var_doc = {
        'junction__elevation':
            'elevation of the node junction relative to some datum',
        'hydraulic__diameter':
            'equivalent diameter of the conduit segment (D_H=4A_c/P_w)',
        'friction__factor':
            'friction factor used in flow equation',
        'hydraulic__head':
            'elevation head at node',
        'conduit__discharge':
            'flow discharge through conduit in direction of link',
    }
    
    def __init__(self, grid, flow_eqn='Darcy-Weisbach', **kwds):
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
        self.flow_eqn = flow_eqn
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
        if 'friction__factor' in grid.at_link:
            self.f = grid.at_link['friction__factor']
        else:
            raise FieldError(
                'Friction factors of conduits are required as a component input!')
        if 'hydraulic__head' in grid.at_node:
            self.h = grid.at_node['hydraulic__head']
        else:
            self.h = grid.add_zeros('node', 'hyraulic__head')
        if 'conduit__discharge' in grid.at_link:
            self.Q = grid.at_link['conduit__discharge']
        else:
            self.Q = grid.add_ones('link', 'conduit__discharge')

    def run_one_step(self, **kwds):
        #Calculate flow in network
        max_tol = 0.001
        tol = 1.
        niter = 1
        links = self.grid.links_at_node
        n_core = self.grid.number_of_core_nodes
        while tol>max_tol:
            a=2.
            r=1.
            Q_ij = self.Q[links][self.grid.core_nodes]*self.grid.active_link_dirs_at_node[self.grid.core_nodes]
            F = Q_ij.sum(axis=1)*(1.-1./a)
            A_ij = np.zeros([n_core, n_core])
            for i, this_node in enumerate(self.grid.core_nodes):
                #calculate diagonal term in A matrix
                node_Qs = Q_ij[i]
                A_ii = 0.
                for link_Q in node_Qs:
                    if link_Q != 0:
                        A_ii += 1./(a*r*np.fabs(link_Q)**(a-1.))
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
                            A_ij[i][j] = -1./(a*r*np.fabs(self.Q[this_link])**(a-1))
                        else:#if boundary add term to F
                            if self.Q[this_link]!=0:#Check whether this is a closed boundary
                                F[i] += 1./(a*r*np.fabs(self.Q[this_link])**(a-1))*self.h[neighbor_node]
            #Turn A into sparse matrix
            A_csr = csr_matrix(A_ij)
            ##Solve linear system for new approximation of heads
            self.h[self.grid.core_nodes] = spsolve(A_csr,F)
            dQ= -(1./a)*self.Q[self.grid.active_links] - 1./(a*r*np.fabs(self.Q[self.grid.active_links])**(a-1)) * (self.grid.calc_diff_at_link(self.h)[self.grid.active_links])
            self.Q[self.grid.active_links] +=  dQ
            tol = sum(np.fabs(dQ))/sum(np.fabs(self.Q[self.grid.active_links]))
            print "Number of iterations =", niter, "tolerance =", tol
            niter += 1
