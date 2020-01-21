from landlab import Component, FieldError
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.optimize import fsolve, brentq
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

    Todini and Rossmann (2013), Unified Framework for Deriving
    Simultaneous Equation Algorithms for Water Distribution
    Networks, Journal of Hydraulic Engineering, 139, 5, 511-526.



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
                    shape_factor=np.pi, g=9.81, f=0.1, solutes=None,
                    transfer_func=None,
                    transfer_func_link_variables=[],
                    transfer_kwd_args={},
                    **kwds):
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
        self.Theta = 0.5 #Relaxation factor for flow solution iterations
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
            self.h = grid.add_zeros('node', 'hydraulic__head')
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
        self.transfer_func = transfer_func
        self.transfer_func_link_variables = transfer_func_link_variables
        self.transfer_kwd_args = transfer_kwd_args

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

    def run_one_step(self, transport=False, transfer_kwd_args=None):#, **kwds):
        self.calculate_flow()
        if transport:
            #Allow updating of kwd args to transfer function during run_one_step()
            if transfer_kwd_args!= None:
                self.transfer_kwd_args = transfer_kwd_args
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

                if self.transfer_func==None:
                    #Calculate output link concentrations (for conservative transport)
                    self.grid.at_link[solute+'__conc_out'][node_outflow_links] = self.grid.at_link[solute+'__conc_in'][node_outflow_links]
            if self.transfer_func != None:
                #Calculate output concentrations according to provided transfer function
                #First extract needed link concentrations
                concs_in = {}
                for solute in self.solutes:
                    concs_in[solute] = self.grid.at_link[solute+'__conc_in'][node_outflow_links]
                link_vars = {}
                for link_variable in self.transfer_func_link_variables:
                    link_vars[link_variable] = self.grid.at_link[link_variable][node_outflow_links]
                transfer_dict = self.transfer_func(**concs_in, **link_vars, **self.transfer_kwd_args, length = self.grid.length_of_link[node_outflow_links])
                solutes_out = transfer_dict['solutes']
                for solute in solutes_out:
                    self.grid.at_link[solute+'__conc_out'][node_outflow_links] = solutes_out[solute]
                new_link_value_dict = transfer_dict['new_link_values']
                #Update link values (e.g. rates for calculating conduit evolution)
                for link_var in new_link_value_dict:
                    self.grid.at_link[link_var][node_outflow_links] = new_link_value_dict[link_var]


    def dyn_wave_solution(self, dt=500., outflow_bnd_type = 'head'):
        """
        Dynamic wave solver based on algorithm used in SWMM. Inertial terms are neglected.
        Here we use this solver only as a means of calculating steady flow.
        """
        #Determine flow depths at upstream and downstream ends of links
        FUDGE = 0.000001
        active_links = self.grid.active_links
        head_nodes = self.grid.node_at_link_head[active_links]
        tail_nodes = self.grid.node_at_link_tail[active_links]
        #Store values at current time
        Q_old = self.Q.copy()
        h_old = self.h.copy()
        Q_sum_old = self.grid.calc_net_flux_at_node(self.Q)[self.grid.core_nodes]/self.grid.dx + self.grid.at_node['input__discharge'][self.grid.core_nodes]
        converged = False
        num_iterations = 0
        max_iterations = 50
        while not converged and num_iterations<max_iterations:
            h_head = self.h[head_nodes]
            h_tail = self.h[tail_nodes]
            #Calculate flow depths using offset and junction elevations
            y_head = h_head - self.grid.at_node['junction__elevation'][head_nodes] \
                            - self.grid.at_link['conduit_head__offset'][active_links]
            y_tail = h_tail - self.grid.at_node['junction__elevation'][tail_nodes] \
                            - self.grid.at_link['conduit_tail__offset'][active_links]
            #Check for cases where flow depth is above conduit ceiling and reset to max flow depth.
            y_head[y_head>self.grid.at_link['maximum__depth'][active_links]] = self.grid.at_link['maximum__depth'][active_links][y_head>self.grid.at_link['maximum__depth'][active_links]]
            y_tail[y_tail>self.grid.at_link['maximum__depth'][active_links]] = self.grid.at_link['maximum__depth'][active_links][y_tail>self.grid.at_link['maximum__depth'][active_links]]
            #Check for negative cases
            y_tail[y_tail<FUDGE] = FUDGE
            y_head[y_head<FUDGE] = FUDGE
            #head_is_higher_than_tail = h_head>h_tail
            #y_avg = np.zeros(len(y_tail))
            #y_avg[head_is_higher_than_tail] = y_head[head_is_higher_than_tail]
            #y_avg[~head_is_higher_than_tail] = y_tail[~head_is_higher_than_tail]
            #y_avg[y_avg<FUDGE] = FUDGE
            y_avg = 0.5*(y_head + y_tail) #This is used in SWMM, I think upstream y may work better for steady flow.
            #Using the upstream value seems to rid of strange downstream boundary effects that impact entire network
            #Calculate flow XC area for square XCs
            A_avg = self.grid.at_link['width'][active_links] * y_avg
            A_avg[A_avg<FUDGE] = FUDGE
            #print('y=',y_avg)
            #print('A=',A_avg)
            #Calculate hydraulic diameters and write into grid values
            self.d_h[active_links] = self.d_h_square(self.grid.at_link['width'][active_links], y_avg)
            #print('D_H=',self.d_h)
            U = self.Q[active_links]/A_avg
            #print('U=',U)
            dQ_pres = self.g*A_avg*(h_head - h_tail)/self.grid.length_of_link[active_links]*dt#Original had a negative sign in front of this term. seems to work without. (check this)
            dQ_fric = self.f*np.abs(U)/(2.*self.d_h[active_links])*dt
            #print('dQ_pres=',dQ_pres)
            #print('dQ_fric=',dQ_fric)
            Q_new = (Q_old[active_links] + dQ_pres) / (1. + dQ_fric)
            self.Q[active_links] = (1. - self.Theta)*self.Q[active_links] + self.Theta*Q_new
            dQ = abs(Q_new - Q_old[active_links])
            #max_percent_change = np.max((Q_old[active_links] - Q_new))
            #print ('Iteration: ', num_iterations, '  (Q_old - Q_new)= ',max_percent_change )
            #if max_percent_change<0.0001:
            #    converged = True
            #num_iterations += 1
            #Iteratively solve momentum equation using values of area, velocity, and d_h from
            #previous head and discharge values.
            #print('max h_old=',max(h_old))
            #print('max Q_old =', max(Q_old))
        #   Head iteration

            #Note: Changed the sign in discharge input in Q_sum (both old and new).
            #This seemed to fix test case 2. Maybe there is something funny with discharge
            #signs to begin with. Go back through landlab grid specs to check my procedure.
            Q_sum_new = self.grid.calc_net_flux_at_node(self.Q)[self.grid.core_nodes]/self.grid.dx + self.grid.at_node['input__discharge'][self.grid.core_nodes]
            dh = 0.5*dt*(Q_sum_new + Q_sum_old)/self.grid.at_node['storage'][self.grid.core_nodes]
            h_new = h_old[self.grid.core_nodes] + dh
            h_new = (1. - self.Theta)*self.h[self.grid.core_nodes] + self.Theta*h_new

            if outflow_bnd_type!='head':
                #Adjust boundary heads on open boundaries
                upwind_links = self.grid.upwind_links_at_node(self.Q)[self.grid.open_boundary_nodes]
                for i, row in enumerate(upwind_links):
                    bnd_node = self.grid.open_boundary_nodes[i]
                    boundary_link = row[row>0]
                    if np.size(boundary_link==1):
                        link_nodes = self.grid.nodes_at_link[boundary_link]
                        upstream_node = link_nodes[link_nodes != bnd_node]
                        equiv_upstream_flow_depth = self.grid.at_node['hydraulic__head'][upstream_node] - self.grid.at_node['junction__elevation'][upstream_node]
                        width = self.grid.at_link['width'][boundary_link]
                        bnd_Q = self.Q[boundary_link]
                        #If full pipe set head to ceiling of conduit
                        if outflow_bnd_type=='normal':
                            if equiv_upstream_flow_depth>=self.grid.at_link['maximum__depth'][boundary_link]:
                                self.grid.at_node['hydraulic__head'][bnd_node] = self.grid.at_link['maximum__depth'][boundary_link] + self.grid.at_node['junction__elevation'][bnd_node]#equiv_upstream_flow_depth*0.95
                            else:
                                print("Entering normal flow calc.")
                                #Set depth to that of normal flow given the current discharge
                                slope = (self.grid.at_node['junction__elevation'][upstream_node] - self.grid.at_node['junction__elevation'][bnd_node])/self.grid.length_of_link[boundary_link]
                                print('slope=',slope, '  Q=',self.Q[boundary_link])
                                if slope >0 and abs(self.Q[boundary_link])>FUDGE: #This fails for flat conduits or conduits with zero Q
                                    print("equiv_upstream_flow_depth=",equiv_upstream_flow_depth)
                                    max_y = self.grid.at_link['maximum__depth'][boundary_link]
    #                                print(slope, width, bnd_Q)
                                    print('f(a)=',self.normal_flow_residual(equiv_upstream_flow_depth/2., slope, bnd_Q,width), '  f(b)=',self.normal_flow_residual(max_y, slope, bnd_Q,width))
                                    if self.normal_flow_residual(FUDGE, slope, bnd_Q,width)*self.normal_flow_residual(max_y, slope, bnd_Q,width)<0:
                                        y_norm = brentq(self.normal_flow_residual, FUDGE, max_y, args=(slope,bnd_Q,width))
                                        print('y_norm=',y_norm)
                                        self.grid.at_node['hydraulic__head'][bnd_node] = y_norm + self.grid.at_node['junction__elevation'][bnd_node]
                        elif outflow_bnd_type=='outfall':
                            if equiv_upstream_flow_depth>self.grid.at_link['maximum__depth'][boundary_link]:
                                equiv_upstream_flow_depth=self.grid.at_link['maximum__depth'][boundary_link]
                            #This works for square conduit only
                            y_crit = (bnd_Q**2/width**2/self.g)**(1./3.)
                            self.grid.at_node['hydraulic__head'][bnd_node] = y_crit + self.grid.at_node['junction__elevation'][bnd_node]

            #Check for convergence
            if num_iterations>0:
                max_change_h = max(h_new - self.h[self.grid.core_nodes])
                print('max change in h: ', max_change_h)
                if max_change_h < 0.0001:
                    converged=True
            self.h[self.grid.core_nodes] = h_new
            num_iterations+=1
        print("average dh=",np.mean(dh), '  average abs(dQ)=', np.mean(dQ))



        #return self.Q

    def d_h_square(self, width, flow_depth):
        d_H = np.zeros(np.size(width))
        is_full_pipe = np.isclose(width,flow_depth)
        if np.size(width) > 1:
            #print("width =",width, " flow_depth=",flow_depth)
            d_H[is_full_pipe] = width[is_full_pipe]
            d_H[~is_full_pipe] = 4.*width[~is_full_pipe]*flow_depth[~is_full_pipe] / (2.*flow_depth[~is_full_pipe] + width[~is_full_pipe])
        else:
            if is_full_pipe:
                d_H = width
            else:
                d_H = 4.*width*flow_depth/(2.*flow_depth + width)
        return d_H

    def normal_flow_residual(self, y, slope, Q, width):
        return slope - self.f *Q**2./(2.*self.g * self.d_h_square(width,y) * width**2 * y**2)
