from landlab import Component, FieldError
from landlab.grid.mappers import map_mean_of_link_nodes_to_link
import numpy as np

class MeltCreep(Component):
    """
    Calculate melt and creep within a glacial conduit network.

    Landlab component that solves for melt and creep in a 2D pressurized
    conduit network.

    Construction:

        MeltCreep(grid, [stuff to add later])

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid object.

    [add later]

    """
    _name = 'MeltCreep'

    _input_var_names = (
        'junction__elevation',
        'hydraulic__diameter',
        'hydraulic__head',
        'conduit__discharge',
        'ice__thickness',
    )

    _output_var_names = (
        'hydraulic__diameter'
    )

    _var_units = {
        'junction__elevation':'m',
        'hydraulic__diameter':'m',
        'hydraulic__head':'m',
        'conduit__discharge':'m3/s',
        'ice__thickness':'m',
    }

    _var_mapping = {
        'junction__elevation':'node',
        'hydraulic__diameter':'link',
        'hydraulic__head':'node',
        'conduit__discharge':'link',
        'ice__thickness':'node',
    }

    _var_doc = {
        'junction__elevation':
            'elevation of the node junction relative to some datum',
        'hydraulic__diameter':
            'equivalent diameter of the conduit segment (D_H=4A_c/P_w)',
        'hydraulic__head':
            'elevation head at node',
        'conduit__discharge':
            'flow discharge through conduit in direction of link',
        'ice__thickness':
            'thickness of ice at node',
    }

    def __init__(self, grid, shape_factor=np.pi, g=9.81, f=0.1,
                rho_w=998., rho_ice=900.,L_v=334000., n=3., B=58000000.,
                dt = 3600., **kwds):

        """
        Initialize the MeltCreep

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object

        [add more later]

        """
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        self.g = g
        self.f = f
        self.shape_factor = shape_factor
        self.rho_w = rho_w
        self.rho_ice = rho_ice
        self.L_v = L_v
        self.n = n
        self.B = B
        self.dt = dt
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
        if 'hydraulic__head' in grid.at_node:
            self.h = grid.at_node['hydraulic__head']
        else:
            raise FieldError(
                'Hydraulic head at nodes is a required component input!')
        if 'conduit__discharge' in grid.at_link:
            self.Q = grid.at_link['conduit__discharge']
        else:
            raise FieldError(
                'Conduit discharge is a required component input!'
            )
        if 'ice__thickness' in grid.at_node:
            self.link_ice_thickness = map_mean_of_link_nodes_to_link(self.grid, 'ice__thickness')
        else:
            raise FieldError(
                'Ice thickness at nodes is a required component input!'
            )

    def calc_melt(self):
        """
        Calculate ice melt within conduits
        """
        L = self._grid.length_of_link
        f = self.f
        g = self.g
        d_h = self.d_h
        r_s = self.shape_factor
        #Can make this more general later
        #Calculate melt
        self.melt = 16.*self.rho_w*f/(np.pi**3.*self.rho_ice*self.L_v)*abs(self._grid.at_link['conduit__discharge'])**3./d_h**6.

    def calc_creep(self):
        dP = self.g*(self.rho_ice*self.link_ice_thickness - self.rho_w*(self.link_head - self.link_elevation))
        self.creep = (1./(self.n*self.B))**3. * self.d_h*dP*np.fabs(dP**(self.n - 1.))

    def calc_a(self):
        """
        Calculate and set head loss exponent.
        """
        if self.flow_eqn == 'Darcy-Weisbach':
            self.a = 2.0*np.ones(self._grid.number_of_links)

    def run_one_step(self, **kwds):
        #Calculate mean heads in links
        self.link_head = map_mean_of_link_nodes_to_link(self._grid, 'hydraulic__head')
        self.link_elevation = map_mean_of_link_nodes_to_link(self._grid, 'junction__elevation')
        self.calc_melt()
        self.calc_creep()
#        print("melt = ", self.melt)
#        print("creep = ", self.creep)
        ddh = (self.melt - self.creep)*self.dt
        print(["mean ddh = ", ddh[self._grid.active_links].mean()])
        self.d_h[self._grid.active_links] += ddh[self._grid.active_links]
