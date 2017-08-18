from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
from landlab.components import PresFlowNetwork

mg = RasterModelGrid((100,150),100)
junc_elev = mg.add_zeros('node', 'junction__elevation')
d_h = mg.add_ones('link', 'hydraulic__diameter')
h = mg.add_zeros('node', 'hydraulic__head')
#h = np.random.rand(mg.number_of_nodes)
Q = mg.add_zeros('link', 'conduit__discharge')
#net_node_flux = mg.add_ones('node', 'net_node_flux', noclobber=False)
#set boundary node head

h[mg.nodes_at_left_edge] = 10.
h[mg.nodes_at_right_edge] = 0.
h[mg.nodes_at_top_edge] = 0.
h[mg.nodes_at_bottom_edge] = 0.
#mg.set_closed_boundaries_at_grid_edges(False,True,False,True)
Q[mg.active_links] = 1.
n_core = mg.number_of_core_nodes
links = mg.links_at_node
print "Number of links = ", mg.number_of_links
print "Number of nodes = ", mg.number_of_nodes
print "Number of active links = ", mg.number_of_active_links
print "Number of core nodes = ", mg.number_of_core_nodes
pfn = PresFlowNetwork(mg)
