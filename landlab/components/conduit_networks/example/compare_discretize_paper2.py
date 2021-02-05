from conduits1D import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition
import numpy as np

res1 = run1Dsim(nsteps=1000,
        nx=50, #number of nodes
        dt=500.,
        L=20e3,
        Ztype = 'sqrt', 
        Zmax = None, # linear change in ice thickness between max and min
        Zmin = None,
        Zslope = None, # or can set a slope
        Qpeak = 1.,
        Qbase = 6.,
        Qsteady=False,
        period = 24.*60.*60., 
        bedslope =0.,
        A_R = np.pi*2**2, # moulin cross-section area
        D0 = 2, # initial hydraulic diameter of the conduit
        hin=400., # initial upstream head
        hout=0.)

plot_3panels_singlestep(res1,500)

plot_2panel_overtime_multiposition(res1)

# path = 'figures_celia/3panels'
# plot_3panels_movie(res1,path,every=10)