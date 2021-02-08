from conduits1D import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition
import numpy as np
from collections import defaultdict
import pandas as pd

secinday = 24*3600
timestep = 300
supraglacial_baseflow = 0.1

def TimeStamps(time_start, time_end, timestep):     
    return np.arange(time_start, time_end+1, timestep)

def Qin_real(time, Qin_data, Qtime_data):
    return np.interp(time, Qtime_data, Qin_data)

def Fill_dict(Q_csv_name,head_csv_name,timestep):
    dictionnary = defaultdict(list)
    
    tmp1 = pd.read_csv(Q_csv_name)
    tmp1 = tmp1.dropna()
    Qin = tmp1.Qm3s.to_numpy() + supraglacial_baseflow
    Qtime = tmp1.SOY.to_numpy()
        
    #time array in seconds
    dictionnary['meltwater_time'] = TimeStamps(Qtime[0],Qtime[-1],timestep)
    dictionnary['meltwater_data'] = Qin_real(dictionnary['meltwater_time'], Qin, Qtime)
    
    tmp2 = pd.read_csv(head_csv_name)
    tmp2 = tmp2.dropna()
    dictionnary['h_real'] = tmp2.head_bed.to_numpy()
    dictionnary['t_real'] = tmp2.soy.to_numpy()
    return dictionnary

jeme = Fill_dict('surface_melt_jeme.csv','head_jeme.csv',timestep)
baseflow = 1 #m3/s
Qin_data = jeme['meltwater_data'][1000:15000]+baseflow
Qtime_data = jeme['meltwater_time'][1000:15000]-jeme['meltwater_time'][1000]
dt = 300
nsteps = len(Qin_data)
#%%
#Simulation for Jeme

res1 = run1Dsim(nsteps=nsteps,
        nx=50, #number of nodes
        dt=dt,
        L=20e3,
        Ztype = 'sqrt', 
        Zmax = None, # linear change in ice thickness between max and min
        Zmin = None,
        Zslope = None, # or can set a slope
        Qin_type = 'custom_array',
        Qpeak = 1,
        Qbase = 6,
        Qsteady = False,
        Qtime_data = Qtime_data, 
        Qin_data = Qin_data,
        period = 24.*60.*60., 
        bedslope =0.,
        A_R = np.pi*2**2, # moulin cross-section area
        D0 = 2, # initial hydraulic diameter of the conduit
        hin=400., # initial upstream head
        hout=0.)

plot_3panels_singlestep(res1,500)

from conduits1D import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition
plot_2panel_overtime_multiposition(res1)

#%%
path = 'figures_celia/3panels'
plot_3panels_movie(res1,path,every=10)

