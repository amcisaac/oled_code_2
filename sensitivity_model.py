import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
from mfss_helper import *
import sys
import matplotlib.pyplot as plt

problem = { 'num_vars': 8,
            'names': ['rec','kr','eea','eca','sigma','kbi','PLQY','brec'],
            'bounds': [[1e+06, 1e+10],
                       [1e+03, 1e+07],
                       [1e+05, 1e+11],
                       [1e+05, 1e+11],
                       [0.0, 10.0],
                       [1e+02, 1e+06],
                       [0.0,1.0],
                       [100, 1e+06]] }

param_values = saltelli.sample(problem,10)
Y = np.full([param_values.shape[0]],np.zeros(100))

inputfile = sys.argv[1]

jobname,exptfile,jobtyp,disr,Vol,shift,rates,plot_title,plot_label = parse_input(inputfile)
Vol = Vol**3
expt_x, expt_y, x_min, x_max, fit = read_file(exptfile)

for i,X in enumerate(param_values):
    x,y,eea,eca = mfss(X, jobtyp,disr,Vol,shift)
    fitparams = np.polyfit(x,y,5)
    stdx = np.linspace(1e-2,10,100)
    Y[i] = fitparams[0]*stdx**5+fitparams[1]*stdx**4+fitparams[2]*stdx**3+fitparams[3]*stdx**2+fitparams[4]*stdx+fitparams[5]
    plt.figure()
    plt.plot(x,y,label='mfss')
    plt.plot(stdx,Y[i],label='fit')
    plt.legend()
    plt.show()

#Si = sobol.analyze(problem,Y,print_to_console=True) 
