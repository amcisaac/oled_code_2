import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
from mfss_helper import *
import sys


inputfile = sys.argv[1]
N = int(sys.argv[2]) # number of param sets to generate

jobname,exptfile,jobtyp,disr,Vol,shift,rates,plot_title,plot_label = parse_input(inputfile)

# set up bounds for rates depending on optimized value (is this bad?)
bounds = []
for i,rate in enumerate(rates):
#    if i == 4:
#        bounds.append([0,rate + 0.5*rate])
#    elif i == 6:
#        bounds.append([rate - rate*0.5, rate + rate*0.5])
#    else:
#        bounds.append([rate*10e-1,rate*10])
    bounds.append([rate/2,rate*2])

Vol = Vol**3
expt_x, expt_y, x_min, x_max, fit = read_file(exptfile)

# optional: can specify x_min and x_max to only test
# response of roll-off part of curve
try:
    x_min = float(sys.argv[3])
    x_max = float(sys.argv[4])
    print("Considering sensitivity of roll-off only")
    print("x min: ",x_min)
    print("x max: ",x_max)

except IndexError:
    pass

# set up problem to analyze
problem = { 'num_vars': 8,
            'names': ['rec','kr','eea','eca','sigma','kbi','PLQY','brec'],
            'bounds': bounds}

# generate parameter values using saltelli sampling
param_values = saltelli.sample(problem,N)
Y = np.zeros([param_values.shape[0]])

#print(objective(rates, jobtyp,fit,x_min,x_max,disr,Vol,shift))

# generate outputs (lsq) based on random params
for i,X in enumerate(param_values):
    Y[i] = objective(X, jobtyp,fit,x_min,x_max,disr,Vol,shift)
#    if N <= 10:
#        print('sample ', i)
#        print(X)
#        print(Y[i])

print("Coefficients determined with N={}") 
# analyze outputs
# returns first, second order, and total sensitivity params
Si = sobol.analyze(problem,Y,print_to_console=True) 
