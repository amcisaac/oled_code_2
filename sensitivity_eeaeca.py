import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
from mfss_helper import *
import sys


inputfile = sys.argv[1]
N = int(sys.argv[2]) # number of param sets to generate

jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
if jobtyp == 'both':
    jobtyp = 'eqe'
    exptfile += '_eqe.csv'
rates = get_rates(rate_dict,jobtyp)

try:
    lb = float(sys.argv[3])
    ub = float(sys.argv[4])
except IndexError:
    lb = 0.9
    ub = 1.1

# set up bounds for rates depending on optimized value (is this bad?)
bounds = [[rates[2]*lb,rates[2]*ub],[rates[3]*lb,rates[3]*ub]]

Vol = Vol**3
expt_x, expt_y, x_min, x_max, fit = read_file(exptfile)

# optional: can specify x_min and x_max to only test
# response of roll-off part of curve
#try:
#    x_min = float(sys.argv[3])
#    x_max = float(sys.argv[4])
#    print("Considering sensitivity of roll-off only")
#    print("x min: ",x_min)
#    print("x max: ",x_max)

#except IndexError:
#    pass

print("Sobol analysis from input file " + inputfile)
print("Coefficients determined with N={}".format(N))
print("Considering effect of keea, keca only. Bounds: {}k to {}k".format(lb,ub))

# set up problem to analyze
problem = { 'num_vars': 2,
            'names': ['eea','eca'],
            'bounds': bounds}

# generate parameter values using saltelli sampling
param_values = saltelli.sample(problem,N)
#print(param_values)

def obj_2(param_value,rates,jobtyp,fit,x_min,x_max,disr,Vol,shift):

    all_params = np.zeros(len(rates))
    all_params[0] = rates[0]
    all_params[1] = rates[1]
    all_params[2] = param_value[0]
    all_params[3] = param_value[1]
    all_params[4] = rates[4]
    all_params[5] = rates[5]
    all_params[6] = rates[6]
    all_params[7] = rates[7]
    lsq = objective(all_params, jobtyp, fit, x_min,x_max,disr,Vol,shift)
    return lsq

#print(rates)
#for x in param_values: 
#    print(x)
#    print(obj_2(x,rates,jobtyp,fit,x_min,x_max,disr,Vol,shift))
Y = np.zeros([param_values.shape[0]])

#print(objective(rates, jobtyp,fit,x_min,x_max,disr,Vol,shift))

# generate outputs (lsq) based on random params
for i,X in enumerate(param_values):
    Y[i] = obj_2(X,rates, jobtyp,fit,x_min,x_max,disr,Vol,shift)

#    if N <= 10:
#        print('sample ', i)
#        print(X)
#        print(Y[i])

# analyze outputs
# returns first, second order, and total sensitivity params
Si = sobol.analyze(problem,Y,print_to_console=True)

var_tot = np.var(Y)
mean_tot = np.mean(Y)
print('Mean lsq value in the system: ',mean_tot)
print('Total variance in the system: ',var_tot)
