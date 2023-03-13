import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
from mfss_helper import *
import sys

inputfile = sys.argv[1]
try:
    order = sys.argv[2]
except IndexError:
    order='second' # first order by default


print("{} order sensitivity analysis on {}".format(order,inputfile))


jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
if jobtyp == 'both':
    jobtyp = 'eqe'
    exptfile += "_eqe.csv"

rates = get_rates(rate_dict,jobtyp)
#Vol = Vol**3
expt_x, expt_y, x_min, x_max, fit = read_file(exptfile)
#print(expt_x)
names=['rec','kr','eea','eca','sigma','kbi','PLQY','brec']
#names = ['rec','kr','eea','eca','kbi','brec']
print("Considering " + jobtyp + " curve")

#if order == 'first':
#    chis = []
#    rates = np.array(rates)
#    for j,kj in enumerate(rates):
#        dj = kj * 0.000001
#        #print(dj)
#        dj_vec = np.zeros(len(rates))
#        dj_vec[j] = dj
#        chi_j = objective(rates+dj_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift) - objective(rates-dj_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift)
#        chi_j /= 2*dj
#        #chi_j /= kj
#        chis.append(chi_j)
#        print("S1 for {}: ".format(names[j]), chi_j**2)

r2 = 0.999710054455 # br
#r2 = 0.998575933508 # cbp
d = 0.01   # change for larger d
norm = False
print("Using d = " + str(d))
if norm: print("Normalizing derivatives")
else: print("Using unnormalized derivatives")
print('')
if order == 'second':
    chis = np.zeros((len(names)-3,len(names)-3))  # change dims for eliminating params
    rates = np.array(rates)
    for j,kj in enumerate(rates):
        if j != 4 and j != 6 and j != 7: # change bound on 6 for brec/not
            # dj = change in rate j
            dj = kj * d
            #print(dj)
            print(kj,dj)
            dj_vec = np.zeros(len(rates))
            if kj > -1e-6 and kj < 1e-6:
                print(kj)
                dj = d
                dj_vec[j] = dj
                chi_j = objective_vector(rates+dj_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift) - objective_vector(rates,jobtyp,fit,x_min,x_max,disr,Vol,shift)
                chi_j /= dj
            else:
                dj_vec[j] = dj

                # take derivative using middle point method-- chi_j is a vector with a value for each point
                chi_j1 = objective_vector(rates+dj_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift)
                chi_j2 = objective_vector(rates-dj_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift)

                chi_j =  chi_j1 - chi_j2
                #print(chi_j1)
                #print(chi_j2)
                #print(2*dj)
                chi_j /= 2*dj
            if norm: chi_j /= kj
            chi_j /= 1-r2 # change for no weighting


            for i,ki in enumerate(rates[j:len(rates)]):
                i += j # matrix will be symmetric--> only calculate upper part
                if i != 4 and i != 6 and i != 7: # change bound on 6 for brec/not
                    # di = change in rate i
                    di = ki * d
                    di_vec = np.zeros(len(rates))
                    if ki > -1e-6 and ki < 1e-6:
                        di = d
                        di_vec[i] = di
                        chi_i = objective_vector(rates+di_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift) - objective_vector(rates,jobtyp,fit,x_min,x_max,disr,Vol,shift)
                        chi_i /= di
                    else:
                        di_vec[i] = di

                        #calculate derivative with middle point method-- chi_i is a vector
                        chi_i1 = objective_vector(rates+di_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift)
                        chi_i2=objective_vector(rates-di_vec,jobtyp,fit,x_min,x_max,disr,Vol,shift)
                        #if i == 5 and j == 0:
                        #    chi_i1 = chi_i1[1:163]


                        chi_i = chi_i1 - chi_i2
                        chi_i /= 2*di
                        #print(chi_i)
                    if norm: chi_i /= ki
                    chi_i /= 1-r2 # change for no weighting
                    #print(j,i)
                    chi_ij = np.dot(chi_i,chi_j) # np.abs to do abs value
                    print("S2 for {} and {}: ".format(names[j],names[i]),chi_ij)
                    x = i
                    y = j
                    if i >4: x -= 1
                    if i > 6: x -= 1

                    if j >4: y -= 1
                    if j > 6: y -=1
                    chis[x][y] = chi_ij
                    chis[y][x] = chi_ij

    print('')
    print(chis.shape)
    #chis *= 10e14
    print(chis)
    c_inv = np.linalg.inv(chis)
    #c_inv *= 10e14
    for i in range(0,len(rates)):
        if i != 4 and i != 6 and i != 7: # change bound on 6 for brec/not
            x = i
            if i >4: x -= 1
            if i > 6: x -= 1
            print("Variance in {}: ".format(names[i]), c_inv[x][x])
    for i in range(0,len(rates)):
        if i != 4 and i != 6 and i != 7: # change bound on 6 for brec/not
            x = i
            if i >4: x -= 1
            if i > 6: x -= 1
            print("Error in {}: ".format(names[i]), np.sqrt(c_inv[x][x]))

    print('')
    print("Calculating eigenvalues and eigenvectors of covariance matrix")
    eigval,eigvec = np.linalg.eigh(chis)
    for i in range(0,len(eigval)):
        print(eigval[i], eigvec[:,i])
