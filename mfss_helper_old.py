import csv
import subprocess as sp
import numpy as np
import scipy.optimize as spo
from scipy.integrate import simps
import sys
import json
from inspect import getsourcefile
import os

def mfss(rates,typ,disr=1,length=1e-7,shift=False):
    '''
    python wrapper for the C++ MFSS codes. runs the SCF.

    Inputs:
        rates (list): rates with the rates
        typ (string): determines the type of calculation. options are 'eqe' or 'plqy'
        disr (int): which rate to put the disorder in. 1 = eea, 3 = eca
                    (the integer corresponds to the index of the rate in the C++)
        Vol (float): box side length in cm (only matters for PLQY). cube to get volume
        shift (bool): if True, will set disorder to 0 and shift the mean by sigma/2.
                      used to test effect of disorder vs shifted mean of distribution.

    Outputs:
        x_array (np array): current density (for eqe) or excitation density (for plqy)
        y_array (np array): EQE or PLQY corresponding to the values in x_array
        eea_array (np array): EEA rate as a function of x_array. keea * [exciton]**2
        eca_array (np array): ECA rate as a function of x_array. keea*[charge][exciton]
    '''
    Area = length**2
    Vol = length**3
    rec,kr,eea,eca,sigma,kbi,PLQY,brec = rates # SIGMA IS VARIANCE not std
    path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
    de = 100
    if sigma < 0:
        sigma = 0

    if disr == 1:
        eea = -np.log(eea)
        meanx = eea
        if shift:
            meanx -= sigma/2
            eea -= sigma/2
            sigma = 0
    elif disr == 3:
        eca = -np.log(eca)
        meanx = eca
        if shift:
            meanx -= sigma/2
            eca -= sigma/2
            sigma = 0


    output_array=[]
    if typ == 'eqe':
        spacing = np.linspace(-12,1,100)
        #spacing = np.linspace(-12,.01,40)
        for ki in spacing:
            output = sp.check_output([path+"/MF","--rec %f"%(rec), "--kr %f"%(kr),"--eea %f"%(eea),"--eca %f"%(eca),"--meanY %f"%(-ki),"--sigmaY %f"%(0), "--meanX %f"%(meanx),"--sigmaX %f"%(sigma),"--kbi %f"%(kbi), "--beea %f"%(0), "--beca %f"%(0), "--brec %f"%(brec),"--ki %f"%(1000000), "--dc %f"%(de), "--bdc %f"%(de), "--de %f"%(de), "--bde %f"%(de), "--disr %d"%(disr)], stderr=None)
            output_array.append(output.split(b" "))
    elif typ == 'plqy':
        spacing = np.linspace(-10,5,40)
        #spacing = np.linspace(-10,5,40)
        for ki in spacing:
            output = sp.check_output([path+"/MFPL","--rec %f"%(rec), "--kr %f"%(kr),"--eea %f"%(eea),"--eca %f"%(eca),"--meanY %f"%(-ki),"--sigmaY %f"%(0), "--meanX %f"%(meanx),"--sigmaX %f"%(sigma),"--kbi %f"%(kbi), "--beea %f"%(0), "--beca %f"%(0), "--brec %f"%(brec),"--ki %f"%(1000000), "--dc %f"%(de), "--bdc %f"%(de), "--de %f"%(de), "--bde %f"%(de), "--disr %d"%(disr)], stderr=None)
            output_array.append(output.split(b" "))

    data = np.array(output_array, np.float)


    V = data[:,0]
    C = data[:,1]
    E = data[:,2]

    V_array = np.array(V, np.float)
    C_array = np.array(C, np.float)
    E_array = np.array(E, np.float)
    ki_array = np.exp(spacing)

    if typ == 'eqe':
        x_array = V_array * ki_array * kbi   # current = ki * phi ; ki = kbi * exp(spacing) as defined in C++ code
        lum_array = kr * E_array * PLQY
        y_array = lum_array/x_array
        x_array *= 1.60217*np.power(10.0,-16)/Area # unit change (e/s -> mA/s), divide by area to get current density

    elif typ == 'plqy':
        emit = kr * E_array * PLQY
        absorbed = kr*E_array/PLQY + eea*E_array*E_array + eca*C_array*E_array
        rel_plqy = emit/ (np.exp(spacing)*V_array*kbi) #absorbed
        x_array = E_array/Vol
        y_array = rel_plqy/rel_plqy[0]

    eea_array = np.exp(-eea)*E_array*E_array
    eca_array = eca*E_array*C_array

    return x_array,y_array,eea_array,eca_array

def exp_fit(x_array,typ,fit,xmin,xmax):
    '''
    Returns y values of a polynomial fit for the raw y data, to make optimization smoother.

    Inputs:
        x_array: list or array of x data for which to get y values.
        typ: 'eqe' or 'plqy'. only matters for deciding which data is too big/small to fit.
        fit: list of strings, corresponding to the coefficients to use in the polynomial fit.
             will choose polynomial order based on length of fit.
        xmin: minimum value of the x_array to include in the fit.
        xmax: maximum value of the x_array to include in the fit.

    Returns:
        ef: list containing
    '''
    ef = []
    for x in x_array:
        if x < xmin or x > xmax:
            ef.append(-1)
        #elif typ == 'eqe' and (x < 1e-3 or x > 15): # throw out data too big/small to be reliable
        #    ef.append(-1)
        elif typ == 'plqy' and x < 1e17:
            ef.append(-1)
        else:
            x = np.log10(x)
            xfit = 0
            for i,param in enumerate(fit):
                xfit += float(param)*x**i
            ef.append(xfit)
    return ef



def objective(rates,typ,fit,xmin,xmax,disr,Vol,shift):
    #data = mfss(rates,typ,disr)
    x_array, y_array, eea_array, eca_array = mfss(rates,typ,disr,Vol,shift)

    expf = np.array(exp_fit(x_array,typ,fit,xmin,xmax)) #evaluates 'exp curve' at comparable currents
    resy = np.abs(expf - y_array)
    resy = resy[expf != -1]
    resx = x_array[expf != -1]
    if 'typ' == 'plqy':
        resx /= 1e17
    lsq = simps(resy,resx)
    return lsq

def objective_vector(rates,typ,fit,xmin,xmax,disr,Vol,shift):
    x_array, y_array, eea_array, eca_array = mfss(rates,typ,disr,Vol,shift)

    expf = np.array(exp_fit(x_array,typ,fit,xmin,xmax)) #evaluates 'exp curve' at comparable currents
    resy = expf - y_array       # do we want absolute value?!
    resy = resy[expf != -1] * 100
    test = np.divide(resy,y_array[expf!=-1]) # divide by experimental y data
    return test


def objective_both(rates,typ,fit,xmin,xmax,disr,Vol,shift):
    '''
    Objective function to minimize in optimize.py. Returns the difference between the (fit to)
    experimental data and the MFSS model for given parameters.

    Inptuts:
        rates: list of rates
        typ: 'eqe' or 'plqy'
        fit: list of strings corresponding to coefficients in polynomial fit to expt data
        xmin: minimum value of x data to include in fit
        xmax: maximum value of x data to include in fit
        disr: 1 if disorder in EEA, 3 if in ECA
        Vol: box side length (for PLQY)
        shift: True --> do shifted mean calc
    '''
    if typ == 'both':
        types = ['eqe','plqy']
    else:
        types = [typ]
    print(types)
    lsq = 0
    for i,jobtyp in enumerate(types):
        fiti = fit[i]
        x_mini = xmin[i]
        x_maxi = xmax[i]

        x_array, y_array, eea_array, eca_array = mfss(rates,jobtyp,disr,Vol,shift)
        print(x_array)
        expf = np.array(exp_fit(x_array,jobtyp,fiti,x_mini,x_maxi)) #evaluates 'exp curve' at comparable currents
        print(expf)
        resy = np.abs(expf - y_array)
        resy = resy[expf != -1]
        resx = x_array[expf != -1]
        #print(resy)
        #print(resx)
        lsq += simps(resy,resx)
        print(lsq)
    return lsq,x_array,expf


def parse_input(filename):
    '''
    converts json input file into parameters for the calculation

    Inputs:
        filename: can be string corresponding to name of input file, or dictionary
                  with multiple input files
    Outputs:
        jobname: name for any files produced
        exptfile: data file with experimental data
        jobtyp: what type of curve to calculate. options are 'eqe' or 'plqy'
        disr: which variable has disorder (1 = eea, 3 = eca)
        Vol: box side length for PLQY calculation. will be cubed to get volume
        shift: True --> do shifted mean calculation
        rates: list of rates for the calculation
        plot_title: title for plot
        plot_label: label for legend for plot
    '''
    if type(filename) == str:
        with open(filename,'r') as jsonfile:
            inputs = json.load(jsonfile)
    #elif type(filename) == dict:
    #    inputs = filename

    jobname = inputs['figname']
    exptfile = inputs['datafile']
    jobtyp = inputs['type']

    if inputs['disr'] == 'eea':
        disr = 1
    elif inputs['disr'] == 'eca':
        disr = 3

    Vol = inputs['Vol']
    shift = inputs['shift']=="True"

    rate_dict = inputs['rate_dict']
    #rates = [rate_dict['rec'],rate_dict['kr'],rate_dict['eea'],rate_dict['eca'],
    #         rate_dict['sigma'],rate_dict['kbi'],rate_dict['PLQY'],rate_dict['brec'],rate_dict['kbie']]
    plot_title = inputs['plot_title']
    plot_label = inputs['plot_label']
    return jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label

def get_rates(rate_dict,typ):
    '''
    Defines rate list using rate_dict and the job type. If eqe,
    defines kbi = kbi (for charge injection), if plqy, defines kbi = 0 (for exciton injection)

    Inputs:
        rate_dict: dictionary of rates, output from parse_inputs
        typ: job type ('eqe' or 'plqy') to define kbi
    Returns:
        rates: list of rates, in the order [rec,kr,eea,eca,sigma,kbi,PLQY,brec]
    '''
    if typ == 'eqe':
        rates = [rate_dict['rec'],rate_dict['kr'],rate_dict['eea'],rate_dict['eca'],
                 rate_dict['sigma'],rate_dict['kbi'],rate_dict['PLQY'],rate_dict['brec']]
    elif typ == 'plqy':
        try:
           kbie = rate_dict['kbie']
        except KeyError:
           kbie = 1
        rates = [rate_dict['rec'],rate_dict['kr'],rate_dict['eea'],rate_dict['eca'],
                 rate_dict['sigma'],kbie,rate_dict['PLQY'],rate_dict['brec']]
    return rates

def read_file(dat,imin=0,imax=-1):
    '''
    Reads the experimental data file and returns the data in lists.

    Inputs:
        dat: string with data file name
        imin (optional): index of first data point to include
        imax (optional): index of last data point to include

    Returns:
        exp_x: experimental x data
        exp_y: experimental y data
        x_min: lowest x data value
        x_max: highest x data value
        fit: coefficients for polynomial fit. order of polynomial
             determined by number of coefficients.
    '''
    exp_x = []
    exp_y = []
    with open(dat, 'r') as csvfile:
        rateReader = csv.reader(csvfile, delimiter=',')
        next(rateReader)
        fit = next(rateReader)
        fit.reverse()
        for rowS in rateReader:
            row = np.array([float(x) for x in rowS])
            exp_y.append(row[1])
            exp_x.append(row[0])

    x_min = exp_x[imin]
    x_max = exp_x[imax]
    #ef = exp_fit(exp_x,fit,x_min,x_max)

    return exp_x,exp_y,x_min,x_max,fit


def param_sweep(rates):
    for i,rate in enumerate(rates):
        # EQE SWEEP
        xarray,yarray,eeaarray,ecaarray = mfss(rates,jobtyp,disr,Vol,shift)
        rates[i] *= 10   # high rate = 10*rate
        xarray_high,yarray_high,eeaarray_high,ecaarray_high = mfss(rates,jobtyp,disr,Vol,shift)
        rates[i] /= 100  # low rate = rate/10 = high rate/100
        xarray_low,yarray_low,eeaarray_low,ecaarray_low = mfss(rates,jobtyp,disr,Vol,shift)
        rates[i] *= 10   # put back to the normal rate

        plt.figure()
        plt.title(ratelabel + ' EQE sweep')

        plt.semilogx(expt_x,expt_y,'.',label="Experiment")
        plt.semilogx(mfss_x_high,mfss_y_high,label=str(rate*10),color='red')
        plt.semilogx(mfss_x, mfss_y,label=str(rate),color='blue')
        plt.semilogx(mfss_x_low,mfss_y_low,label=str(rate/10),color='green')
        plt.savefig(ratelabel+'_eqe_sweep.pdf')


if __name__ == '__main__':
    inputfile = sys.argv[1]

    jobname,exptfile1,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
    Vol = Vol**3
    rates = get_rates(rate_dict,jobtyp)
    xarray,yarray,eeaarray,ecaarray = mfss(rates,jobtyp,disr,Vol,shift)
    print("MFSS " + jobtyp + " values: ", yarray)
    if jobtyp == 'eqe':
        xlab = 'current density'
    else:
        xlab = 'excitation density'
    print("MFSS " + xlab + ": ", xarray)
