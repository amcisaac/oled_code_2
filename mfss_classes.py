import time
import csv
import subprocess as sp
import numpy as np
import scipy.optimize as spo
from scipy.integrate import simps
import sys
import json
from inspect import getsourcefile
import os
import matplotlib.pyplot as plt


class oled:

    def __init__(self,filename):
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
        # open input
        with open(filename,'r') as jsonfile:
            inputs = json.load(jsonfile)

        # parse experimental data
        self.exptfile = inputs['datafile']
        self.eqe_file = self.exptfile.split('_')[0] + '_eqe.csv'
        self.plqy_file = self.exptfile.split('_')[0] + '_plqy.csv'

        # random parameters that don't really need to be used
        self.jobname = inputs['figname']
        self.jobtyp = inputs['type']
        self.plot_title = inputs['plot_title']
        self.plot_label = inputs['plot_label']
        self.shift = inputs['shift']=="True"

        # parse box size parameters
        self.length = inputs['Vol']
        self.Area = self.length**2
        self.Volume = self.length**3

        # parse rates
        self.rate_labels = ['rec','kr','eea','eca','sigma','kbi','PLQY','brec','kbie']
        self.rate_dict = inputs['rate_dict']
        self.eqe_rates = [self.rate_dict['rec'],self.rate_dict['kr'],self.rate_dict['eea'],self.rate_dict['eca'],
                         self.rate_dict['sigma'],self.rate_dict['kbi'],self.rate_dict['PLQY'],self.rate_dict['brec']]
        try:
            kbie = self.rate_dict['kbie']
        except KeyError:
            kbie = 1
        self.plqy_rates = [self.rate_dict['rec'],self.rate_dict['kr'],self.rate_dict['eea'],self.rate_dict['eca'],
                           self.rate_dict['sigma'],kbie,self.rate_dict['PLQY'],self.rate_dict['brec']]
        if inputs['disr'] == 'eea':
            self.disr = 1
        elif inputs['disr'] == 'eca':
            self.disr = 3
        elif inputs['disr'] == 'kr':
            self.disr = 0

        self.pl_spacing = []
        self.eqe_spacing = []

        print(self.eqe_rates[1]*(1-self.eqe_rates[6]))
        return

    def read_file(self,imin=0,imax=-1):
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
        eqe_exp_x = []
        eqe_exp_y = []
        with open(self.eqe_file, 'r') as csvfile:
            rateReader = csv.reader(csvfile, delimiter=',')
            next(rateReader)
            fit = next(rateReader)
            fit.reverse()
            self.eqe_fit = fit
            for rowS in rateReader:
                row = np.array([float(x) for x in rowS])
                eqe_exp_y.append(row[1])
                eqe_exp_x.append(row[0])

        self.eqe_exp_x = eqe_exp_x
        self.eqe_exp_y = eqe_exp_y
        self.eqe_x_min = eqe_exp_x[imin]
        self.eqe_x_max = eqe_exp_x[imax]


        plqy_exp_x = []
        plqy_exp_y = []
        with open(self.plqy_file, 'r') as csvfile:
            rateReader = csv.reader(csvfile, delimiter=',')
            fit = next(rateReader)
            fit.reverse()
            self.plqy_fit = fit
            for rowS in rateReader:
                row = np.array([float(x) for x in rowS])
                plqy_exp_y.append(row[1])
                plqy_exp_x.append(row[0])

        self.plqy_exp_x = np.array(plqy_exp_x)
        self.plqy_exp_y = np.array(plqy_exp_y)
        self.plqy_x_min = plqy_exp_x[imin]
        self.plqy_x_max = plqy_exp_x[imax]
        return

    def exp_fit(self,mfss_x,typ):
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
        if typ == 'eqe':
            eqe_ef = []
            for x in mfss_x:
                if x < self.eqe_x_min or x > self.eqe_x_max:
                    eqe_ef.append(-1)
                #elif (x < 1e-3 or x > 15): # throw out data too big/small to be reliable

                elif x > 40:
                    eqe_ef.append(-1)
                else:
                    x = np.log10(x)
                    xfit = 0
                    for i,param in enumerate(self.eqe_fit):
                        xfit += float(param)*x**i
                    eqe_ef.append(xfit)
            self.eqe_ef = eqe_ef

        elif typ == 'plqy':
            plqy_ef = []
            for x in mfss_x:
                if x < self.plqy_x_min or x > self.plqy_x_max:
                    plqy_ef.append(-1)
                #elif x < 1e17:
                #    plqy_ef.append(-1)
                else:
                    x = np.log10(x)
                    xfit = 0
                    for i,param in enumerate(self.plqy_fit):
                        xfit += float(param)*x**i
                    plqy_ef.append(xfit)
            self.plqy_ef = plqy_ef
        return


def mfss(oled,typ,verbose=False):
    '''
    python wrapper for the C++ MFSS codes. runs the SCF.

    Inputs:
        oled (instance of oled class): controls parameters such as rates, etc
        typ (string): determines the type of calculation. options are 'eqe' or 'plqy'

    Outputs:
        x_array (np array): current density (for eqe) or excitation density (for plqy)
        y_array (np array): EQE or PLQY corresponding to the values in x_array
        eea_array (np array): EEA rate as a function of x_array. keea * [exciton]**2
        eca_array (np array): ECA rate as a function of x_array. keea*[charge][exciton]
    '''
    # define which rates to use based on what type of calculation
    if typ == 'eqe':
        rates = oled.eqe_rates
    elif typ == 'plqy':
        rates = oled.plqy_rates

    rec,kr,eea,eca,sigma,kbi,PLQY,brec = rates # SIGMA IS VARIANCE not std
    de = 100

    # deal with disorder
    # disordered parameter has to be -ln(k) to work with C++/log normal stuff
    if sigma < 0:
        sigma = 0
    if oled.disr == 1:
        eea = -np.log(eea)
        meanx = eea
        if oled.shift:
            meanx -= sigma/2
            eea -= sigma/2
            sigma = 0
    elif oled.disr == 3:
        eca = -np.log(eca)
        meanx = eca
        if oled.shift:
            meanx -= sigma/2
            eca -= sigma/2
            sigma = 0
    elif oled.disr == 0:
        kr = -np.log(kr)
        meanx = kr
        if oled.shift:
            meanx -= np.sqrt(sigma)/2
            kr -= np.sqrt(sigma)/2
            sigma=0


    # run MFSS using C++ code
    path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
    output_array=[]
    if typ == 'eqe':
        if len(oled.eqe_spacing) > 0:
            spacing = oled.eqe_spacing
        else:
            spacing = np.linspace(-13,2,200)
        #print(oled.Area/((1e-7)**2))
        #eqe_corr = np.log(oled.Area/((1e-7)**2))
        #print(eqe_corr)
        #spacing = np.linspace(-12,.01,40)
        for ki in spacing:
            #ki -= eqe_corr
            output = sp.check_output([path+"/MF","--rec %f"%(rec), "--kr %f"%(kr),"--eea %f"%(eea),"--eca %f"%(eca),"--meanY %f"%(-ki),"--sigmaY %f"%(0), "--meanX %f"%(meanx),"--sigmaX %f"%(sigma),"--kbi %f"%(kbi), "--beea %f"%(0), "--beca %f"%(0), "--brec %f"%(brec),"--ki %f"%(1000000), "--dc %f"%(de), "--bdc %f"%(de), "--de %f"%(de), "--bde %f"%(de), "--disr %d"%(oled.disr)], stderr=None)
            output_array.append(output.split(b" "))
    elif typ == 'plqy':
        if len(oled.pl_spacing) > 0:
            spacing=oled.pl_spacing
        else:
            spacing = np.linspace(-5,7,200)
        #spacing = np.linspace(-10,5,40)
        for ki in spacing:
            output = sp.check_output([path+"/MFPL","--rec %f"%(rec), "--kr %f"%(kr),"--eea %f"%(eea),"--eca %f"%(eca),"--meanY %f"%(-ki),"--sigmaY %f"%(0), "--meanX %f"%(meanx),"--sigmaX %f"%(sigma),"--kbi %f"%(kbi), "--beea %f"%(0), "--beca %f"%(0), "--brec %f"%(brec),"--ki %f"%(1000000), "--dc %f"%(de), "--bdc %f"%(de), "--de %f"%(de), "--bde %f"%(de), "--disr %d"%(oled.disr)], stderr=None)
            output_array.append(output.split(b" "))

    data = np.array(output_array, np.float)


    # calculate experimental observables
    V = data[:,0]
    C = data[:,1]
    E = data[:,2]

    V_array = np.array(V, np.float)
    C_array = np.array(C, np.float)
    E_array = np.array(E, np.float)
    ki_array = np.exp(spacing)

    if verbose:
        print(np.all(V_array > 0),'vacancy')
        print(np.all(C_array > 0),'charge')
        print(np.all(E_array > 0),'exciton')

    if oled.disr == 0:
        kr = np.exp(-kr)
    if typ == 'eqe':
        x_array = V_array * ki_array * kbi   # current = ki * phi ; ki = kbi * exp(spacing) as defined in C++ code
        lum_array = kr * E_array * PLQY
        y_array = lum_array/x_array
        x_array *= 1.60217*np.power(10.0,-16)/oled.Area # unit change (e/s -> mA/s), divide by area to get current density

    elif typ == 'plqy':
        emit = kr * E_array * PLQY
        absorbed = kr*E_array/PLQY + eea*E_array*E_array + eca*C_array*E_array
        rel_plqy = emit/ (np.exp(spacing)*V_array*kbi) #absorbed
        x_array = E_array/oled.Volume
        y_array = rel_plqy/rel_plqy[0]

    # calculate actual rates
    if oled.disr ==1:
        eea_array = np.exp(-eea)*E_array*E_array
        eca_array = eca*E_array*C_array
    elif oled.disr == 3:
        eea_array = eea*E_array*E_array
        eca_array = np.exp(-eca)*E_array*C_array
    else:
        eea_array = eea*E_array*E_array
        eca_array = eca*E_array*C_array

    if verbose:
        print('k_i,curr_dens,V_array,C_array,E_array')
        for i in range(0,len(V_array)):
            print(kbi*ki_array[i],',',x_array[i],',',V_array[i],',',C_array[i],',',E_array[i])

    return x_array,y_array,eea_array,eca_array


def objective(eqe_rates,oled,typ):

    oled.eqe_rates = eqe_rates

    x_array, y_array, eea_array, eca_array = mfss(oled,typ)

    # calculate experimental curve at MFSS x values
    oled.exp_fit(x_array,typ)
    if typ == 'eqe':
        expf = np.array(oled.eqe_ef)
    elif typ == 'plqy':
        expf = np.array(oled.plqy_ef)

    resy = np.abs(expf - y_array)
    resy = resy[expf != -1]  # remove points where MFSS x values not in experimental curve
    resx = x_array[expf != -1]
    if typ == 'plqy':
        resx /= 1e17
    lsq = simps(resy,resx)
    #print(resx,resy)
    print(lsq,eqe_rates)
    # plt.figure()
    # plt.plot(resx,resy,'.-')
    # plt.show()
    return lsq

def objective_vector(oled,typ):
    x_array, y_array, eea_array, eca_array = mfss(oled,typ)

    # calculate experimental curve at MFSS x values
    oled.exp_fit(x_array,typ)
    if typ == 'eqe':
        expf = np.array(oled.eqe_ef)
    elif typ == 'plqy':
        expf = np.array(oled.plqy_ef)

    resy = expf - y_array       # do we want absolute value?!
    resy = resy[expf != -1] * 100  # where did *100 come from...
    test = np.divide(resy,expf[expf!=-1]) # divide by experimental y data
    return test


def objective_both(oled):
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
    types = ['eqe','plqy']
    lsq = 0
    for i,jobtyp in enumerate(types):

        x_array, y_array, eea_array, eca_array = mfss(oled,jobtyp)

        oled.exp_fit(x_array,jobtyp)
        if jobtyp == 'eqe':
            expf = np.array(oled.eqe_ef)
        elif jobtyp == 'plqy':
            expf = np.array(oled.plqy_ef)

        resy = np.abs(expf - y_array)
        resy = resy[expf != -1]
        resx = x_array[expf != -1]
        #print(resy)
        #print(resx)
        lsq += simps(resy,resx)
        print(lsq)
    return lsq,x_array,expf


def eqe_param_sweep(oled):
    jobtyp = 'eqe'
    rates = oled.eqe_rates
    for i,rate in enumerate(oled.eqe_rates):
        xarray,yarray,eeaarray,ecaarray = mfss(oled,jobtyp)
        oled.eqe_rates[i] *= 10   # high rate = 10*rate
        xarray_high,yarray_high,eeaarray_high,ecaarray_high = mfss(oled,jobtyp)
        oled.eqe_rates[i] /= 100  # low rate = rate/10 = high rate/100
        xarray_low,yarray_low,eeaarray_low,ecaarray_low = mfss(oled,jobtyp)
        oled.eqe_rates[i] *= 10   # put back to the normal rate

        plt.figure()
        plt.title(oled.rate_labels[i] + ' '+jobtyp+' sweep')

        plt.semilogx(oled.eqe_exp_x,oled.eqe_exp_y,'.',label="Experiment")
        plt.semilogx(xarray_high,yarray_high,label=str(rates[i]*10),color='red')
        plt.semilogx(xarray, yarray,label=str(rates[i]),color='blue')
        plt.semilogx(xarray_low,yarray_low,label=str(rates[i]/10),color='green')
        plt.legend()
        plt.savefig(oled.rate_labels[i]+'_'+jobtyp+'_sweep.pdf')
    return

def plqy_param_sweep(oled):
    jobtyp = 'plqy'
    rates = oled.plqy_rates
    for i,rate in enumerate(oled.plqy_rates):
        xarray,yarray,eeaarray,ecaarray = mfss(oled,jobtyp)
        oled.plqy_rates[i] *= 10   # high rate = 10*rate
        xarray_high,yarray_high,eeaarray_high,ecaarray_high = mfss(oled,jobtyp)
        oled.plqy_rates[i] /= 100  # low rate = rate/10 = high rate/100
        xarray_low,yarray_low,eeaarray_low,ecaarray_low = mfss(oled,jobtyp)
        oled.plqy_rates[i] *= 10   # put back to the normal rate

        plt.figure()
        plt.title(oled.rate_labels[i] + ' '+jobtyp+' sweep')

        plt.semilogx(oled.plqy_exp_x,oled.plqy_exp_y,'.',label="Experiment")
        plt.semilogx(xarray_high,yarray_high,label=str(rates[i]*10),color='red')
        plt.semilogx(xarray, yarray,label=str(rates[i]),color='blue')
        plt.semilogx(xarray_low,yarray_low,label=str(rates[i]/10),color='green')
        plt.legend()
        plt.savefig(oled.rate_labels[i]+'_'+jobtyp+'_sweep.pdf')


if __name__ == '__main__':
    starttime=time.time()
    inputfile = sys.argv[1]
    syst = oled(inputfile)
    syst.read_file()
    eqe_x_array,eqe_y_array,eqe_eea_array,eqe_eca_array=mfss(syst,'plqy',verbose=True)
    endtime = time.time()
    print(starttime-endtime)
    for i in range(0,len(eqe_x_array)):
        print(eqe_x_array[i],eqe_y_array[i])
    #jobname,exptfile1,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
    #Vol = Vol**3
    #rates = get_rates(rate_dict,jobtyp)
    #xarray,yarray,eeaarray,ecaarray = mfss(rates,jobtyp,disr,Vol,shift)
    #print("MFSS " + jobtyp + " values: ", yarray)
    #if jobtyp == 'eqe':
    #    xlab = 'current density'
    #else:
    #    xlab = 'excitation density'
    #print("MFSS " + xlab + ": ", xarray)
