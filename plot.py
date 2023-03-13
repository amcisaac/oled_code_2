import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from mfss_helper import *
import sys

def plot_all(jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label):
    '''
    plots the experiment vs mfss curve for eqe or plqy
    TO DO:
        allow this function to add mfss curves from multiple inputs
        could do by externally plotting expt curve, then having this function just do 
        the mfss part.
    '''
    rates = get_rates(rate_dict,jobtyp) # translate rate_dict --> list of rates; depends on job type
    print(jobtyp,rates[5])
    expt_x, expt_y, x_min, x_max, ef = read_file(exptfile)
    mfss_x, mfss_y, mfss_eea, mfss_eca = mfss(rates,jobtyp,disr,Vol,shift)
    if jobtyp == 'eqe':
        mfss_y = 100*np.array(mfss_y)
        expt_y = 100*np.array(expt_y)

    plt.figure()
    plt.semilogx(expt_x,expt_y,'.',label="Experiment")

    plt_x = mfss_x[mfss_x>x_min] # filter out small values of x
    plt_y = mfss_y[mfss_x>x_min]
    plt.semilogx(plt_x,plt_y,label=plot_label)

    plt.title(plot_title+ ' ' + jobtyp)
    if jobtyp == 'eqe':
        plt.xlabel("Current density (mA/cm^2)")
        plt.ylabel("External Quantum Efficiency (%)")
    
    if jobtyp == 'plqy':
        plt.xlabel("Excitation density (cm^-3)")
        plt.ylabel("rel. PLQY")
    

    plt.legend()
    plt.savefig(jobname+'_'+jobtyp+".pdf")
    plt.show()

inputfile = sys.argv[1]

jobname,exptfile1,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
Vol = Vol**3

if jobtyp == 'eqe':
    exptfile = exptfile1 + '_eqe.csv'
    plot_all(jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label)

elif jobtyp == 'plqy':
    exptfile = exptfile1 + '_plqy.csv'
    plot_all(jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label)

else: # both jobs
    exptfile = exptfile1 + '_eqe.csv'
    plot_all(jobname,exptfile,'eqe',disr,Vol,shift,rate_dict,plot_title,plot_label)

    exptfile = exptfile1 + '_plqy.csv'
    plot_all(jobname,exptfile,'plqy',disr,Vol,shift,rate_dict,plot_title,plot_label)
