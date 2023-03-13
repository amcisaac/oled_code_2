import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from mfss_helper import *
import sys

jobdict = {'eqe':['eqe'],'plqy':['plqy'],'both':['eqe','plqy']}


inputfile = sys.argv[1]

jobname,exptfiles,jobtyps,disr,Vol,shift,rates,plot_title,plot_label = parse_input(inputfile)
Vol = Vol**3
rates = get_rates(rates)

jobs = jobdict[jobtyps]

#print(jobs)
for i,jobtyp in enumerate(jobs):
    if type(exptfiles) == list: 
        exptfile = exptfiles[i]
    else:
        exptfile=exptfiles        
        #raise ValueError('Only one data file supplied, two jobs requested')

    expt_x, expt_y, x_min, x_max, ef = read_file(exptfile)
    mfss_x, mfss_y, mfss_eea, mfss_eca = mfss(rates,jobtyp,disr,Vol,shift)
    if jobtyp == 'eqe':
        mfss_y = 100*np.array(mfss_y)
        expt_y = 100*np.array(expt_y)

    plt.figure()
    plt.semilogx(expt_x,expt_y,'.',label="Experiment")
    plt.semilogx(mfss_x,mfss_y,label=plot_label)

    plt.title(plot_title+ ' ' + jobtyp)
    if jobtyp == 'eqe':
        plt.xlabel("Current density (mA/cm^2)")
        plt.ylabel("External Quantum Efficiency (%)")
    
    if jobtyp == 'plqy':
        plt.xlabel("Excitation density (cm^-3)")
        plt.ylabel("rel. PLQY")
    
    
    plt.legend()
    plt.savefig(jobname+'_'+jobtyp+".pdf")
