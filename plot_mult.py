import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from mfss_helper import *
import sys


Ninp = len(sys.argv)

plt.figure()
for N in range(1,Ninp):
    inputfile = sys.argv[N]

    jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
    Vol = Vol**3
    rates = get_rates(rate_dict,jobtyp)
    if N == 1:
        expt_x, expt_y, x_min, x_max, ef = read_file(exptfile)

    mfss_x, mfss_y, mfss_eea, mfss_eca = mfss(rates,jobtyp,disr,Vol,shift)

    if jobtyp == 'eqe':
        mfss_y = 100*np.array(mfss_y)
        expt_y = 100*np.array(expt_y)

    if N == 1:
        plt.semilogx(expt_x,expt_y,'.',label="Experiment")
    plt.semilogx(mfss_x,mfss_y,label=plot_label)



plt.title(plot_title)
if jobtyp == 'eqe':
    plt.xlabel("Current density (mA/cm^2)")
    plt.ylabel("External Quantum Efficiency (%)")

if jobtyp == 'plqy':
    plt.xlabel("Excitation density (cm^-3)")
    plt.ylabel("rel. PLQY")


plt.legend()
#plt.show()
plt.savefig(jobname+".pdf")
