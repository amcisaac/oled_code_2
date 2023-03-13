from mfss_classes import *
import numpy as np
import matplotlib.pyplot as plt
import sys

inputfile = sys.argv[1]

syst = oled(inputfile)
syst.read_file()




eqe_x_array,eqe_y_array,eqe_eea_array,eqe_eca_array=mfss(syst,'eqe')
#syst.exp_fit(eqe_x_array,'eqe')
#expf = np.array(syst.eqe_ef)
#print(np.array(eqe_x_array)[expf != -1])

plqy_x_array,plqy_y_array,plqy_eea_array,plqy_eca_array=mfss(syst,'plqy')

#eqe_y_array = 100*eqe_y_array
#syst.eqe_exp_y = 100*np.array(syst.eqe_exp_y)

lsq=objective(syst.eqe_rates,syst,'eqe')
print('eqe fit error ',lsq)
#print(syst.eqe_rates)

#x,y = eqe_x_array,eqe_y_array
##print('x,',job)
##print('current density',',','eea rate',',','eca rate')
#print('current density',',','y value')
#for i in range(0,len(x)):
#    print(x[i],',',y[i])

plt.figure()
plt.semilogx(syst.eqe_exp_x[5:42],syst.eqe_exp_y[5:42],'.',label="Experiment")
plt.semilogx(eqe_x_array[eqe_x_array < syst.eqe_exp_x[-1]],eqe_y_array[ eqe_x_array < syst.eqe_exp_x[-1]],label="MFSS")
plt.title('EQE from '+inputfile)
plt.legend()

plt.figure()
plt.semilogx(eqe_x_array,eqe_eea_array/eqe_eca_array)
#plt.show()

#plt.figure()
#plt.plot(syst.eqe_exp_x[syst.eqe_exp_x > 1e-3],syst.eqe_exp_y[syst.eqe_exp_x > 1e-3],'.',label="Experiment")
#plt.plot(eqe_x_array[eqe_x_array < syst.eqe_exp_x[-1]],eqe_y_array[eqe_x_array < syst.eqe_exp_x[-1]],label="MFSS")
#plt.title('EQE from '+inputfile)
#plt.legend()
##plt.show()

#plt.figure()
#plt.plot(syst.plqy_exp_x,syst.plqy_exp_y,'.',label="Experiment")
#plt.plot(plqy_x_array[plqy_x_array < syst.plqy_exp_x[-1]],plqy_y_array[plqy_x_array < syst.plqy_exp_x[-1]],label="MFSS")
#plt.title('PLQY from '+inputfile)
#plt.legend()
#plt.show()



plt.figure()
plt.semilogx(syst.plqy_exp_x,syst.plqy_exp_y,'.',label="Experiment")
plt.semilogx(plqy_x_array[plqy_x_array < syst.plqy_exp_x[-1]],plqy_y_array[plqy_x_array < syst.plqy_exp_x[-1]],label="MFSS")
plt.title('PLQY from '+inputfile)
plt.legend()
plt.show()
