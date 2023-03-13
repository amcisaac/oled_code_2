import os
from oled_code.mfss_classes import *
import scipy.optimize as spo
import sys

Nfeval = 0
def callbackF(Xi):
    global Nfeval
    Nfeval += 1
    if Nfeval % 1 == 0:
        rate_dict=  {
            "rec": Xi[0],
            "kr": Xi[1],
            "eea": Xi[2],
            "eca": Xi[3],
            "sigma": Xi[4],
            "kbi": Xi[5],
            "PLQY": Xi[6],
            "brec": Xi[7]
            }
        rates = [rate_dict['rec'],rate_dict['kr'],rate_dict['eea'],rate_dict['eca'],
                 rate_dict['sigma'],rate_dict['kbi'],rate_dict['PLQY'],rate_dict['brec']]
        print(Nfeval,rate_dict,objective(rates,syst,'eqe'))
    return


inputfile = sys.argv[1]
syst = oled(inputfile)
syst.read_file()

print(objective(syst.eqe_rates,syst,'eqe'))
res = spo.minimize(objective,syst.eqe_rates, args = (syst,'eqe'), method='Nelder-Mead',options={'disp': True}, callback=callbackF)
#res = spo.minimize(objective,rates, args = (typ,fit,xmin,xmax,disr,Vol,shift), method='Powell',options={'disp': True}, callback=callbackF)
#res = spo.minimize(objective,rates, args = (typ,fit,xmin,xmax,disr,Vol,shift), method='CG',jac=False,options={'disp': True}, callback=callbackF)
print(res.x)
