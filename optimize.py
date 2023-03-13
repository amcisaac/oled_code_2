import os
from oled_code.mfss_helper import *
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
        rates = get_rates(rate_dict,jobtyp)
        print(Nfeval,rate_dict,objective(rates,jobtyp,fit,x_min,x_max,disr,Vol,shift))
    return 

inputfile = sys.argv[1]
jobname,exptfile,jobtyp,disr,Vol,shift,rate_dict,plot_title,plot_label = parse_input(inputfile)
rates=get_rates(rate_dict,jobtyp)
expt_x, expt_y, x_min, x_max, fit = read_file(exptfile)

print(objective(rates,jobtyp,fit,x_min,x_max,disr,Vol,shift))


res = spo.minimize(objective,rates, args = (jobtyp,fit,x_min,x_max,disr,Vol,shift), method='Nelder-Mead',options={'disp': True}, callback=callbackF)
#res = spo.minimize(objective,rates, args = (typ,fit,xmin,xmax,disr,Vol,shift), method='Powell',options={'disp': True}, callback=callbackF)
#res = spo.minimize(objective,rates, args = (typ,fit,xmin,xmax,disr,Vol,shift), method='CG',jac=False,options={'disp': True}, callback=callbackF)
print(res.x)

filepathlist = exptfile.split('/')
eqefile = filepathlist.pop()
plqyfilename = eqefile.split('_')[0] + "_plqy.csv"
plqyfile = '/'.join(filepathlist)+'/'+plqyfilename

if disr == 1:
    disr_str = 'eea'
elif disr == 3:
    disr_str = 'eca'
output_dict = {'figname': jobname+"_optimized",
               'datafile': [exptfile,plqyfile],
               'type': jobtyp,
               'disr': disr_str,
               'Vol': Vol1,
               'shift':shift,
               'rate_dict':
                   {
                   "rec": res.x[0],
                   "kr": res.x[1],
                   "eea": res.x[2],
                   "eca": res.x[3],
                   "sigma": res.x[4],
                   "kbi": res.x[5],
                   "PLQY": res.x[6],
                   "brec": res.x[7]
                   },
               'plot_title': plot_title,
               'plot_label': plot_label
               }


outputfile = inputfile.split('.')[0] + ".out"
with open(outputfile,'w') as jsonfile:
    json.dump(output_dict,jsonfile)

