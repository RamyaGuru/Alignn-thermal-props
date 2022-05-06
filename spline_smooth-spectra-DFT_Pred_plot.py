# spline smoothing

from jarvis.db.jsonutils import loadjson
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

run='run17'
max_freq=1800
min_freq=-200
step=5


fact_smooth=0.2
spline_degree=2


f_in=str(run)+'/temp/multi_out_predictions.json'
d=loadjson(f_in)
print(d[0].keys(), len(d))
#print("length spectral data=",len(d[0]['target']))

llen=len(d[0]['predictions'])
print("Make sure the spectral parameters are correct for this data!!!!")
print(" ")
print("Input file=",f_in)
print("Number of channels=",llen)
print("Spectral parameters=: max_freq=",max_freq,"min_freq=",min_freq,"step=",step)
nchan=(max_freq-min_freq)/step + 1
if (nchan != llen):
    print("Wrong values for max_freq, min_freq and step!!")
    exit()

"""
print(" ")
print("Spectral parameters=: max_freq=",max_freq,"min_freq=",min_freq,"step=",step)
print("Make sure the spectral parameters are correct for this data!!!!")
print(" ")
"""

for index in range(60):
#for index in range(len(d)):
   #print(index,'jid=',d[index]['jid'])
   jid0=(d[index]['id'].split("-"))[2]
   jid_number=(jid0.split("."))[0]
   jid="JVASP-"+str(jid_number)
   #print("")
   #print("pdos_elast")
   #print(d[index]['pdos_elast'])
   f_plot="DFT-Pred_"+run+"_"+jid+"_"+"smooth_"+str(fact_smooth)+"-deg_"+str(spline_degree)+".png"
   print("file name=",f_plot)
   y_DFT=[]
   io=0
   for y in d[index]['target']:
         y_DFT.append(y)
         io = io+1

   x1=[]
   y1=[]
   io=0
   for y in d[index]['predictions']: 
         x1.append(min_freq+io*step)
         y1.append(y)
         io = io+1

   spl = UnivariateSpline(x1, y1, k=spline_degree)
   spl.set_smoothing_factor(fact_smooth)
   fig,ax=plt.subplots(figsize=(6,4))
   ax.plot(x1,y_DFT,c="b",lw=3,label='original DFT')
   ax.plot(x1,y1,"--k",label='original PREDICTED')
   #ax.plot(x1,spl(x1),c="r", label='smoothed DFT')
   ax.plot(x1,spl(x1),c="r", lw=3, label='spline smoothing')
   ax.set_xlabel("Frequency (cm-1)",fontsize=16)
   ax.xaxis.set_tick_params(labelsize=14,width=1.5)
   ax.yaxis.set_tick_params(labelsize=14,width=1.5)
   for axis in ['top','bottom','left','right']:
         ax.spines[axis].set_linewidth(1.5)
   ax.set_xlim(min_freq,max_freq)
   ax.set_ylabel("Intensity (arb. units)",fontsize=16)
   ax.legend()
   titolo='Smooth_factor='+str(fact_smooth)+'  spline_degree='+str(spline_degree)
   plt.title(titolo)
   #plt.title('Smooth_factor='+str(fact_smooth)+'  degree_spline=3')
   plt.tight_layout()
   #plt.plot(d[index]['target'])
   #plt.plot(d[index]['predictions'])
   plt.savefig(f_plot)
   plt.clf()
   #plt.savefig('pred_spectrum_1.png')
   #plt.show()
