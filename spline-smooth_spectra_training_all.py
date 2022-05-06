# It outputs a new json file, smoother that the initial one
# Algorithm:
# - when it finds a y=0 inbetween two channels with y >0, it substitute the zero with 
#   the average of the two nn channels
# - when it finds two y=0 inbetween two channels with y >0, it substitute the zero with
#   the weighted average of the two closest non-zero channels

from jarvis.db.jsonutils import loadjson
import numpy as np
import matplotlib.pyplot as plt
from jarvis.db.jsonutils import dumpjson
from scipy.interpolate import UnivariateSpline

name1="Freq-200_1800_5"
max_freq=1800
min_freq=-200
step=5

fact_smooth=0.2
spline_degree=3


f_out='PhononData-'+name1+'-splines_smooth.json'
f_in='PhononData-'+name1+'.json'
d=loadjson(f_in)
print(d[0].keys(), len(d))
#print("length spectral data=",len(d[0]['target']))

llen=len(d[0]['pdos_elast'])
print("Make sure the spectral parameters are correct for this data!!!!")
print(" ")
print("Input file=",f_in)
print("Number of channels=",llen)
print("Spectral parameters=: max_freq=",max_freq,"min_freq=",min_freq,"step=",step)
nchan=(max_freq-min_freq)/step + 1
if (nchan != llen):
    print("Wrong values for max_freq, min_freq and step!!")
    exit()

mem = []
#for index in range(1):
for index in range(len(d)):
   #print(index,'jid=',d[index]['jid'])
   jid=d[index]['jid']
   info = {}
   info["jid"] = jid
   info["atoms"] = d[index]['atoms']
   #print("")
   #print("pdos_elast")
   #print(d[index]['pdos_elast'])
   x1=[]
   y1=[]
   y2=[]
   io=0
   for y in d[index]['pdos_elast']: 
         x1.append(min_freq+io*step)
         y1.append(y)
         y2.append(y)
         io = io+1

   zero = 0.000001
   for ik in range(2,len(y1)-2):
         if (y1[ik] < zero):
             if ((y1[ik-1] > zero) and (y1[ik+1] > zero) ):
                 value = (y1[ik-1]+y1[ik+1])/2.0
             elif ((y1[ik-1] < zero) and (y1[ik+1] > zero) ):
                 value = (y1[ik-2]+2*y1[ik+1])/3.0
             elif ((y1[ik-1] > zero) and (y1[ik+1] < zero) ):
                 value = (2*y1[ik-1]+y1[ik+2])/3.0
             else:
                 value = y1[ik]
         else:
             value = y1[ik]
         y2[ik]=value
         #print("ik=",ik,"y1[ik]=",y1[ik],"y2[ik]=",y2[ik])
         
   spl = UnivariateSpline(x1, y2, k=spline_degree)
   spl.set_smoothing_factor(fact_smooth)
   # Alignn requires normalization
   mmax = np.max(sp1(x1))
   y3 = np.array(sp1(x1)/mmax)
   info['pdos_elast']=list(y3)
   mem.append(info)

dumpjson(data=mem, filename=f_out)
