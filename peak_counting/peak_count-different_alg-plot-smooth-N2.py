# It only works for ONE run at the time (no averaging over runs)

# It finds and plots the numebr of peaks in the DOS for both DFt and Pred
#
#
#
# It needis info on the original binning to determine the prequencies at which the peaks occurs.
# signal.find_peaks_cwt outputs the indices of the locations in the vector where peaks were found.
# The list is sorted.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html

# Function "peakdet":
# https://gist.github.com/endolith/250860
# A point is considered a maximum peak if it has the maximal value, and was preceded (to the left)
# by a value lower by DELTA
# ==> delta = difference of y value below which 2 maxima are considered part of the same peak

from jarvis.db.jsonutils import loadjson
from jarvis.core.atoms import Atoms
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_absolute_error
from scipy.stats import kde
from scipy.stats import gaussian_kde

#from matplotlib.colors import Normalize, LogNorm, to_hex
#from matplotlib.cm import plasma, ScalarMappable

from scipy import signal
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array
from scipy.interpolate import UnivariateSpline

rr=['run11']
max_freq=1000
min_freq=0
step=5

fact_smooth=0.4
spline_degree=4
dd=[0.05]

#dd=[0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]

def main():
  name="Peak_count-DFT_Pred-"
  title="ALIGNN Phonons-"
  for it in rr:
      nn=it.split("n")
      name=name+"R"+str(nn[1])
      title=title+"R"+str(nn[1])
  title=title+"-splines_"+str(fact_smooth)+"-deg_"+str(spline_degree)+"-"+str(dd[0])
  f_name0=name+"-splines_"+str(fact_smooth)+"-deg_"+str(spline_degree)+"-"+str(dd[0])+".png"
  f_name=name+"-splines_"+str(fact_smooth)+"-deg_"+str(spline_degree)+"-"+str(dd[0])+".plot"
  print(f_name)
  f_out1=open(f_name,"w")

  line0="Counter  Delta   Jid     Npeaks_DFT    Npeaks_Pred\n"
  f_out1.write(line0)
  for run in rr:
    f_in='../'+str(run)+'/temp/multi_out_predictions.json'
    print(f_in)
    data=loadjson(f_in)
    n_materials = len(data)
    print("Number of materials in dataset=",n_materials)

    llen=len(data[0]['predictions'])
    print("Make sure the spectral parameters are correct for this data!!!!")
    print(" ")
    print("Input file=",f_in)
    print("Number of channels=",llen)
    print("Spectral parameters=: max_freq=",max_freq,"min_freq=",min_freq,"step=",step)
    nchan=(max_freq-min_freq)/step + 1
    if (nchan != llen):
        print("Wrong values for max_freq, min_freq and step!!")
        exit()

    for delta in dd:
      n_peaks_DFT=[]
      n_peaks_Pred=[]
      MAE_Npeaks=0.0
      io=0
      for ii in data:
           jid=ii["id"]
           jid0=jid.split(".")
           jid1=jid0[0].split("-")
           name_jid="JVASP-"+jid1[2]
           #if (jid == "POSCAR-JVASP-8024.vasp"):
           #if (jid == "POSCAR-JVASP-12852.vasp"):
           DFT=ii["target"]
           # Spline smoothing
           x1=[]
           y1_DFT=[]
           io=0
           for y in DFT:
               x1.append(min_freq+io*step)
               y1_DFT.append(y)
               io = io+1
           spl = UnivariateSpline(x1, y1_DFT, k=spline_degree)
           spl.set_smoothing_factor(fact_smooth)
           # Peak count
           peakind_DFT, mintab_DFT = peakdet(spl(x1), delta)
           #peakind_DFT = signal.find_peaks_cwt(DFT, width)
           np_DFT = len(peakind_DFT)
           n_peaks_DFT.append(np_DFT)

           Pred=ii["predictions"]
           # Spline smoothing
           x1=[]
           y1_Pred=[]
           io=0
           for y in Pred:
               x1.append(min_freq+io*step)
               y1_Pred.append(y)
               io = io+1
           spl = UnivariateSpline(x1, y1_Pred, k=spline_degree)
           spl.set_smoothing_factor(fact_smooth)
           # Peak count
           peakind_Pred, mintab_Pred = peakdet(spl(x1), delta)
           np_Pred = len(peakind_Pred)
           n_peaks_Pred.append(np_Pred)
           line0=str(io)+"  "+str(delta)+"  "+str(name_jid)+"  "+str(np_DFT)+"  "+str(np_Pred)+"\n"
           f_out1.write(line0)
           io = io+1
           #if (io > 2): break
           #break

  #histogram the data
  n_bins = 10
  bins = [n_bins, n_bins] # number of bins
  range1 = [[1,n_bins+1],[1,n_bins+1]]  # range inside which to calculate the histogram
  hh, locx, locy = np.histogram2d(n_peaks_DFT, n_peaks_Pred, bins=bins, range=range1)
  print("locx=",locx)
  print("locy=",locy)
  print("Histogram")
  print(hh)

  zz = np.array(hh)
  zz1 = zz.flatten()
  print(zz1)
  x0 = []
  y0 = []
  z0 = []
  ik=0
  for ii in range(n_bins):
      for jj in range(n_bins):
          #print(ii,jj,locx[ii],locy[jj],zz1[ik])
          if (int(zz1[ik]) > 0):
             x0.append(locx[ii])
             y0.append(locy[jj])
             z0.append(zz1[ik])
          ik = ik+1
  x1=np.array(x0)
  y1=np.array(y0)
  z1=np.array(z0)

  # Plot definitions
  xmax=n_bins+1
  ymax=xmax
  xmin=0
  ymin=xmin
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.xticks(np.arange(xmin, xmax+0.01, 1))
  plt.yticks(np.arange(ymin, ymax+0.01, 1))
  plt.plot([xmin,xmax], [ymin,ymax], 'k-', alpha=0.75, zorder=0)
  plt.xlabel('Number of DFT Peaks')
  plt.ylabel('Number of PREDICTED Peaks')
  #plt.ylabel('Predicted n-PowerFactor ($\\mu$W/mK$^2$)')
  plt.title(title)

  #plt.scatter(x3, y3, c=z3, s=50, cmap='jet', edgecolor='', marker='o')
  #plt.scatter(n_peaks_DFT, n_peaks_Pred, c='r', s=100, marker='.')
  #plt.scatter(x3, y3, c=z3, cmap='jet', edgecolor='', marker='o', s=50)
  plt.scatter(x1, y1, c=z1, cmap='jet', edgecolor='', marker='o', s=50)
  #plt.colorbar(ticks=np.linspace(0,10,2), label='density')
  #plt.colorbar(label='Error on individual prediction')
  plt.colorbar(label='Counts')
  plt.clim(0,100)

  plt.savefig(f_name0)
  #plt.show()



def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)


main()



