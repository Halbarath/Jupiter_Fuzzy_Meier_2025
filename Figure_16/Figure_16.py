#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd

fs_label = 16
fs = 13
fs_legend = 13
fw = 'normal'

def main():
    B1 = np.loadtxt('541_timeseries.txt')
    B2 = np.loadtxt('542_timeseries.txt')
    B3 = np.loadtxt('544_timeseries.txt')
    xmin = np.min([np.min(B1[0,:]),np.min(B2[0,:]),np.min(B3[0,:])])
    xmax = np.max([np.max(B1[0,:]),np.max(B2[0,:]),np.max(B3[0,:])])
    ymin = 0.0
    
    colormap = plt.cm.tab10 
    
    fig, axs = plt.subplots(1,3,figsize=(14*0.85, 4*0.85))

    k, value = 0, 2
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:])])
    axs[k].plot(B1[0,:],B1[value,:],'r-',label=r'50:50, $10^8$')
    axs[k].plot(B2[0,:],B2[value,:],'g-',label=r'50:50, $10^9$')
    axs[k].plot(B3[0,:],B3[value,:],'b-',label=r'18:82, $10^9$')
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 1, 4
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:])])
    axs[k].plot(B1[0,:],B1[value,:],'r-',label=r'50:50, $10^8$')
    axs[k].plot(B2[0,:],B2[value,:],'g-',label=r'50:50, $10^9$')
    axs[k].plot(B3[0,:],B3[value,:],'b-',label=r'18:82, $10^9$')
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 2, 7
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:])])
    axs[k].plot(B1[0,:],B1[value,:],'r-',label=r'50:50, $10^8$')
    axs[k].plot(B2[0,:],B2[value,:],'g-',label=r'50:50, $10^9$')
    axs[k].plot(B3[0,:],B3[value,:],'b-',label=r'18:82, $10^9$')
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_yscale("log")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta\gamma} \ (M_{\oplus})$', fontsize=fs_label)
    axs[k].legend(loc='best',fontsize=fs_legend)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("jupiter_demixing_timeseries.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
