#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd

fs_label = 16
fs = 13
fs_legend = 13
fw = 'normal'

def main():
    fig, axs = plt.subplots(1,2,figsize=(14.0*0.85, 5*0.85))
    
    colormap = plt.cm.tab10 #nipy_spectral, Set1,Paired    
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]

    value = 2

    k = 0
    B1 = np.loadtxt('771_timeseries.txt')
    B2 = np.loadtxt('772_timeseries.txt')
    B3 = np.loadtxt('773_timeseries.txt')
    B4 = np.loadtxt('774_timeseries.txt')
    B5 = np.loadtxt('775_timeseries.txt')
    B6 = np.loadtxt('776_timeseries.txt')
    B7 = np.loadtxt('777_timeseries.txt')
    B8 = np.loadtxt('778_timeseries.txt')
    B9 = np.loadtxt('779_timeseries.txt')
    B10 = np.loadtxt('7710_timeseries.txt')
    B1[value,:] /= (2.0/3.0*0.5)
    B2[value,:] /= (2.0/3.0*0.5)
    B3[value,:] /= (2.0/3.0*0.5)
    B4[value,:] /= (2.0/3.0*0.5)
    B5[value,:] /= (2.0/3.0*0.5)
    B6[value,:] /= (20.0/3.0*0.5)
    B7[value,:] /= (20.0/3.0*0.5)
    B8[value,:] /= (20.0/3.0*0.5)
    B9[value,:] /= (20.0/3.0*0.5)
    B10[value,:] /= (20.0/3.0*0.5)
    xmin = np.min([np.min(B1[0,:]),np.min(B2[0,:]),np.min(B3[0,:]),np.min(B4[0,:]),np.min(B5[0,:]),np.min(B6[0,:]),np.min(B7[0,:]),np.min(B8[0,:]),np.min(B9[0,:]),np.min(B10[0,:])])
    xmax = np.max([np.max(B1[0,:]),np.max(B2[0,:]),np.max(B3[0,:]),np.max(B4[0,:]),np.max(B5[0,:]),np.max(B6[0,:]),np.max(B7[0,:]),np.max(B8[0,:]),np.max(B9[0,:]),np.max(B10[0,:])])
    ymin = 0.0
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:]),np.max(B4[value,:]),np.max(B5[value,:]),np.max(B6[value,:]),np.max(B7[value,:]),np.max(B8[value,:]),np.max(B9[value,:]),np.max(B10[value,:])])
    ymax = 2.0
    axs[k].plot(B1[0,:],B1[value,:],'-',color=colors[0],label=r'$1 M_{\oplus},\ 10^{7.0}$')
    axs[k].plot(B2[0,:],B2[value,:],'-',color=colors[1],label=r'$1 M_{\oplus},\ 10^{7.5}$')
    axs[k].plot(B3[0,:],B3[value,:],'-',color=colors[2],label=r'$1 M_{\oplus},\ 10^{8.0}$')
    axs[k].plot(B4[0,:],B4[value,:],'-',color=colors[3],label=r'$1 M_{\oplus},\ 10^{8.5}$')
    axs[k].plot(B5[0,:],B5[value,:],'-',color=colors[4],label=r'$1 M_{\oplus},\ 10^{9.0}$')
    axs[k].plot(B6[0,:],B6[value,:],'--',color=colors[0],label=r'$10 M_{\oplus},\ 10^{7.0}$')
    axs[k].plot(B7[0,:],B7[value,:],'--',color=colors[1],label=r'$10 M_{\oplus},\ 10^{7.5}$')
    axs[k].plot(B8[0,:],B8[value,:],'--',color=colors[2],label=r'$10 M_{\oplus},\ 10^{8.0}$')
    axs[k].plot(B9[0,:],B9[value,:],'--',color=colors[3],label=r'$10 M_{\oplus},\ 10^{8.5}$')
    axs[k].plot(B10[0,:],B10[value,:],'--',color=colors[4],label=r'$10 M_{\oplus},\ 10^{9.0}$')
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} / M_{iron,tot}$', fontsize=fs_label)
    axs[k].legend(loc='best',ncols=2,fontsize=fs_legend)
    
    k = 1
    B1 = np.loadtxt('775_timeseries.txt')
    B2 = np.loadtxt('791_timeseries.txt')
    B3 = np.loadtxt('792_timeseries.txt')
    B1[value,:] /= (2.0/3.0*0.5)
    B2[value,:] /= (0.99*1.0/3.0)
    B3[value,:] /= (0.99*0.18)
    xmin = np.min([np.min(B1[0,:]),np.min(B2[0,:]),np.min(B3[0,:])])
    xmax = np.max([np.max(B1[0,:]),np.max(B2[0,:]),np.max(B3[0,:])])
    ymin = 0.0
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:])])
    ymax = 2.25
    axs[k].plot(B1[0,:],B1[value,:],'r-',label='50:50')
    axs[k].plot(B2[0,:],B2[value,:],'g-',label='33:67')
    axs[k].plot(B3[0,:],B3[value,:],'b-',label='18:82')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} / M_{iron,tot}$', fontsize=fs_label)
    axs[k].legend(loc='best',fontsize=fs_legend)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("demixing_resolution_ratio_study.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
