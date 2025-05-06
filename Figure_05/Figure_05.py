#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import matplotlib.lines as mlines

fs_label = 16
fs = 13
fs_legend = 12
fw = 'normal'

def main():
    B1 = np.loadtxt('311_timeseries.txt')
    B2 = np.loadtxt('312_timeseries.txt')
    B3 = np.loadtxt('513_timeseries.txt')
    B4 = np.loadtxt('514_timeseries.txt')
    B5 = np.loadtxt('518_timeseries.txt')
    B6 = np.loadtxt('511_timeseries.txt')
    B7 = np.loadtxt('512_timeseries.txt')
    B8 = np.loadtxt('519_timeseries.txt')
    B9 = np.loadtxt('5110_timeseries.txt')
    xmin = np.min([np.min(B1[0,:]),np.min(B2[0,:]),np.min(B3[0,:]),np.min(B4[0,:]),np.min(B5[0,:]),np.min(B6[0,:]),np.min(B7[0,:]),np.min(B8[0,:]),np.min(B9[0,:])])
    xmax = np.max([np.max(B1[0,:]),np.max(B2[0,:]),np.max(B3[0,:]),np.max(B4[0,:]),np.max(B5[0,:]),np.max(B6[0,:]),np.max(B7[0,:]),np.max(B8[0,:]),np.max(B9[0,:])])
    ymin = 0.0
    
    colormap = plt.cm.tab10
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    fig, axs = plt.subplots(1,3,figsize=(14*0.85, 4*0.85))
    lw = 1.75
    alpha = 0.7

    k, value = 0, 2
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:]),np.max(B4[value,:]),np.max(B5[value,:]),np.max(B6[value,:]),np.max(B7[value,:]),np.max(B8[value,:]),np.max(B9[value,:])])
    ymax = 6
    axs[k].plot(B2[0,:],B2[value,:],'-.',color=colors[0],label=r'Oblique, $25.8\times10^6$',linewidth=lw, alpha=alpha)
    axs[k].plot(B4[0,:],B4[value,:],'--',color=colors[0],label=r'Oblique, $10^8$',linewidth=lw, alpha=alpha)
    axs[k].plot(B7[0,:],B7[value,:],'-',color=colors[0],label=r'Oblique, $10^9$',linewidth=lw, alpha=alpha)
    axs[k].plot(B1[0,:],B1[value,:],'-.',color=colors[1],label=r'Head-on, $25.8\times10^6$',linewidth=lw, alpha=alpha)
    axs[k].plot(B3[0,:],B3[value,:],'--',color=colors[1],label=r'Head-on, $10^9$',linewidth=lw, alpha=alpha)
    axs[k].plot(B6[0,:],B6[value,:],'-',color=colors[1],label=r'Head-on, $10^9$',linewidth=lw, alpha=alpha)
    axs[k].plot(B9[0,:],B9[value,:],':',color=colors[1],label=r'Head-on, $2.1\times10^9$',linewidth=lw, alpha=alpha)
    axs[k].plot(B5[0,:],B5[value,:],'--',color=colors[2],label=r'Intermediate, $10^8$',linewidth=lw, alpha=alpha)
    axs[k].plot(B8[0,:],B8[value,:],'-',color=colors[2],label=r'Intermediate, $10^9$',linewidth=lw, alpha=alpha)
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 1, 4
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:]),np.max(B4[value,:]),np.max(B5[value,:]),np.max(B6[value,:]),np.max(B7[value,:]),np.max(B8[value,:]),np.max(B9[value,:])])
    ymax = 1.1
    axs[k].plot(B2[0,:],B2[value,:],'-.',color=colors[0],label='Oblique, 25.8M',linewidth=lw, alpha=alpha)
    axs[k].plot(B4[0,:],B4[value,:],'--',color=colors[0],label='Oblique, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B7[0,:],B7[value,:],'-',color=colors[0],label='Oblique, 1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B1[0,:],B1[value,:],'-.',color=colors[1],label='Head-on, 25.8M',linewidth=lw, alpha=alpha)
    axs[k].plot(B3[0,:],B3[value,:],'--',color=colors[1],label='Head-on, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B6[0,:],B6[value,:],'-',color=colors[1],label='Head-on, 1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B9[0,:],B9[value,:],':',color=colors[1],label='Head-on, 2.1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B5[0,:],B5[value,:],'--',color=colors[2],label='Intermediate, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B8[0,:],B8[value,:],'-',color=colors[2],label='Intermediate, 1B',linewidth=lw, alpha=alpha)
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 2, 7
    ymax = np.max([np.max(B1[value,:]),np.max(B2[value,:]),np.max(B3[value,:]),np.max(B4[value,:]),np.max(B5[value,:]),np.max(B6[value,:]),np.max(B7[value,:]),np.max(B8[value,:]),np.max(B9[value,:])])
    ymax = 2.0
    axs[k].plot(B2[0,:],B2[value,:],'-.',color=colors[0],label='Oblique, 25.8M',linewidth=lw, alpha=alpha)
    axs[k].plot(B4[0,:],B4[value,:],'--',color=colors[0],label='Oblique, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B7[0,:],B7[value,:],'-',color=colors[0],label='Oblique, 1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B1[0,:],B1[value,:],'-.',color=colors[1],label='Head-on, 25.8M',linewidth=lw, alpha=alpha)
    axs[k].plot(B3[0,:],B3[value,:],'--',color=colors[1],label='Head-on, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B6[0,:],B6[value,:],'-',color=colors[1],label='Head-on, 1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B9[0,:],B9[value,:],':',color=colors[1],label='Head-on, 2.1B',linewidth=lw, alpha=alpha)
    axs[k].plot(B5[0,:],B5[value,:],'--',color=colors[2],label='Intermediate, 100M',linewidth=lw, alpha=alpha)
    axs[k].plot(B8[0,:],B8[value,:],'-',color=colors[2],label='Intermediate, 1B',linewidth=lw, alpha=alpha)
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('Time [h]', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta\gamma} \ (M_{\oplus})$', fontsize=fs_label)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
            
    handles, labels = [], []
    h, l = axs[0].get_legend_handles_labels()
    handles.extend(h)
    labels.extend(l)
    
    empty_handle = mlines.Line2D([], [], color='white', label=" ")
    handles.insert(3, empty_handle)  # Insert at index 1 (adjust as needed)
    labels.insert(3, " ")  # A space makes it appear empty
    fig.legend(handles, labels, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.17),frameon=False,fontsize=fs_legend)

    plt.tight_layout()
    plt.savefig("jupiter_impacts_mixing_timeseries.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
