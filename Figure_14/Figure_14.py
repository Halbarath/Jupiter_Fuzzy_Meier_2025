#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd

fs_label = 16
fs = 13
fs_legend = 10
fw = 'normal'

def main():
    RockIce1 = np.array([[0.039375152, 0, 0.000726714],
        [0.028463577, 0, 0.001036104],
        [0.021210556, 0, 0.001409956],
        [0.016191813, 1.61E-06, 0.00111644]])
    RockIce10 = np.array([[0.600269378, 0.034567383, 0.090025125],
        [0.410080662, 0.00445708, 0.082388145],
        [0.363122413, 0.01433777, 0.099358368],
        [0.407098374, 0.086806304, 0.140481525]])
    IronRock1 = np.array([[0.094375044, 0.004741485, 0.022715267],
        [0.098138788, 0.00845836, 0.026262413],
        [0.090595012, 0.006534724, 0.028992199],
        [0.097949926, 0.00914909, 0.033515148]])
    IronRock10 = np.array([[1.151824885, 0.124599872, 0.265896441],
        [1.147571344, 0.248976443, 0.332024864],
        [1.069940043, 0.310596919, 0.376148855],
        [1.103586833, 0.419616887, 0.433643254]])
    IronIce1 = np.array([[0.042539306, 0.000202103, 0.009566497],
        [0.037329421, 0.000731311, 0.011948808],
        [0.036396813, 0.003099101, 0.014773732],
        [0.036168805, 0.003573887, 0.015995357]])
    IronIce10 = np.array([[0.662574454, 0.117472196, 0.178991781],
        [0.643356852, 0.174557362, 0.220444119],
        [0.560450694, 0.14504676, 0.214414745],
        [0.548940236, 0.149402407, 0.223600057]])
    
    xmin = 10**7
    xmax = 10**8.5
    ymin = 0.0
    
    colormap = plt.cm.tab10 #nipy_spectral, Set1,Paired    
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    
    fig, axs = plt.subplots(1,3,figsize=(14*0.85, 4*0.85))

    k, value = 0, 0
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 1, 1
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} \ (M_{\oplus})$', fontsize=fs_label)
    axs[k].legend(loc='best',ncols=2,fontsize=fs_legend)
    
    k, value = 2, 2
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta\gamma} \ (M_{\oplus})$', fontsize=fs_label)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("mixing_resolution_study_2mat_dense.pdf",bbox_inches='tight')
    
    RockIce1 /= 0.5
    RockIce10 /= 5.0
    IronRock1 /= 0.5
    IronRock10 /= 5.0
    IronIce1 /= 0.5
    IronIce10 /= 5.0
    
    fig, axs = plt.subplots(1,3,figsize=(14*0.85, 4*0.85))

    k, value = 0, 0
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    ymax = 0.25
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} / M_{core,tot}$', fontsize=fs_label)
    
    k, value = 1, 1
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    ymax = 0.1
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} / M_{core,tot}$', fontsize=fs_label)
    axs[k].legend(loc='upper left',fontsize=fs_legend)
    
    k, value = 2, 2
    ymax = np.max([np.max(RockIce1[:,value]),np.max(RockIce10[:,value]),np.max(IronRock1[:,value]),np.max(IronRock10[:,value]),np.max(IronIce1[:,value]),np.max(IronIce10[:,value])])
    ymax = 0.1
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce1[:,value],label=r'Rock/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],RockIce10[:,value],label=r'Rock/Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock1[:,value],label=r'Iron/Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronRock10[:,value],label=r'Iron/Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce1[:,value],label=r'Iron/Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],IronIce10[:,value],label=r'Iron/Ice $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta\gamma} / M_{core,tot}$', fontsize=fs_label)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("mixing_resolution_study_2mat_dense_scaled.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
