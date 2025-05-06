#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd

fs_label = 16
fs = 13
fs_legend = 12
fw = 'normal'

def main():
    ValuesRock1 = np.array([[0.645155865, 0.157131681, 0.083871187],
        [0.716983788, 0.235274781, 0.136296385],
        [0.792670029, 0.291068948, 0.165763905],
        [0.894526417, 0.344376411, 0.209856586]])
    ValuesRock10 = np.array([[5.709150396, 2.134165046, 0.891512833],
        [5.894065817, 2.894000563, 1.377435145],
        [6.059736508, 3.639149313, 1.685611769],
        [6.644959034, 4.259164139, 2.063188611]])
    ValuesIce1 = np.array([[0.530288226, 0.14051745, 0.069333538],
        [0.631938246, 0.182414056, 0.094263689],
        [0.651875403, 0.247197199, 0.131685812],
        [0.6591228, 0.291847286, 0.157434104]])
    ValuesIce10 = np.array([[5.227090879, 1.602459414, 0.744019438],
        [5.406548027, 2.27499856, 1.07786601],
        [5.739712682, 2.874283755, 1.42341497],
        [6.51606576, 3.356420624, 1.810987935]])
    ValuesIron1 = np.array([[0.598725053, 0.245780388, 0.126652362],
        [0.64180853, 0.342893375, 0.182723797],
        [0.665150818, 0.405431541, 0.213632715],
        [0.779290187, 0.465935396, 0.264419434]])
    ValuesIron10 = np.array([[5.61091894, 2.029550287, 0.904103773],
        [5.951464312, 2.685908806, 1.264396197],
        [6.123808172, 3.750248945, 1.881524754],
        [7.073454712, 4.360174929, 2.549513609]])
    
    xmin = 10**7
    xmax = 10**8.5
    ymin = 0.0
    
    colormap = plt.cm.tab10   
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]   
    
    fig, axs = plt.subplots(1,3,figsize=(14*0.85, 4*0.85))

    k, value = 0, 0
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} \ (M_{\oplus})$', fontsize=fs_label)
    axs[k].legend(loc='best',fontsize=fs_legend)
    
    k, value = 1, 1
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, value = 2, 2
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
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
    plt.savefig("mixing_resolution_study_pure_dense.pdf",bbox_inches='tight')
    
    ValuesRock1 /= 0.5
    ValuesRock10 /= 5.0
    ValuesIce1 /= 0.5
    ValuesIce10 /= 5.0
    ValuesIron1 /= 0.5
    ValuesIron10 /= 5.0
    
    fig, axs = plt.subplots(1,3,figsize=(14, 4))

    k, value = 0, 0
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    ymax = 2.0
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\beta} / M_{core,tot}$', fontsize=fs_label)
    axs[k].legend(loc='lower right',fontsize=fs_legend)
    
    k, value = 1, 1
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    ymax = 1.0
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
    for i,j in enumerate(axs[k].lines):
        j.set_color(colors[i])
    axs[k].set_xscale("log")
    axs[k].set_yscale("linear")
    axs[k].set_xlim(xmin=xmin,xmax=xmax)
    axs[k].set_ylim(ymin=ymin,ymax=ymax)
    axs[k].set_xlabel('# Particles', fontsize=fs_label)
    axs[k].set_ylabel(r'$M_{mix,\delta,0.2} / M_{core,tot}$', fontsize=fs_label)
    
    k, value = 2, 2
    ymax = np.max([np.max(ValuesRock1[:,value]),np.max(ValuesRock10[:,value]),np.max(ValuesIce1[:,value]),np.max(ValuesIce10[:,value]),np.max(ValuesIron1[:,value]),np.max(ValuesIron10[:,value])])
    ymax = 0.6
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock1[:,value],label=r'Rock $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesRock10[:,value],label=r'Rock $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce1[:,value],label=r'Ice $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIce10[:,value],label=r'Ice $10 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron1[:,value],label=r'Iron $1 M_{\oplus}$')
    axs[k].plot([10**7,10**7.5,10**8,10**8.5],ValuesIron10[:,value],label=r'Iron $10 M_{\oplus}$')
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
    plt.savefig("mixing_resolution_study_pure_dense_scaled.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
