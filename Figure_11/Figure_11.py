#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd

fs_label = 16
fs = 13
fs_legend = 9
fw = 'normal'

def main():
    Values711_4 = np.array([np.ones(4) * 1.0,[1.29E+00,1.43E+00,1.59E+00,1.79E+00]])
    Values715_8 = np.array([np.ones(4) * 1.0,[1.14E+00,1.18E+00,1.21E+00,1.33E+00]])
    Values721_4 = np.array([np.ones(4) * 1.0,[1.06E+00,1.26E+00,1.30E+00,1.32E+00]])
    Values725_8 = np.array([np.ones(4) * 1.0,[1.05E+00,1.08E+00,1.15E+00,1.30E+00]])
    Values731_4 = np.array([np.ones(4) * 1.842461683,[7.88E-02,5.69E-02,4.24E-02,3.24E-02]])
    Values735_8 = np.array([np.ones(4) * 1.79184914,[1.20E-01,8.20E-02,7.26E-02,8.14E-02]])
    Values741_4 = np.array([np.ones(4) * 1.0,[1.20E+00,1.28E+00,1.33E+00,1.56E+00]])
    Values745_8 = np.array([np.ones(4) * 1.0,[1.12E+00,1.19E+00,1.22E+00,1.41E+00]])
    Values751_4 = np.array([np.ones(4) * 2.114794959,[1.89E-01,1.96E-01,1.81E-01,1.96E-01]])
    Values755_8 = np.array([np.ones(4) * 1.854889906,[2.30E-01,2.30E-01,2.14E-01,2.21E-01]])
    Values761_4 = np.array([np.ones(4) * 3.996823013,[8.51E-02,7.47E-02,7.28E-02,7.23E-02]])
    Values765_8 = np.array([np.ones(4) * 3.403352267,[1.33E-01,1.29E-01,1.12E-01,1.10E-01]])
    ValuesJupiterHeadon = np.array([np.ones(3) * 2.392734121,[9.08E-02,6.64E-02,4.28E-02]])
    ValuesJupiterOblique = np.array([np.ones(3) * 2.392734121,[1.98E-01,1.81E-01,1.39E-01]])
    ValuesJupiterCoreHit = np.array([np.ones(1) * 2.392734121,[1.06E-01]])

    fig, ax = plt.subplots(1,1,figsize=(7*0.85, 5*0.85))
    colormap = plt.cm.tab10   
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    
    ax.scatter(Values711_4[0,:]*0.95,Values711_4[1,:],color=colors[0],marker='o',label=r'Rock/Rock $1M_\oplus$')
    ax.scatter(Values715_8[0,:]*0.97,Values715_8[1,:],color=colors[0],marker='x',label=r'Rock/Rock $10M_\oplus$')
    ax.scatter(Values721_4[0,:]*0.99,Values721_4[1,:],color=colors[1],marker='o',label=r'Ice/Ice $1M_\oplus$')
    ax.scatter(Values725_8[0,:]*1.01,Values725_8[1,:],color=colors[1],marker='x',label=r'Ice/Ice $10M_\oplus$')
    ax.scatter(Values741_4[0,:]*1.03,Values741_4[1,:],color=colors[2],marker='o',label=r'Iron/Iron $1M_\oplus$')
    ax.scatter(Values745_8[0,:]*1.05,Values745_8[1,:],color=colors[2],marker='x',label=r'Iron/Iron $10M_\oplus$')
    ax.scatter(Values731_4[0,:],Values731_4[1,:],color=colors[3],marker='o',label=r'Rock/Ice $1M_\oplus$')
    ax.scatter(Values735_8[0,:],Values735_8[1,:],color=colors[3],marker='x',label=r'Rock/Ice $10M_\oplus$')
    ax.scatter(Values751_4[0,:],Values751_4[1,:],color=colors[4],marker='o',label=r'Iron/Rock $1M_\oplus$')
    ax.scatter(Values755_8[0,:],Values755_8[1,:],color=colors[4],marker='x',label=r'Iron/Rock $10M_\oplus$')
    ax.scatter(Values761_4[0,:],Values761_4[1,:],color=colors[5],marker='o',label=r'Iron/Ice $1M_\oplus$')
    ax.scatter(Values765_8[0,:],Values765_8[1,:],color=colors[5],marker='x',label=r'Iron/Ice $10M_\oplus$')
    ax.scatter(ValuesJupiterHeadon[0,:],ValuesJupiterHeadon[1,:],color=colors[6],marker='o',label=r'Jupiter Head-on')
    ax.scatter(ValuesJupiterCoreHit[0,:],ValuesJupiterCoreHit[1,:],color=colors[6],marker='x',label=r'Jupiter Intermediate')
    ax.scatter(ValuesJupiterOblique[0,:],ValuesJupiterOblique[1,:],color=colors[6],marker='+',label=r'Jupiter Oblique')
    ax.set_yscale("log")
    ax.set_xlim(xmin=0.5,xmax=4.5)
    ax.set_ylim(ymin=3e-2,ymax=2)
    ax.set_xlabel('CMB density contrast', fontsize=fs_label)
    ax.set_ylabel(r'$M_{mix,\beta} / M_{core,tot}$', fontsize=fs_label)
    ax.legend(loc='best',ncols=2,fontsize=fs_legend,frameon=False)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontweight(fw)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("Mbeta_vs_CMB_density_contrast.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
