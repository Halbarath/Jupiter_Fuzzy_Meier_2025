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
    cases=['01_fully_segregated.std.data.txt','02_central_sphere_even_odd.std.data.txt','03_central_sphere_random.std.data.txt','04_full_sphere_random.std.data.txt','05_fully_segregated_moved.std.data.txt','06_central_sphere_even_odd_moved.std.data.txt','07_central_sphere_random_moved.std.data.txt','08_full_sphere_random_moved.std.data.txt']
    cases_label=['Segregated core','Core even/odd mixed','Core random mixed','Sphere random mixed','Segregated Core, noise','Core even/odd mixed, noise','Core random mixed, noise','Sphere random mixed, noise']
    cases_label=['Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7','Case 8']
    resolutions=['01_1M','02_3.16M','03_10M','04_31.6M','05_100M', '06_316M', '07_1B']
    resolutions_data=[1e6, 3.16e6, 10e6, 31.6e6, 100e6, 316e6, 1e9]
    
    fig, axs = plt.subplots(3,3,figsize=(14*0.85, 10*0.85))
    colormap = plt.cm.tab10 #nipy_spectral, Set1,Paired    
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]

    k, l, value = 0, 0, 1
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\alpha} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 0, 1, 2
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\beta} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 0, 2, 3
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\gamma} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 1, 0, 4
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\delta,0.2} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 1, 1, 5
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\delta,0.4} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 1, 2, 6
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\delta,0.6} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 2, 0, 7
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\beta\gamma} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 2, 1, 8
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("linear")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=-0.1,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel(r'$M_{mix,\epsilon} \ (M_{\oplus})$', fontsize=fs_label)
    
    k, l, value = 2, 2, 9
    maxval = 0
    for i, case in enumerate(cases):
        data = []
        for resolution in resolutions:
            A = np.loadtxt('{}/{}'.format(resolution,case))
            data.append(A)
        data = np.array(data)
        plotdata = data[:,value]
        if i < 4:
            axs[k,l].plot(resolutions_data,plotdata,color=colors[i],label=cases_label[i])
        else:
            axs[k,l].plot(resolutions_data,plotdata,'--',color=colors[i],label=cases_label[i])
        maxval = np.max([maxval , np.max(plotdata)])
    axs[k,l].set_xscale("log")
    axs[k,l].set_yscale("log")
    axs[k,l].set_xlim(xmin=np.min(resolutions_data),xmax=np.max(resolutions_data))
    axs[k,l].set_ylim(ymin=10,ymax=1.1*maxval)
    axs[k,l].set_xlabel('# Particles', fontsize=fs_label)
    axs[k,l].set_ylabel('$A_I\ (R_{\oplus}^2)$', fontsize=fs_label)
    axs[k,l].legend(loc='best',ncols=2,fontsize=fs_legend)

    for ax in axs.flatten():
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
            tick.label1.set_fontweight(fw)
    plt.tight_layout()
    plt.savefig("mixing_resolution_study.pdf",bbox_inches='tight')

if __name__ == '__main__':
    main()
