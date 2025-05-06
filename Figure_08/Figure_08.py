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
    A = np.loadtxt('target.model')
    B = np.loadtxt('511.00100.profile.txt')
    C = np.loadtxt('511.02000.profile.txt')
    D = np.loadtxt('512.00110.profile.txt')
    E = np.loadtxt('512.02000.profile.txt')

    F = np.array([[0.0004962779156327543, 4.4088397790055245],
[0.04962779156327543, 4.4088397790055245],
[0.09975186104218363, 4.359116022099447],
[0.14590570719602977, 4.2707182320441985],
[0.19404466501240694, 4.143646408839778],
[0.2456575682382134, 3.9723756906077345],
[0.29727047146401986, 3.779005524861878],
[0.3543424317617866, 3.5303867403314917],
[0.41786600496277915, 3.2265193370165743],
[0.4674937965260546, 2.917127071823204],
[0.5196029776674937, 2.5635359116022096],
[0.5751861042183622, 2.2209944751381214],
[0.630272952853598, 1.889502762430939],
[0.6883374689826303, 1.6022099447513811],
[0.7245657568238213, 1.4143646408839778],
[0.747394540942928, 1.287292817679558],
[0.798014888337469, 1.0441988950276242],
[0.8689826302729529, 0.712707182320442],
[0.9245657568238214, 0.38121546961325964],
[0.9687344913151364, 0.143646408839779],
[1.0, 0.0]])

    QQ = pd.read_csv('CMS19+HG23.csv', sep=',')
    G = QQ[['R_RTOT', 'RHO_GCC']].replace({'D': 'E'}, regex=True).apply(pd.to_numeric).to_numpy()
    
    QQ = np.loadtxt('out10967_layers_mrt_with_header.txt',skiprows=30)
    H = np.array([QQ[:,1], QQ[:,6]]).T
   
    labelA='Initial'
    labelB='Head-on, Intermediate'
    labelC='Head-on, Final'
    labelD='Oblique, Intermediate'
    labelE='Oblique, Final'
    labelG='Howard et al. 2023'
    labelF='Militzer et al. 2022'
    labelH='Militzer & Hubbard 2024'
    
    colormap = plt.cm.tab10 #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 1, 10)]
    lw = 1.75
    alpha = 0.7

    fig, ax = plt.subplots(1,1,figsize=(7*0.85, 5*0.85))
    
    ax.plot(A[:,0],A[:,1]*0.368477421278,color=colors[0],label=labelA,linewidth=lw, alpha=alpha)
    ax.plot(D[:,0],D[:,2],color=colors[1],label=labelD,linewidth=lw, alpha=alpha)
    ax.plot(E[:,0],E[:,2],color=colors[2],label=labelE,linewidth=lw, alpha=alpha)
    ax.plot(B[:,0],B[:,2],color=colors[3],label=labelB,linewidth=lw, alpha=alpha)
    ax.plot(C[:,0],C[:,2],color=colors[4],label=labelC,linewidth=lw, alpha=alpha)
    ax.plot(G[:,0]*11.209,G[:,1],'--',color='k',label=labelG,linewidth=lw, alpha=alpha)
    ax.plot(H[:,0]*11.209,H[:,1],'-.',color='k',label=labelH,linewidth=lw, alpha=alpha)
    ax.plot(F[:,0]*11.209,F[:,1],':',color='k',label=labelF,linewidth=lw, alpha=alpha)
    ax.set_yscale("linear")
    ax.set_xlim(xmin=0,xmax=12)
    ax.set_ylim(ymin=0,ymax=20)
    ax.set_xlabel(r'Radius ($R_{\oplus}$)', fontsize=fs_label)
    ax.set_ylabel(r'Density ($gcm^{3}$)', fontsize=fs_label)
    ax.legend(loc='best',fontsize=fs_legend,frameon=False)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontweight(fw)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontweight(fw)

    plt.tight_layout()
    plt.savefig("jupiter_density_profiles.pdf",bbox_inches='tight')
    
if __name__ == '__main__':
    main()
