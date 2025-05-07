#!/usr/bin/env python3

import matplotlib.cm as cm
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 12
plt.style.use('seaborn-colorblind')
cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
markers = ["x","s", "o", "v", "^", "<", ">", "1", "2", "3"]

M_ear_cgs       = 5.972e27
M_ear_mks       = 5.972e24
M_Mars_cgs      = 6.417e26
M_Mars_mks      = 6.417e23
M_jup_cgs       = 1.898e30
M_jup_mks       = 1.898e27
M_Saturn_cgs    = 5.683e29
M_Saturn_mks    = 5.683e26
M_Uranus_cgs    = 8.6810e28
M_Uranus_mks    = 8.6810e25
M_Neptune_cgs   = 1.02413e29
M_Neptune_mks   = 1.02413e26 
M_sun_cgs       = 1.989e33
M_sun_mks       = 1.989e30

AU_cgs          = 1.49597870e13
AU_mks          = 1.49597870e11

#Mjup=c.M_jup_mks/c.M_sun_mks
#Rjup=c.R_jup_mks/c.AU_mks
#rho_jup=Mjup/(4./3.*np.pi*Rjup*Rjup*Rjup)

def disk_Gauss(r,Mdisk=40.*M_ear_mks/M_sun_mks,r0=7.,mu=1.,Sigma0=10.*AU_cgs*AU_cgs/M_sun_cgs):
    sigma=Mdisk/((2.*np.pi)**1.5 *mu*r0*r0*Sigma0)
    f0=np.sqrt(2.*np.pi)*sigma*Sigma0
    #Mdisk/(2.*np.pi *r0*r0 *mu)
    #mu=Mdisk/(2.*np.pi *r0*r0 *Sigma0)
    x=r/r0
    return f0/(np.sqrt(2.*np.pi)*sigma) *np.exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))

rc_Gauss=6.
mu_Gauss=1.
Mdisk_Gauss=40. *M_ear_mks/M_sun_mks
Sigma_Gauss=10. *AU_cgs*AU_cgs/M_sun_cgs
sig=Mdisk_Gauss/((2.*np.pi)**1.5 *mu_Gauss*rc_Gauss*rc_Gauss*Sigma_Gauss)
print(sig)

def Sigma(r):
    return disk_Gauss(r,Mdisk=Mdisk_Gauss,r0=rc_Gauss,mu=mu_Gauss,Sigma0=Sigma_Gauss)

def M_isolation(Sigma,r,Cfz=2.*np.sqrt(3.),Ms=1.):
    return 8.*(Cfz*np.pi*Sigma)**1.5 *r*r*r /(np.sqrt(3.*Ms))

def Hill(mu,ap):
    return (0.3333333333333333*mu)**(0.3333333333333333)*ap

def fx(x,XX):
    a_i = XX[0]
    M_i = XX[1]
    k   = XX[2]
    Cfz = XX[3]

    M_ii = M_isolation(Sigma(x),x,Cfz=Cfz,Ms=1.)
    a_ave = 0.5*(a_i+x)

    return np.abs(a_i-x) - k *Hill(M_i+M_ii, a_ave)

def bisection(XX, fx, epsl=1.e-3,x0=0.,x1=1.e1, nmax=1000, debug=False):
    n=0
    while True:
        x=0.5*(x0+x1)

        if fx(x0,XX)*fx(x,XX) < 0.:
            x1=x
        elif fx(x,XX)*fx(x1,XX) < 0.:
            x0=x
        else:
            print("Error, incorret x0 and x1: x0, x1, fx0, fx1, n = {:}, {:}, {:}, {:}, {:}".format(x0, x1, fx(x0,XX), fx(x,XX), n))
            x=np.nan
            break

        n+=1
        if abs(fx(x,XX))<epsl:
            break

        if n>nmax:
            print("Could not reach the equilibrium point, ", n)
            break
        
    if debug:
        print("fx converges by n = {:d}".format(n))
        
    return x


def main():
    rp_Jup=6.5

    ylim=(1.e-2,3e1)
        
    nx_plot=1
    ny_plot=1
    dx_plot=5
    dy_plot=3

    k=0
    fig = plt.figure(figsize=(nx_plot*dx_plot,ny_plot*dy_plot))
    ax = fig.add_subplot(ny_plot,nx_plot,k+1)

    ax.set_xlabel('Semi-Major Axis [au]')
    ax.set_ylabel('Mass [$M_\oplus$]')
    ax.set_yscale("log")
    ax.set_xlim(3.,10.)
    ax.set_ylim(ylim)

    x0=0.1
    x1=15.
    Nx=1000
    dx=(x1-x0)/Nx

    lines=[]
    colors=[cmap[0],cmap[1]]#["black", "red"]
    for jj, k_Hill in enumerate([5.,8.]):

        Cfz=0.5*k_Hill #2.*np.sqrt(3.)

        #>>>
        x,y=[],[]
        for r in np.arange(x0,x1,dx):
            M_iso = M_isolation(Sigma(r),r,Cfz=Cfz,Ms=1.) *M_sun_cgs/M_ear_cgs
            x.append(r)
            y.append(M_iso)
        line2, = ax.plot(x, y, color=colors[jj], linestyle="--", label="Isolation mass") #cmap[1])
        lines.append(line2)
        #ax.vlines(7.,ax.get_ylim()[0],ax.get_ylim()[1])


        x,y=[],[]
        x_error=[]

        r0=rp_Jup
        M0 = M_isolation(Sigma(r0),r0,Cfz=Cfz,Ms=1.)
        delta_axi = 0.5*k_Hill*Hill(M0,r0)
        x.append(r0)
        y.append(M0)
        x_error.append(delta_axi)
        for i in range(3):

            XX = [r0,M0,k_Hill,Cfz]
            r1 = bisection(XX, fx, epsl=1.e-3, x0=r0, x1=r0*2., nmax=1000)
            M1 = M_isolation(Sigma(r1),r1,Cfz=Cfz,Ms=1.)
            delta_axi = 0.5*k_Hill*Hill(M1,r1)
            print(r1,M1*M_sun_cgs/M_ear_cgs, (r1-r0)/Hill(M0+M1,0.5*(r1+r0)))

            x.append(r1)
            y.append(M1)
            x_error.append(delta_axi)

            M0=M1
            r0=r1

        r0=rp_Jup
        M0 = M_isolation(Sigma(r0),r0,Cfz=Cfz,Ms=1.)
        for i in range(6):

            XX = [r0,M0,k_Hill,Cfz]
            r1 = bisection(XX, fx, epsl=1.e-3, x0=0.5*r0, x1=r0, nmax=1000)
            M1 = M_isolation(Sigma(r1),r1,Cfz=Cfz,Ms=1.)
            delta_axi = 0.5*k_Hill*Hill(M1,r1)
            print(r1,M1*M_sun_cgs/M_ear_cgs)

            x.append(r1)
            y.append(M1)
            x_error.append(delta_axi)

            M0=M1
            r0=r1

        zip_lists = zip(x,y)
        zip_sort = sorted(zip_lists)
        x,y = zip(*zip_sort)
        xerr_l,xerr_r=[],[]
        for i in range(len(x)):
            if i==0:
                left  = 0.
            else:
                left  = 0.5*k_Hill*Hill(y[i]+y[i-1],0.5*(x[i]+x[i-1]))
            xerr_l.append(left)

            if i==len(x)-1:
                right = 0.
            else:
                right = 0.5*k_Hill*Hill(y[i]+y[i+1],0.5*(x[i]+x[i+1]))
            xerr_r.append(right)

        line3 = ax.errorbar(x,np.array(y)*M_sun_cgs/M_ear_cgs, xerr=[xerr_l,xerr_r], capsize=5, fmt='o', markersize=5, color=colors[jj],label="Embryos")#cmap[2])
        lines.append(line3)

    ax2 = ax.twinx()
    #ax.set_title("rc={:1.2f}, sig={:1.3f}, Mdisk={:2.1f}".format(rc_Gauss, sig, Mdisk_Gauss*c.M_sun_cgs/c.M_ear_cgs))
    ax2.set_xlabel('Semi-Major Axis [au]')
    ax2.set_ylabel('Surface density [g cm$^{-2}$]', rotation=-90, labelpad=15) #$F_\mathrm{col}$')#Fraction of collided planetesimals')
    #ax.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_ylim(ylim)
    ax2.set_xlim(3.,10.)

    x0=0.1
    x1=15.
    Nx=1000
    dx=(x1-x0)/Nx

    x=np.arange(x0,x1,dx)
    line1, = ax2.plot(x,[Sigma(r)*M_sun_cgs/AU_cgs**2. for r in x], color="black", linestyle=":", label="Gaussian ring")#cmap[0])
    k+=1


    lines[0].set_label("Isolation mass")
    lines[1].set_label("$k=5$")
    lines[3].set_label("$k=8$")
    lines = [lines[1], lines[3], lines[0], line1]  # Collect line objects
    labels = [line.get_label() for line in lines]  # Extract labels
    ax.legend(lines, labels, loc="lower center", fontsize=10)

    k+=1


    fig.savefig("disk.png",dpi=300,bbox_inches='tight')
    
if __name__ == '__main__':
    main()
