#!/usr/bin/env python3

import os
import glob as glob
import matplotlib.pyplot as plt
import numpy as np
import math as m
import importlib
import pandas as pd

dx_plot=5
dy_plot=4
mass_star=1.

plt.rcParams["font.size"] = 12
plt.style.use('seaborn-v0_8-colorblind')
cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']

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

class Particle:
    """
    Class to store data of particle
    """
    #def __init__(self, mass, x, y, z, vx, vy, vz):
    #    self.mass = mass
    #    self.x  = x
    #    self.y  = y
    #    self.z  = z
    #    self.vx = vx
    #    self.vy = vy
    #    self.vz = vz

    def __init__(self, id, time, mass, radius, x, y, z, vx, vy, vz, massSP, exts=False):
        self.id = id
        self.time = time
        self.mass = mass
        self.radius = radius
        self.x  = x
        self.y  = y
        self.z  = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.massSP = massSP
        if exts!=False:
            self.exts=exts
        
    def pos(self, target=None):
        if target is None:
            return [self.x, self.y, self.z]
        else:
            return [self.x-target.x, self.y-target.y, self.z-target.z]

    def pos2(self, target=None):
        x, y, z = self.pos(target)
        return x**2 + y**2 + z**2
        
    def vel(self, target=None):
        if target is None:
            return [self.vx, self.vy, self.vz]
        else:
            return [self.vx-target.vx, self.vy-target.vy, self.vz-target.vz]

    def vel2(self, target=None):
        vx, vy, vz = self.vel(target)
        return vx**2 + vy**2 + vz**2

    def r_dot_v(self, target=None):
        if target is None:
            return self.x*self.vx + self.y*self.vy + self.z*self.vz
        else:
            return (self.x-target.x)*(self.vx-target.vx) + (self.y-target.y)*(self.vy-target.vy) + (self.z-target.z)*(self.vz-target.vz)

    def r_x_v(self, target=None):
        if target is None:
            return [self.y*self.vz - self.z*self.vy,\
                    self.z*self.vx - self.x*self.vz,\
                    self.x*self.vy - self.y*self.vx\
                    ]
        else:
            return  [(self.y-target.y)*(self.vz-target.vz) - (self.z-target.z)*(self.vy-target.vy),\
                     (self.z-target.z)*(self.vx-target.vx) - (self.x-target.x)*(self.vz-target.vz),\
                     (self.x-target.x)*(self.vy-target.vy) - (self.y-target.y)*(self.vx-target.vx)\
                    ]


    def getOrbitalElement(self, mu=1., target=None):
        if target is None:
            r = [self.x, self.y, self.z]
        else:
            r = [self.x-target.x, self.y-target.y, self.z-target.z]

        rn = m.sqrt(self.pos2(target=target))
        v2 = self.vel2(target=target)
        rdv = self.r_dot_v(target=target)
        rxv = self.r_x_v(target=target)
        hn = m.sqrt(rxv[0]*rxv[0] +rxv[1]*rxv[1] + rxv[2]*rxv[2])
        dr_dt=rdv/rn
        
        self.axi  = 1./(2./rn - v2/mu)
        self.ecc = m.sqrt( (1.-rn/self.axi)**2 + (rdv)**2/(mu*self.axi) )
        self.inc = m.atan2(m.sqrt(rxv[0]**2+rxv[1]**2), rxv[2])

        #--- Longitude of Ascending node
        #if rxv[2]==0.e0 :
        if self.inc==0.e0 :
            self.loa = 0.e0
        elif rxv[2]>0 :
            self.loa = m.atan2(rxv[0],-rxv[1])
        else:
            self.loa = m.atan2(-rxv[0],rxv[1])

        #--- theta
        if self.inc==0:
            theta=m.atan2( r[1] , r[0] )

        else:
            sin_theta   =   r[2]/(rn*m.sin(self.inc))
            cos_theta   =   ( self.x/rn +m.sin(self.loa)*sin_theta*m.cos(self.inc) )/m.cos(self.loa)
            theta       =   m.atan2( sin_theta , cos_theta )

        #--- longitude of pericenter and true anomaly
        if self.ecc==0:
            self.lop=0.e0
            self.tra=theta
        else:
            self.tra=m.atan2(self.axi*(1.e0-self.ecc*self.ecc)/hn *dr_dt, self.axi*(1.e0-self.ecc*self.ecc)/rn -1.e0 )
            omega=theta-self.tra
            omega=omega+10.e0*m.pi
            self.lop=omega % (2.e0*m.pi)

        return [self.axi, self.ecc, self.inc, self.loa, self.lop, self.tra]

    def setInitialCondition(self,ptcl0):
        self.axi0=ptcl0.axi
        self.ecc0=ptcl0.ecc
        self.inc0=ptcl0.inc
        self.mass0=ptcl0.mass
        self.radius0=ptcl0.radius
        
        return 0

    def getSuperParticleMass(self, n, Sigma):
        self.Msp = Sigma(self.axi0)/n(self.axi0)
        return self.Msp
    

class LogCollision:
    """
    Class to store data of particle
    """

    def __init__(self, id0, id1, time, m0, m1, \
                x0, y0, z0, vx0, vy0, vz0, \
                x1, y1, z1, vx1, vy1, vz1, \
                exts \
                #xs, ys, zs, vxs, vys, vzs \
        ):
        self.id0 = id0
        self.id1 = id1
        self.time = time
        self.m0 = m0
        self.m1 = m1
        
        self.x0  = x0
        self.y0  = y0
        self.z0  = z0
        self.vx0 = vx0
        self.vy0 = vy0
        self.vz0 = vz0
        
        self.x1  = x1
        self.y1  = y1
        self.z1  = z1
        self.vx1 = vx1
        self.vy1 = vy1
        self.vz1 = vz1
    
        self.exts = exts

        #>>>
        self.ptcl0 = Particle(0, self.time, self.m0, np.nan, self.x0, self.y0, self.z0, self.vx0, self.vy0, self.vz0, np.nan, exts=False)
        self.ptcl1 = Particle(1, self.time, self.m1, np.nan, self.x1, self.y1, self.z1, self.vx1, self.vy1, self.vz1, np.nan, exts=False)

    def pos(self, id_impactor):
        if id_impactor==0:
            # impactor _0
            x0=self.x0
            y0=self.y0
            z0=self.z0

            # target _1
            x1=self.x1
            y1=self.y1
            z1=self.z1
        elif id_impactor==1:
            # impactor _0
            x0=self.x1
            y0=self.y1
            z0=self.z1

            # target _1
            x1=self.x0
            y1=self.y0
            z1=self.z0
        else:
            print("Error:")

        return [x0-x1,y0-y1,z0-z1]

    def pos2(self, id_impactor):
        x, y, z = self.pos(id_impactor)
        return x**2 + y**2 + z**2
        
    def vel(self, id_impactor):
        if id_impactor==0:
            # impactor _0
            vx0=self.vx0
            vy0=self.vy0
            vz0=self.vz0

            # target _1
            vx1=self.vx1
            vy1=self.vy1
            vz1=self.vz1
        elif id_impactor==1:
            # impactor _0
            vx0=self.vx1
            vy0=self.vy1
            vz0=self.vz1

            # target _1
            vx1=self.vx0
            vy1=self.vy0
            vz1=self.vz0
        else:
            print("Error:")

        return [vx0-vx1,vy0-vy1,vz0-vz1]

    def vel2(self, id_impactor):
        vx, vy, vz = self.vel(id_impactor)
        return vx**2 + vy**2 + vz**2

    def r_dot_v(self, id_impactor):
        if id_impactor==0:
            # impactor _0
            x0=self.x0
            y0=self.y0
            z0=self.z0
            vx0=self.vx0
            vy0=self.vy0
            vz0=self.vz0

            # target _1
            x1=self.x1
            y1=self.y1
            z1=self.z1
            vx1=self.vx1
            vy1=self.vy1
            vz1=self.vz1
            
        elif id_impactor==1:
            # impactor _0
            x0=self.x1
            y0=self.y1
            z0=self.z1
            vx0=self.vx1
            vy0=self.vy1
            vz0=self.vz1

            # target _1
            x1=self.x0
            y1=self.y0
            z1=self.z0
            vx1=self.vx0
            vy1=self.vy0
            vz1=self.vz0

        else:
            print("Error:")

        return (x0-x1)*(vx0-vx1) + (y0-y1)*(vy0-vy1) + (z0-z1)*(vz0-vz1)

    def r_x_v(self, id_impactor, target=None):
        if id_impactor==0:
            # impactor _0
            x0=self.x0
            y0=self.y0
            z0=self.z0
            vx0=self.vx0
            vy0=self.vy0
            vz0=self.vz0

            # target _1
            x1=self.x1
            y1=self.y1
            z1=self.z1
            vx1=self.vx1
            vy1=self.vy1
            vz1=self.vz1
            
        elif id_impactor==1:
            # impactor _0
            x0=self.x1
            y0=self.y1
            z0=self.z1
            vx0=self.vx1
            vy0=self.vy1
            vz0=self.vz1

            # target _1
            x1=self.x0
            y1=self.y0
            z1=self.z0
            vx1=self.vx0
            vy1=self.vy0
            vz1=self.vz0

        else:
            print("Error:")


        return  [   (y0-y1)*(vz0-vz1) - (z0-z1)*(vy0-vy1), \
                    (z0-z1)*(vx0-vx1) - (x0-x1)*(vz0-vz1), \
                    (x0-x1)*(vy0-vy1) - (y0-y1)*(vx0-vx1) \
                ]

    def getCollisionProperty(self):
        if self.m0 > self.m1:
            id1=0 # target
            id0=1 # impactor
            self.mass_target   = self.m0
            self.mass_impactor = self.m1
        else:
            id1=1 # target
            id0=0 # impactor
            self.mass_target   = self.m1
            self.mass_impactor = self.m0

        self.rcol=self.pos(id0)
        self.rncol=np.linalg.norm(self.rcol)
        self.vcol=self.vel(id0)
        self.vncol=np.linalg.norm(self.vcol)
        self.vesc = np.sqrt(2.*(self.m0+self.m1)/self.rncol)
        self.vesc_type2 = np.sqrt(2.*(np.max([self.m0,self.m1]))/self.rncol)
        self.cosAng = -self.r_dot_v(id0)/(self.rncol*self.vncol)
        return 0


def Load_MyLogcolFile(filename, l0=2, status="default"):
    pp = []
    with open(filename) as f:
        l=-1
        for line in f:
            l+=1
            #>> header
            if l<l0:
                continue

            #>> header
            part = np.array([p.strip() for p in line.split()])
            if len(part)<=17:
                ptcl = LogCollision( \
                        # id0, id1
                        int(part[0]),\
                        int(part[1]),\
                        # time, m0, ma
                        float(part[2]),\
                        float(part[3]),\
                        float(part[4]),\
                        # x,y,z,vx,vy,vz
                        float(part[5]),\
                        float(part[6]),\
                        float(part[7]),\
                        float(part[8]),\
                        float(part[9]),\
                        float(part[10]), \
                        # x,y,z,vx,vy,vz
                        float(part[11]),\
                        float(part[12]),\
                        float(part[13]),\
                        float(part[14]),\
                        float(part[15]),\
                        float(part[16]) \
                        # exts
                        )
            else:
                ptcl = LogCollision( \
                        # id0, id1
                        int(part[0]),\
                        int(part[1]),\
                        # time, m0, ma
                        float(part[2]),\
                        float(part[3]),\
                        float(part[4]),\
                        # x,y,z,vx,vy,vz
                        float(part[5]),\
                        float(part[6]),\
                        float(part[7]),\
                        float(part[8]),\
                        float(part[9]),\
                        float(part[10]), \
                        # x,y,z,vx,vy,vz
                        float(part[11]),\
                        float(part[12]),\
                        float(part[13]),\
                        float(part[14]),\
                        float(part[15]),\
                        float(part[16]), \
                        # exts
                        [ x for x in part[17:len(part)]] \
                        )                

            if status=="default":
                pp.append(ptcl)
            elif status=="collision":
                if ptcl.id0>=0 and ptcl.id1>=0:
                    pp.append(ptcl)

    return pp

def getOrbitalElement( r, v, mu=1.):
    rx, ry, rz = r[0], r[1], r[2]
    vx, vy, vz = v[0], v[1], v[2]
    r = m.sqrt(rx*rx +ry*ry +rz*rz)
    v2 = vx*vx +vy*vy +vz*vz
    rv = rx*vx + ry*vy + rz*vz
    rxv = [ry*vz - rz*vy,\
           rz*vx - rx*vz,\
           rx*vy - ry*vx \
    ]
    
    ax  = 1./(2./r - v2/mu)
    ecc = m.sqrt( (1.-r/ax)**2 + (rv)**2/(mu*ax) )
    inc = m.atan2(m.sqrt(rxv[0]**2+rxv[1]**2), rxv[2])
    return ax,ecc,inc

def main():
    Nsim=520

    WORK_DIRS=[]#['Incp0','Incp1e-4','Incp1e-3','Incp1e-2','Incp1e-1']
    WORK_DIRS.append("./data/psIncp/Incp0")
    WORK_DIRS.append("./data/psIncp/Incp1e-4")
    WORK_DIRS.append("./data/psIncp/Incp1e-2")
    WORK_DIRS.append("./data/psdAxi/kHill5")

    vals_inc=[0.,1.e-4,1.e-2,1.e-2]
    vals_kHill=[8.,8.,8.,5.]

    datasets=[]
    for j,WORK_DIR in enumerate(WORK_DIRS):

        # log collision
        cols=[]
        filenames = sorted(glob.glob(WORK_DIR+"/Log/ID?????logcol.data"))
        for filename in filenames:
            pp = Load_MyLogcolFile(filename, l0=2, status="default")
            cols.append(pp)

        datasets.append(np.array([np.array(cols), np.array(cols)]).T)

    rhop=0.125 * 1./M_sun_cgs / (1./AU_cgs)**3.

    Mmin=0.5*M_ear_cgs/M_sun_cgs
    Mtg_th=150.*M_ear_cgs/M_sun_cgs

    rho_J=0.125* AU_cgs*AU_cgs*AU_cgs/M_sun_cgs
    rho_E=5.5* AU_cgs*AU_cgs*AU_cgs/M_sun_cgs

    pp_psInc=[]
    for j,dds in enumerate(datasets):

        tcol=[]
        mpcol=[]
        micol=[]
        vcol=[]
        vesc=[]
        cosAng=[]
        Ncols=[]
        a0imp=[]
        bcol=[]
        rcol=[]
        eccp=[]
        
        for dd in dds:
            pp=dd[1]
            Ncol=0
            for ptcl in pp:
                if (ptcl.id0==1 or ptcl.id0==1):
                    Mimpactor=min(ptcl.m0,ptcl.m1)
                    Mtarget=max(ptcl.m0,ptcl.m1)
                    if Mimpactor>Mmin:
                        ptcl.getCollisionProperty()
                        tcol.append(ptcl.time)
                        mpcol.append(max(ptcl.m0,ptcl.m1))
                        micol.append(Mimpactor)
                        cosAng.append(ptcl.cosAng)
                        if ptcl.m0 > ptcl.m1:
                            rrel=ptcl.rncol
                            r0=(3.*ptcl.m0/(4.*np.pi*rho_J))**(1./3.)
                            r1=(3.*ptcl.m1/(4.*np.pi*rho_E))**(1./3.)
                            rcol.append(rrel/(r0+r1))

                            rr=[ptcl.x0, ptcl.y0, ptcl.z0]
                            vv=[ptcl.vx0, ptcl.vy0, ptcl.vz0]
                            ax,ecc,inc=getOrbitalElement( rr, vv, mu=1.)
                            eccp.append(ecc)

                        else:
                            rrel=ptcl.rncol
                            r0=(3.*ptcl.m0/(4.*np.pi*rho_E))**(1./3.)
                            r1=(3.*ptcl.m1/(4.*np.pi*rho_J))**(1./3.)
                            rcol.append(rrel/(r0+r1))

                            rr=[ptcl.x1, ptcl.y1, ptcl.z1]
                            vv=[ptcl.vx1, ptcl.vy1, ptcl.vz1]

                            ax,ecc,inc=getOrbitalElement( rr, vv, mu=1.)
                            eccp.append(ecc)

                        #vcol.append(ptcl.vncol)
                        vcol_cor=np.sqrt(ptcl.vncol**2. -2.*(ptcl.m0+ptcl.m1)/rrel +2.*(ptcl.m0+ptcl.m1)/(r0+r1))
                            
                        vcol.append(vcol_cor) #ptcl.vncol)
                        #vcol.append(ptcl.vncol)

                        #vesc.append(ptcl.vesc)
                        vesc.append(np.sqrt(2.*(ptcl.m0+ptcl.m1)/(r0+r1)))


                        #>>
                        a0imp.append(np.nan) #max(ifile[ptcl.id0].axi,ifile[ptcl.id1].axi))
                    
                        r0 = (3.*max(ptcl.m0,ptcl.m1) / (4.*np.pi*rhop))**(1./3.)
                        bcol.append( r0 *np.sqrt(1.-ptcl.cosAng*ptcl.cosAng) )
                        Ncol+=1
            Ncols.append(Ncol)
        pp_psInc.append([np.array(tcol),np.array(mpcol),np.array(vcol),np.array(vesc),np.array(cosAng),np.array(micol),np.array(Ncols),np.array(a0imp),np.array(bcol),np.array(rcol),np.array(eccp)])

    IDs=[0,1,2,3]

    lw=1.5
    dx=0.05

    labels=[]
    labels.append("$\sin i_0=0$, $k=8$")
    labels.append("$\sin i_0=10^{-4}$, $k=8$")
    labels.append("$\sin i_0=10^{-2}$, $k=8$")
    labels.append("$\sin i_0=10^{-2}$, $k=5$")
    labels.append("$\sin i_0=10^{-4}$, $k=5$")
    labels.append("$\sin i_0=10^{-4}$, $k=10$")

    nx=1
    ny=1
    k=1
    dy_plot=4
    fig = plt.figure(figsize=(dx_plot*nx,dy_plot*ny))

    ax = fig.add_subplot(ny,nx,k)
    ax.set_xlabel(r"Collision Number")
    ax.set_ylabel(r"Probability") #$km/s$
    #ax.set_yscale('log')
    for j,id in enumerate(IDs):
        cps=pp_psInc[id]
        label=labels[id]#"$\sin i_0$ = {:1.1e}, $k$ = {:d}".format(vals_inc[id], int(vals_kHill[id]))
        #label="$\sin i_0$ = {:1.1e}".format(values[id])

        x=cps[6]
        weights = np.ones(len(x))/float(len(x)) #, density=True
        shift=dx*j
        ax.hist(np.array(x), bins=np.linspace(0,5,6), range=(0,5), histtype="step", weights=weights, align='left', label=label, color=cmap[j], linewidth=lw)
    #ax.legend(fontsize=10,loc="lower left", bbox_to_anchor=(0.05, 0.05))
    ax.legend(fontsize=10)
    #ax.hlines(0.5,0,4)
    k+=1    

    fig.savefig("CollisionNumber_psInc.png",dpi=300,bbox_inches='tight')


    nx=2
    ny=2
    k=1
    fig = plt.figure(figsize=(dx_plot*nx,dy_plot*ny))
    plt.subplots_adjust(wspace=0.25, hspace=0.25)
    x_title=0.05
    y_title=1

    Nbin=10
    k=1
    ax = fig.add_subplot(ny,nx,k)
    ax.set_xlabel(r"Impact Angle [${}^\circ$]")
    ax.set_ylabel(r"Probability") #$km/s$
    ax.set_yscale("log")
    ax.set_title("(a)",x=x_title,y=y_title)
    for j,id in enumerate(IDs):
        cps=pp_psInc[id]
        label=labels[id] #"$\sin i_0$ = {:1.1e}, $k$ = {:d}".format(vals_inc[id], vals_kHill[id])
        x=np.arccos(cps[4])*180/np.pi
        weights = np.ones(len(x))/Nsim #float(len(x)) #, density=True
        bins=np.linspace(0,90,Nbin)#+dx*j
        ax.hist(x, bins=bins, range=(0,90), weights=weights, histtype="step",label=label, color=cmap[j], linewidth=lw)
    #ax.legend()
    k+=1


    ax = fig.add_subplot(ny,nx,k)
    ax.set_xlabel(r"Impact Veloctiy / Escape Velocity")#v_{\rm esc}
    ax.set_ylabel(r"Probability") #$km/s$
    ax.set_yscale("log")
    ax.set_title("(b)",x=x_title,y=y_title)
    for j,id in enumerate(IDs):
        cps=pp_psInc[id]
        label=labels[id] #label="$\sin i_0$ = {:1.1e}".format(values[id])
        x=cps[2]/cps[3]
        weights = np.ones(len(x))/Nsim #float(len(x)) #, density=True
        ax.hist(x, bins=Nbin, range=(0.995,1.07), weights=weights, histtype="step",label=label, color=cmap[j], linewidth=lw)
    ax.legend()
    k+=1    


    ax = fig.add_subplot(ny,nx,k)
    ax.set_xlabel(r"Proto Jupiter's Mass at Collision [$M_\oplus$]")
    ax.set_ylabel(r"Probability") #$km/s$
    ax.set_yscale("log")
    ax.set_title("(c)",x=x_title,y=y_title)
    for j,id in enumerate(IDs):
        cps=pp_psInc[id]
        label=labels[id]#label="$\sin i_0$ = {:1.1e}".format(values[id])
        x=cps[1]*M_sun_mks/M_ear_mks
        weights = np.ones(len(x))/Nsim #float(len(x)) #, density=True
        ax.hist(x, bins=Nbin, range=(30, 320), weights=weights, histtype="step",label=label, color=cmap[j], linewidth=lw)
    #ax.legend()
    k+=1


    ax = fig.add_subplot(ny,nx,k)
    ax.set_xlabel(r"Impactor's mass [$M_\oplus$]")
    ax.set_ylabel(r"Probability") #$km/s$
    ax.set_title("(d)",x=x_title,y=y_title)
    ax.set_yscale("log")
    for j,id in enumerate(IDs):
        cps=pp_psInc[id]
        label=labels[id] #label="$i_0$ = {:1.1e}".format(values[id])
        x=cps[5] *M_sun_cgs/M_ear_cgs
        weights = np.ones(len(x))/Nsim #float(len(x)) #, density=True
        ax.hist(x, bins=5, range=(0,10), weights=weights, histtype="step",label=label, color=cmap[j], linewidth=lw)
    #ax.legend()
    k+=1    

    fig.savefig("ImpactParameters_psInc.png",dpi=300,bbox_inches='tight')
    
if __name__ == '__main__':
    main()
