#!/usr/bin/python
#Compute global biomass from metabolic theory

import datetime as dt
import numpy as np
import netCDF4 as nc
import scipy.io as sio

def main():
    #initialize model run
    print ' '
    print 'Trophic efficiency model started at ', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    sst=14.0
    gpp=.15
    par=24.0
    
    ncat=4
    niter=1000
    
    orgmass=np.zeros([ncat,1],dtype='f8')
    endo=np.zeros([ncat,1],dtype='f8')
    tclass='H' * ncat
    
    # marine size structure
    orgslope=14.0/ncat
    for ii in range(0,ncat):
        orgmass[ii]=10**(orgslope*ii-10)
    
    # terrestrial size structure
    #orgslope=7.0/ncat
    #for ii in range(0,ncat):
    #    orgmass[ii]=10**(orgslope*ii-5)
    #orgmass[0]=1

    tte=np.zeros(ncat)

    #determine array attributes and preallocate memory for data results
    total_mass=np.zeros([niter],dtype='f8')
    mass_dens=np.zeros([niter,ncat],dtype='f8')
    flow_mat=np.zeros([niter,ncat,ncat],dtype='f8')
    tte_mat=np.zeros([niter,ncat],dtype='f8')
    
    #create output file
    outfile = '/Data/cobialab/conserveandclimate/model_data/biomass_cat4_' + dt.datetime.strftime(dt.datetime.now(),'%Y%m%d') + '.nc'

    biodata = create_biomassnc(outfile,niter,ncat)
    print 'Output datafile',outfile[45:], 'created at ', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    #perform model computations
    print 'Starting model computations at', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    tstepp=datetime2matlabdn(dt.datetime.now())
    
    for ii in range(0,niter):
        #generate flow matrix
        tflow=np.zeros([ncat,ncat],dtype='f8')
        #tte=np.ones(ncat)*.1
        tmax=np.random.randint(20,460,size=1)*.01
        tte=np.random.randint(50,tmax*1000,size=ncat)*.0001
        tte[0]=0.0
            
        for jj in range(0,(ncat-1)):
            tflow[jj,jj+1]=1
            if ii==0:
                tte[jj]=0.1

        if ii<(niter/4):
            tflow=tflow
        elif ii>=(niter/4) and ii<(2*niter/4):
            tflow[ncat-3,ncat-2]=0.6
            tflow[ncat-4,ncat-2]=0.3
            tflow[ncat-5,ncat-2]=0.1
        elif ii>=(2*niter/4) and ii<(3*niter/4):
            tflow[ncat-2,ncat-1]=0
            tflow[2,ncat-1]=1
        else:
            tflow[ncat-2,ncat-1]=0
            tflow[2,ncat-1]=1
            tflow[ncat-3,ncat-2]=0.6
            tflow[ncat-4,ncat-2]=0.3
            tflow[ncat-5,ncat-2]=0.1

        #distribute energy and account for trophic transfer
        tflow[0,0]=1
        for jj in range(0,ncat):
            for kk in range(1,ncat):
                if tflow[jj,kk]>0:
                    tflow[jj,kk]=tflow[jj,kk]*np.sum(tflow[:,jj])

        flow_mat[ii,:,:]=tflow
        tte_mat[ii,:]=tte

        mass_dens[ii,:]=compute_sat_biomass(sst,gpp,par,tflow,endo,orgmass,tclass)
        total_mass[ii]=np.sum(mass_dens[ii,:])

        if np.mod(ii+1,(niter/20))==0:
            tstep=datetime2matlabdn(dt.datetime.now())
            print 'Run is',(float((ii+1))/niter*1.0)*100.0,'% done in',round((tstep-tstepp)*86400,2),' seconds'
            tstepp=tstep

    print 'Model computations completed at', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    biodata.variables['org_mass'][:] = orgmass
    biodata.variables['tte'][:,:] = tte_mat
    biodata.variables['total_mass'][:] = total_mass
    biodata.variables['mass_dens'][:,:] = mass_dens
    biodata.variables['flow_matrix'][:,:,:] = flow_mat
    print 'Model data written at', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    print 'DONE...'

#subfunction compute_sat_biomass to compute ecosystem parameters
def compute_sat_biomass(sst1,gpp1,par1,tweb,endo1,orgmass1,tclass1):

    k=8.62e-5 #Boltzman constant
    #set constant variables for marine phytoplankton (Lopez-Urrutia et al 2006)
    Po_p=1.32e9; E_p=0.32; beta_p=1.036;
    #set constant variables for zooplankton
    Po_z=1.09e12; E_z=0.65; beta_z=0.86;
    #set constant variables for higher trophic levels
    Po=2.16e9; E=0.65; beta=0.75;
    
    #trim file to number of trophic compartments
    nn=orgmass1.size
    flow2=np.zeros([nn,nn],dtype='f8')
    T=np.zeros(nn,dtype='f8')
    Pind=np.zeros(nn,dtype='f8')
    gross_prod=np.zeros(nn,dtype='f8')
    net_prod=np.zeros(nn,dtype='f8')
    population=np.zeros(nn,dtype='f8')
    group_massdensity=np.zeros(nn,dtype='f8')
    trop_trans=np.zeros(nn,dtype='f8')
    
    #compute trophic transfer efficiency (TTE)
    #trop_trans=np.ones(nn,dtype='int64')*0.1
    #trop_trans[1:nn]=-.0112*np.log10(orgmass1[1:nn])+.0656; #log-linear varying from 20% to 2%
    #sens_test=trop_trans/.1 # ratio of trophic transfer to 10% TT
    nper=0.5
    
    #compute region specific flow matrix with new TTE
    for ii in range(0,nn):
        for jj in range(0,nn):
           flow2[ii,jj]=tweb[ii,jj]*gpp1

    #print flow2
    #compute phytoplankton size
    #Ps=12.19*np.log(gpp1)+37.248
    #orgmass1[0]=10**(.06*(100-Ps)-15)
    
    #compute fluxout proportion
    #fluxout=10**(.94*np.arctan(6*(np.log10(orgmass1)+1))-1.45);

    #compute ecosystem properties ind production, gross production, net
    #production, population size, and group total biomass for each trophic
    #compartment based on Schramski et al 2015 (PNAS)
    
    #print orgmass[5],endo[5],tclass[5],flow[5,:]
    
    for ii in range(0,nn):
        if endo1[ii]==0:
           T[ii]=273+sst1
        elif endo1[ii]==1:
           T[ii]=273+37
        if tclass1[ii] == 'P':
            Pind[ii]=Po_p*(orgmass1[ii]**beta_p)*np.exp(-E_p/(k*T[ii]))*par1
        elif tclass1[ii] == 'Z':
            Pind[ii]=Po_z*(orgmass1[ii]**beta_z)*np.exp(-E_z/(k*T[ii]))
        elif tclass1[ii] == 'H':
            Pind[ii]=Po*(orgmass1[ii]**beta)*np.exp(-E/(k*T[ii]))
        gross_prod[ii]=np.sum(flow2[:,ii])
        net_prod[ii]=nper*gross_prod[ii]
        population[ii]=net_prod[ii]/Pind[ii]
        group_massdensity[ii]=population[ii]*orgmass1[ii]

    return group_massdensity

#subfunction create_biomassnc to write out output netcdf file
def create_biomassnc(outfile,n,m):
    biodata = nc.Dataset(outfile, 'w')
    biodata.createDimension('niter', n)
    biodata.createDimension('comp', m)
    biodata.createDimension('ones', 1)

    #Total mass density
    biodata.createVariable('org_mass', 'f8', ('comp'))
    biodata.variables['org_mass'].long_name = 'organism mass'
    biodata.variables['org_mass'].units = 'kg C ind^-1'

    #Trophic transfer efficiency
    biodata.createVariable('tte', 'f8', ('niter', 'comp'))
    biodata.variables['tte'].long_name = 'trophic transfer efficiency'
    biodata.variables['tte'].units = 'ratio'

    #Total mass density
    biodata.createVariable('total_mass', 'f8', ('niter'))
    biodata.variables['total_mass'].long_name = 'total biomass density'
    biodata.variables['total_mass'].units = 'kg C m^-2'
    
    #biomass density
    biodata.createVariable('mass_dens', 'f8', ('niter', 'comp'))
    biodata.variables['mass_dens'].long_name = 'biomass density by trophic category'
    biodata.variables['mass_dens'].units = 'kg C m^-2'

    #energy based trophic web
    biodata.createVariable('flow_matrix', 'f8', ('niter', 'comp', 'comp'))
    biodata.variables['flow_matrix'].long_name = 'energy based trophic web'
    biodata.variables['flow_matrix'].units = 'proportion of energy'

    return biodata

def datetime2matlabdn(inp):
    ord = inp.toordinal()
    mdn = inp+dt.timedelta(days=366)
    frac = (inp-dt.datetime(inp.year,inp.month,inp.day,0,0,0)).seconds/(24.0*60.0*60.0)
    return mdn.toordinal() + frac

def perms(n):
    if not n:
        return
    
    for i in xrange(2**n):
        s = bin(i)[2:]
        s = "0" * (n-len(s)) + s
        yield s

main()
















