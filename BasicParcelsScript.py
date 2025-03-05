#!/usr/bin/env python
# coding: utf-8

# # Base script for Ocean Parcels with land avoidance and particle age

# Author: Johanna Wren
# 
# email: johanna.wren@noaa.gov
# 
# Date: June 11, 2021

# This script advects particles in pre-defined locations at a pre-defined depth using HYCOM currents stored on file. The trajectories are saved as a net cdf file. Since particles bunch up on the land border and then diffuse onto land, I added the land avoidance modification from These are the kernels from Erik's [*Delandmeter, P., and van Sebille, E. (2019). The Parcels v2.0 Lagrangian framework: new field interpolation schemes. Geoscientific Model Development 12, 3571–3584*](https://gmd.copernicus.org/articles/12/3571/2019/#section6), copied straight off. The original code can be found [here.](https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py)

# ## Load libraries
# We are not loading some of the basic libraries (e.g. `AdvectionRK4` and `DiffusionUniformKh`) here since we are modifying them for the land avoidance part. Rather, we are placing them in the code as kernels that we define instead. 

# In[1]:


# Load libraries etc.
from parcels import FieldSet, ParticleSet, JITParticle, ErrorCode, plotTrajectoriesFile, Variable, ParcelsRandom, Field, AdvectionRK4, DiffusionUniformKh
import numpy as np
import pandas as pd
import xarray as xr
from datetime import timedelta, datetime
import math

# ## Define custom kenels
# This is where we define other useful custom kenels.

# ### DELETE PARTICLES
# `DeleteParticles` allows for particles to be deleted when they exit the dispersal domain. This is applied to the borders of the domain but NOT to land. 

def DeleteParticle(particle, fieldset, time):
    particle.delete()


# ### AGEING
# The `Ageing` kernel removes particles after a specific time period (PLD) by calling the `DeleteParticles` kernel. This allows you to run simulation for long periods of time while removing particles after a set number of days.

def Ageing(particle, fieldset, time):
    particle.age += particle.dt
    if particle.age >= fieldset.pld:
        particle.delete()


# ## Create a fieldset
# We are creating a fieldset from HYCOM current files downloaded and stored on a drive. The fieldset is what OceanParcels uses to store all the information it needs for the flowfield that it moves particles through.
# 
# The HYCOM data was downloaded using the `HYCOM_download.jnl` script and grabs u and v currents for a region 120E-120W and 0-50N. The surface, 0-25m average, 0-50m average, and 0-100m average currents are calculated and saved. The files are saved as HYCOM_LAPS_????.nc files where the number indexing is days since the start of the dataseries on Jan 1, 1994. The HYCOM data can be found on the [APDRC website](http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_3151_b22a_aae6). 
# 
# Because these files were saved using Ferret the lon and lat variable names are a bit different. You can check the variable names by reading them into this script (see commented out section below) or through `ncdump -h HYCOM_LAPS_2648.nc` command in your terminal (on a unix/linux machine). 

filenames = {'U': "HYCOM_onaga/HYCOM_LAPS_*.nc",
             'V': "HYCOM_onaga/HYCOM_LAPS_*.nc"}
variables = {'U': 'AVG_U_25',
             'V': 'AVG_V_25'}
dimensions = {'lat': 'LATITUDE1001_1626',
              'lon': 'LONGITUDE3751_5251',
              'depth': 'LEV1_1',
              'time': 'TIME'}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)

# ## Add diffusivity
# In order for diffusivity to be incorporated we need to add zonal and meridional Kh fields. This code uses spatially homogenous diffusivity, but if you want spatially varying diffusivity use the `AdvectionDiffusionEM` or `AdvectionDiffusionM1` kernel instead.
# 
# Skip the next cell if you want to run without diffusivity. 

kh = 100   # This is the eddy diffusivity in m2/s

# Add even diffusivity to the fieldset
fieldset.add_constant_field('Kh_zonal', kh, mesh='spherical')
fieldset.add_constant_field('Kh_meridional', kh, mesh='spherical')            


# ## Add PLD to the fieldset
# I'm adding PLD as a constant to the fieldset here. YOU NEED TO DEFINE PLD HERE!!

# Add PLD to fieldset
pld = 365   # in days
fieldset.add_constant('pld', (pld*86400))


# ## Define your particle type
# We are using a JIT particle but we are adding the age variable along with beached and unbeachCount variables to them. These will print to the output file, but it's easy to change to not show in the ouput file. 

class LandAvoidParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)
    #beached = Variable('beached', dtype=np.int32, initial=0.)
    #unbeachCount = Variable('unbeachCount', dtype=np.int32, initial=0.)
    releaseSite = Variable('releaseSite', dtype=np.int32)


# ## Set particle release locations
# Next, we instantiate a `ParticeSet` composed of `JITParticle`. 
# 
# `ParticleSet` holds release locations, the number of particles released, the release depth, and how often a particel should be released (daily, every 6h, etc.)
# 
# Particles are released from locations specified in a file. We use this if we want to release at for example all reef locations in the islands, or any other list of locations. 
# 
# You can add or omit `pdepth`, `time` and `repeatdt` for a one time release at the default `fieldset` depth (in a multi depth file default is top layer) and time. 

# In[20]:


# Set input file
infile = pd.read_csv('LAPS_release_sites.csv')

# Set number of particles you want to release
npart = 10

# Set the second column as longitude and third column as latitude (remember that python indexing starts at 0!)
habilon=np.repeat(infile.Lon, npart)
habilat=np.repeat(infile.Lat, npart)
habisite = np.repeat(infile.Site,npart)

# Time interval between particle release (in seconds)
release_int = 86400*4

# Start date for release (if you want it different from the first day of the currents in the fielset)
#start_date = datetime(2000, 1, 16)

# Define the pset
pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=LandAvoidParticle, 
                             lon=habilon,
                             lat=habilat,
                             releaseSite=habisite,
                             repeatdt=release_int)


# ## Advect particles

# To invoke the land avoidance and custom kernels we wrote earlier we need to combine them into a format Parcels can use. Simply string all the kernels together that you want to use.


# Combine kernels you are using
#kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(BeachTesting_2D) + pset.Kernel(UnBeaching) + pset.Kernel(DiffusionUniformKh) + pset.Kernel(BeachTesting_2D) + pset.Kernel(Ageing)
kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(DiffusionUniformKh) + pset.Kernel(Ageing)


# Execute the advection. Depending on how many particles you release and the length of your release this may take some time. 
# 
# Particles specified in `ParticleSet` are advected using `kernels` at a `model_dt` time step and printed in a netcdf file (`outfile`) at `save_dt` time steps. Particels reaching the limits of the domain are removed through `recovery`. 

# Set output file name
outfile = "LAPS_n10_pld365_25m_nday730_01012000.nc"
# Time step in model
model_dt = timedelta(minutes=15)
# Time step to save to file
save_dt = timedelta(hours=24)
# Length of model run
run_days = 850

# Execute
pset.execute(kernels,
            runtime=timedelta(days=run_days),
            dt=model_dt, 
            recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
            output_file=pset.ParticleFile(name=outfile, outputdt=save_dt))


