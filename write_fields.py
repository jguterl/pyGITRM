"""
Created on Wed Nov 15 17:39:55 2023

@author: nathd
"""
import numpy as np
import netCDF4
import scipy.linalg
import matplotlib
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

CREATE_NETCDF_FILE=0
CREATE_LU_FILE=0
'''
Read Magnetic field File
file_nc='/lore/nathd/DIII-D-GA/DIMES/Fields/bfield_176509.nc'
ds1_all = netCDF4.Dataset(file_nc)
edist1 = ds1_all['r'][:]
'''

#Fields except the mesh
data=np.load('/lore/nathd/DIII-D-GA/DIMES/Fields/test_167196_bgp_mesh.npy', allow_pickle=True).tolist()
key_list=[k for k in data.keys() if k!='mesh']

#Mesh
mesh_array=np.array(data['mesh'])

#Close any previous figure
plt.close('all')

'''
These are the keywords and associated fields and associated names
['KNBS', 'KTEBS', 'KTIBS', 'KVHS', 'KES', 'TEGS', 'TIGS']
'''
field_name=["ne", "te", "ti", "Velocity", "eField","gradTe","gradTi"];
ind=0

for k in key_list:
    if data.get(k) is None:
        data[k] = 0*data['KNBS']
        print("This field has no data:",k)
    
    field_array=np.array(data[k])
    #Figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    patches = [matplotlib.patches.Polygon(v, closed=False) for v in mesh_array]
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1, edgecolor=None)
    p.set_array(field_array)
    ax.add_collection(p)
    ax.set_xlim([1, 2.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_aspect('equal', adjustable='box')
    plt.xticks(np.arange(1, 2.6, step=0.5), fontsize = 20) 
    plt.yticks(fontsize = 20) 
    plt.xlabel("r (m)", fontsize=20)
    plt.ylabel("z (m)", fontsize=20)
    plt.title(field_name[ind],fontsize=25)
    cbar=fig.colorbar(p)
    cbar.ax.tick_params(axis="both", labelsize=20) 
    ax.set(xticks=[], yticks=[])
    cbar.ax.yaxis.get_offset_text().set(size=20)

    ind=ind+1

if CREATE_NETCDF_FILE:
    ne=np.array(data['KNBS'])
    te=np.array(data['KTEBS'])
    ti=np.array(data['KTIBS'])
    v_t=np.array(data['KVHS'])
    e_f=np.array(data['KES'])
    gradTe = 0*e_f
    gradTi = 0*e_f

    # Create NETCDF file for fields
    profiles_filename = "profiles_OEDGE_dimes.nc"
    rootgrp = netCDF4.Dataset(profiles_filename, "w", format="NETCDF4")

    #Dimesnion
    nx=np.shape(ne)[0]
    ny=1
    nrr = rootgrp.createDimension("nR", nx)
    nzz = rootgrp.createDimension("nZ", ny)

    #Electron Density
    nee = rootgrp.createVariable("ne", "f8", ("nZ", "nR"))
    nee[:]=ne

    #Ion Density
    nii = rootgrp.createVariable("ni", "f8", ("nZ", "nR"))
    nii[:]=ne

    #Electron Temperature
    tee = rootgrp.createVariable("te", "f8", ("nZ", "nR"))
    tee[:]=te

    #Ion Temperature
    tii = rootgrp.createVariable("ti", "f8", ("nZ", "nR"))
    tii[:]=ti

    # Velocity
    vrr  = rootgrp.createVariable("vr", "f8", ("nZ", "nR"))
    vrr[:]=0*ti

    vtt  = rootgrp.createVariable("vt", "f8", ("nZ", "nR"))
    vtt[:]=0*ti

    vzz  = rootgrp.createVariable("vz", "f8", ("nZ", "nR"))
    vzz[:]=0*ti
    
    #Electric Field
    Err = rootgrp.createVariable("Er", "f8", ("nZ", "nR"))
    Err[:] =0*ti

    Ett = rootgrp.createVariable("Et", "f8", ("nZ", "nR"))
    Ett[:] =0*ti

    Ezz = rootgrp.createVariable("Ez", "f8", ("nZ", "nR"))
    Ezz[:] = 0*ti

    #GradtiR
    te_grad_rr = rootgrp.createVariable("gradTeR", "f8", ("nZ", "nR"))
    te_grad_rr[:] = 0*ti
    te_grad_tt = rootgrp.createVariable("gradTeT", "f8", ("nZ", "nR"))
    te_grad_tt[:] = 0*ti
    te_grad_zz = rootgrp.createVariable("gradTeZ", "f8", ("nZ", "nR"))
    te_grad_zz[:] = 0*ti

    #GradteR
    ti_grad_rr = rootgrp.createVariable("gradTiR", "f8", ("nZ", "nR"))
    ti_grad_rr[:] = 0*ti
    ti_grad_tt = rootgrp.createVariable("gradTiT", "f8", ("nZ", "nR"))
    ti_grad_tt[:] = 0*ti
    ti_grad_zz = rootgrp.createVariable("gradTiZ", "f8", ("nZ", "nR"))
    ti_grad_zz[:] = 0*ti

    #Flux
    fnax00= rootgrp.createVariable("fnax", "f8", ("nZ", "nR"))
    fnax00[:]=0*ti+1
    
    rootgrp.close()


if CREATE_LU_FILE:
    mesh_array=np.array(data['mesh'])
    mesh_array=mesh_array

    tri=np.zeros((3,3,2*ny*nx)) 
    for i in range(nx*ny): 
        tri[:,:,i*2] =  [[mesh_array[i,0,0],mesh_array[i,1,0],mesh_array[i,2,0]],\
                          [mesh_array[i,0,1],mesh_array[i,1,1],mesh_array[i,2,1]],[1,1,1]] 
            
        tri[:,:,i*2+1] = [[mesh_array[i,0,0],mesh_array[i,2,0],mesh_array[i,3,0]],\
                          [mesh_array[i,0,1],mesh_array[i,2,1],mesh_array[i,3,1]],[1,1,1]]
            
    #First index loops over all elements, second over the 4 vertices and 3rd over x/y
    #For OEDGE data the sequence of vertices of the quad is 0, 1, 2, 3 (SOLPS is different) 

    l=np.zeros((3,3,2*ny*nx)) 
    p=np.zeros((3,3,2*ny*nx)) 
    u=np.zeros((3,3,2*ny*nx)) 
    for i in range(2*ny*nx): 
        P, L, U = scipy.linalg.lu(tri[:,:,i]) 
        l[:,:,i]=L[:,:]; 
        p[:,:,i]=np.transpose(P[:,:]); 
        u[:,:,i]=U[:,:]; 

    mesh_filename = "lu.nc"
    rootgrp = netCDF4.Dataset(mesh_filename, "w", format="NETCDF4")

    h  = rootgrp.createDimension("h", 3)
    w  = rootgrp.createDimension("w", 3)
    m_s= rootgrp.createDimension("m_s", 2*nx*ny)

    pf = rootgrp.createVariable("P", "f8", ("h","w","m_s"))
    pf[:]=p
    lf = rootgrp.createVariable("L", "f8", ("h","w","m_s"))
    lf[:]=l
    uf = rootgrp.createVariable("U", "f8", ("h","w","m_s"))
    uf[:]=u
    rootgrp.close()
    
    print(p[:,:,0])
    print(l[:,:,0])
    print(u[:,:,0])
