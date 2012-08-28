import sys
import getopt

import numpy
import scipy
import tables
import numpy


#z2d is used to read hydro datafiles from ZEUS-2D
def z2d(filename):
    h5file = tables.openFile(filename, 'r')
    datasetr = h5file.getNode('/', 'fakeDim0')
    datasetz = h5file.getNode('/', 'fakeDim1')
    datasetd = h5file.getNode('/', 'Data-Set-2')
    datasete = h5file.getNode('/', 'Data-Set-3')
    datasetvz = h5file.getNode('/', 'Data-Set-4')
    datasetvr = h5file.getNode('/', 'Data-Set-5')
    datasetvtheta = h5file.getNode('/', 'Data-Set-6')

    # Get the time from the attribute string
    attributestring = datasetd._f_getAttr('long_name')
    timestring = attributestring.split("=")[1]
    time = float(timestring)
    
    #Read the data sets
    r = datasetr.read()
    z = datasetz.read()
    d = datasetd.read()
    e = datasete.read()
    vztemp = datasetvz.read()
    vrtemp = datasetvr.read()
    vtheta = datasetvtheta.read()
    h5file.close()

    # Now we have to account for offsets.  There are three ghost zones
    # included in the data, which we trim out.  Also, the coordinates are
    # for *cell-centered* quantities.  The vr and vz components are
    # edge-centered.  So we must average the components on both sides of
    # the cell to yield the correct quantities
    
    r = r[1:-2]
    z = z[1:-2]
    vtheta = vtheta[1:-2,1:-2]
    d = d[1:-2,1:-2]
    e = e[1:-2,1:-2]

    z = z - 5.0

    # vr and vz are initialized from vrtemp and vztemp just to get the size
    # right.  The data is replaced in the for loop below.
    vr = vrtemp[1:-2,1:-2]
    vz = vztemp[1:-2,1:-2]

    for i in range(vrtemp.shape[0]-3):
        for j in range(vrtemp.shape[1]-3):
          vr[i,j]=0.5*(vrtemp[i+1,j+1]+vrtemp[i+2,j+1])
          vz[i,j]=0.5*(vztemp[i+1,j+1]+vztemp[i+1,j+1])


    #Calculate angular momentum
<<<<<<< HEAD
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(l.shape[0]):
        for j in range(l.shape[1]):
            l[i,j]=vtheta[i,j]*r[i]

    #Calculate angular velocity
<<<<<<< HEAD
    rinv = numpy.zeros(r.shape[0],dtype=numpy.float32)
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    rinv = numpy.zeros(r.shape[0])
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(omega.shape[0]):
        for j in range(omega.shape[1]):
            omega[i,j]=vtheta[i,j]*rinv[i]

    #Return a dictionary with all of the data sets
    valuearray = {'r': r, 'z': z, 'd': d, 'e': e, 'vz': vz, 'vr': vr, 'vtheta': vtheta, 'time': time, 'l': l, 'omega': omega}
    return valuearray
    

#zmp is used to read Zeus-MP MHD datafiles
def zmp(filename):
    h5file = tables.openFile(filename, 'r')
    datasetr = h5file.getNode('/', 'fakeDim1')
    datasetz = h5file.getNode('/', 'fakeDim2')
    datasetvz = h5file.getNode('/', 'Data-Set-2')
    datasetvr = h5file.getNode('/', 'Data-Set-3')
    datasetvtheta = h5file.getNode('/', 'Data-Set-4')
    datasetbz = h5file.getNode('/', 'Data-Set-5')
    datasetbr = h5file.getNode('/', 'Data-Set-6')
    datasetbtheta = h5file.getNode('/', 'Data-Set-7')
    datasetd = h5file.getNode('/', 'Data-Set-8')
    datasete = h5file.getNode('/', 'Data-Set-9')

    # Get the time from the attribute string
    attributestring = datasetd._f_getAttr('long_name')
    timestring = attributestring.split("=")[1]
    time = float(timestring)
    
    #Read the data sets
    r = datasetr.read()
    z = datasetz.read()
    d = datasetd.read()
    e = datasete.read()
    vztemp = datasetvz.read()
    vrtemp = datasetvr.read()
    bztemp = datasetbz.read()
    brtemp = datasetbr.read()
    vtheta = datasetvtheta.read()
    btheta = datasetbtheta.read()
    h5file.close()

    # Now we have to account for offsets.  There are three ghost zones
    # included in the data, which we trim out.  Also, the coordinates are
    # for *cell-centered* quantities.  The vr, vz, br, and bz components are
    # edge-centered.  So we must average the components on both sides of
    # the cell to yield the correct quantities
    
    r = r[1:-2]
    z = z[1:-2]
    vtheta = vtheta[0,1:-2,1:-2]
    btheta = btheta[0,1:-2,1:-2]
    d = d[0,1:-2,1:-2]
    e = e[0,1:-2,1:-2]

    # Accounts for artificial 5cm offset in z coordinates 
    z = z - 5.0

    # vr, vz, br, and bz are initialized from vrtemp and vztemp just to
    # get the size right.  The data is replaced in the for loop below.
    vr = vrtemp[0,1:-2,1:-2]
    vz = vztemp[0,1:-2,1:-2]
    br = brtemp[0,1:-2,1:-2]
    bz = bztemp[0,1:-2,1:-2]

    #Now do the averaging
    for i in range(vr.shape[0]):
        for j in range(vr.shape[1]):
          vr[i,j]=0.5*(vrtemp[0,i+1,j+1]+vrtemp[0,i+2,j+1])
          vz[i,j]=0.5*(vztemp[0,i+1,j+1]+vztemp[0,i+1,j+1])
          br[i,j]=0.5*(brtemp[0,i+1,j+1]+brtemp[0,i+2,j+1])
          bz[i,j]=0.5*(bztemp[0,i+1,j+1]+bztemp[0,i+1,j+1])


    #Calculate angular momentum
<<<<<<< HEAD
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(l.shape[0]):
        for j in range(l.shape[1]):
            l[i,j]=vtheta[i,j]*r[i]

    #Calculate angular velocity
<<<<<<< HEAD
    rinv = numpy.zeros(r.shape[0],dtype=numpy.float32)
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    rinv = numpy.zeros(r.shape[0])
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(omega.shape[0]):
        for j in range(omega.shape[1]):
            omega[i,j]=vtheta[i,j]*rinv[i]

    #Return a dictionary with all of the data sets
    valuearray = {'r': r, 'z': z, 'd': d, 'e': e, 'vz': vz, 'vr': vr, 'vtheta': vtheta, 'bz': bz, 'br': br, 'btheta': btheta, 'time': time, 'l': l, 'omega': omega}
    return valuearray



#zmpnobc is used to read Zeus-MP MHD datafiles *without* the partially-conducting
#boundary condition segments
def zmpnobc(filename):
    h5file = tables.openFile(filename, 'r')
    datasetr = h5file.getNode('/', 'fakeDim1')
    datasetz = h5file.getNode('/', 'fakeDim2')
    datasetvz = h5file.getNode('/', 'Data-Set-2')
    datasetvr = h5file.getNode('/', 'Data-Set-3')
    datasetvtheta = h5file.getNode('/', 'Data-Set-4')
    datasetbz = h5file.getNode('/', 'Data-Set-5')
    datasetbr = h5file.getNode('/', 'Data-Set-6')
    datasetbtheta = h5file.getNode('/', 'Data-Set-7')
    datasetd = h5file.getNode('/', 'Data-Set-8')
    datasete = h5file.getNode('/', 'Data-Set-9')

    # Get the time from the attribute string
    attributestring = datasetd._f_getAttr('long_name')
    timestring = attributestring.split("=")[1]
    time = float(timestring)
    
    #Read the data sets
    r = datasetr.read()
    z = datasetz.read()
    d = datasetd.read()
    e = datasete.read()
    vztemp = datasetvz.read()
    vrtemp = datasetvr.read()
    bztemp = datasetbz.read()
    brtemp = datasetbr.read()
    vtheta = datasetvtheta.read()
    btheta = datasetbtheta.read()
    h5file.close()

    # Now we have to account for offsets.  There are three ghost zones
    # included in the data, which we trim out.  Also, the coordinates are
    # for *cell-centered* quantities.  The vr, vz, br, and bz components are
    # edge-centered.  So we must average the components on both sides of
    # the cell to yield the correct quantities

    #We also trim out additional zonesegments from the
    #radial coordinate.  These are the partially conducting inner and outer
    #rings.
    
    r = r[8:-10]
    z = z[1:-2]
    vtheta = vtheta[0,8:-10,1:-2]
    btheta = btheta[0,8:-10,1:-2]
    d = d[0,8:-10,1:-2]
    e = e[0,8:-10,1:-2]

    # Accounts for artificial 5cm offset in z coordinates 
    z = z - 5.0

    # vr, vz, br, and bz are initialized from vrtemp and vztemp just to
    # get the size right.  The data is replaced in the for loop below.
    vr = vrtemp[0,8:-10,1:-2]
    vz = vztemp[0,8:-10,1:-2]
    br = brtemp[0,8:-10,1:-2]
    bz = bztemp[0,8:-10,1:-2]

    #Now do the averaging
    for i in range(vr.shape[0]):
        for j in range(vr.shape[1]):
          vr[i,j]=0.5*(vrtemp[0,i+8,j+1]+vrtemp[0,i+9,j+1])
          vz[i,j]=0.5*(vztemp[0,i+8,j+1]+vztemp[0,i+8,j+2])
          br[i,j]=0.5*(brtemp[0,i+8,j+1]+brtemp[0,i+9,j+1])
          bz[i,j]=0.5*(bztemp[0,i+8,j+1]+bztemp[0,i+8,j+2])

    #Calculate angular momentum - copy vr first to get shape right
<<<<<<< HEAD
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    l = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(l.shape[0]):
        for j in range(l.shape[1]):
            l[i,j]=vtheta[i,j]*r[i]

    #Calculate angular velocity
<<<<<<< HEAD
    rinv = numpy.zeros(r.shape[0],dtype=numpy.float32)
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]),dtype=numpy.float32)
=======
    rinv = numpy.zeros(r.shape[0])
    for i in range(rinv.shape[0]):
        rinv[i]=1.0/r[i]
    
    omega = numpy.zeros((vtheta.shape[0],vtheta.shape[1]))
>>>>>>> 9d1c6097b90e114e95f888c9cce0f2729d8ecd6a
    for i in range(omega.shape[0]):
        for j in range(omega.shape[1]):
            omega[i,j]=vtheta[i,j]*rinv[i]
    
    #Return a dictionary with all of the data sets
    valuearray = {'r': r, 'z': z, 'd': d, 'e': e, 'vz': vz, 'vr': vr, 'vtheta': vtheta, 'bz': bz, 'br': br, 'btheta': btheta, 'time': time, 'l': l, 'omega': omega}
    return valuearray
