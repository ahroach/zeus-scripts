from pylab import *
import numpy
import sys
sys.path.append('/home/aroach/MRI/pythonscripts')
import readzeus
import glob
import math
import scipy.io.array_import
import time
from scipy import optimize

protorin=3.8
protorout=14.9
ion()


def avgr(minfile, maxfile, zindex, component):
    allfiles = glob.glob('*.h5')
    files=[]
    for file in allfiles:
        if (file >= minfile) and (file <= maxfile):
           files.append(file)
    numfiles=len(files)

    data=readzeus.z2d(files[0])
    valuearray = numpy.zeros((numfiles,data[component].shape[0]),
                             dtype=float32)

    i=0
    for file in files:
        data=readzeus.z2d(file)
        valuearray[i,]=data[component][:,zindex]
        print "Now reading " + file
        i=i+1
        
    valuemean=numpy.mean(valuearray,0)
    valuestd=numpy.std(valuearray,0)

    valuereturn = {'r': data['r'], 'z': data['z'][zindex], 'mean': valuemean, 'std': valuestd}
    return valuereturn



def pltavgr(minfile, maxfile, zindex, component):
    data = avgr(minfile, maxfile, zindex, component)
    titlestring = '<' + component + '>'
    labelstring = 'z=' + str(data['z'])
    errorbar(data['r'], data['mean'], yerr=data['std'], label=labelstring)
    title(titlestring)
    xlabel("r [cm]")
    legend()


def pltatheight(filename, component, zindex):
    data = readzeus.z2d(filename)
    titlestring = component + ' at t = ' + str(data['time']) + ' sec.'
    labelstring = 'z = ' + str(data['z'][zindex])
    plot(data['r'], data[component][:,zindex], label=labelstring)
    title(titlestring)
    xlabel("r [cm]")
    legend()


def gencouette(r, omegain, omegaout, rin, rout):
    #Make vtheta with the correct dimensions
    vtheta=numpy.zeros((r.shape))

    vin = 2*math.pi*omegain*rin/60.0
    vout = 2*math.pi*omegaout*rout/60.0

    a = (vin*rin - vout*rout)/(rin**2-rout**2)
    b = (vin*rin - a*rin**2)

    for i in range(vtheta.shape[0]):
        vtheta[i]=(a*r[i]) + (b/r[i])

    valuereturn={'r': r, 'vtheta': vtheta}
    return valuereturn



def pltcouette(minfile, maxfile, zindex, omegain, omegaout):
    avgprofile=avgr(minfile, maxfile, zindex, 'vtheta')
    couette = gencouette(avgprofile['r'], omegain, omegaout, protorin, protorout)

    titlestring = '<V_theta>, ideal Couette profile at z = ' + str(avgprofile['z']) + 'cm'
    plot(couette['r'], couette['vtheta'], label='Ideal Couette')
    errorbar(avgprofile['r'], avgprofile['mean'], yerr=avgprofile['std'], label='Avg. Profile')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("V_theta [cm/sec]")
    legend()


def pltvthetacontour(filename):
    clf()
    data = readzeus.z2d(filename)
    titlestring = 'Azimuthal velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vtheta'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltheightavgprofile(file, omegain, omegaout):
    data = readzeus.z2d(file)
    avgprofile = numpy.mean(data['vtheta'],1)
    stdprofile = numpy.std(data['vtheta'],1)

    couette = gencouette(data['r'], omegain, omegaout, protorin, protorout)

    titlestring = 'Z-averaged V_theta'
    plot(couette['r'], couette['vtheta'], label='Ideal Couette')
    errorbar(data['r'], avgprofile, yerr=stdprofile, label='Avg. Profile')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("V_theta [cm/sec]")
    legend()


def pltvthetacontour(filename):
    clf()
    data = readzeus.z2d(filename)
    titlestring = 'Azimuthal velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vtheta'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltvrcontour(filename):
    clf()
    data = readzeus.z2d(filename)
    titlestring = 'Radial velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vr'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltvzcontour(filename):
    clf()
    data = readzeus.z2d(filename)
    titlestring = 'Axial velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vz'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltstreamcontour(filename):
    clf()
    data = readzeus.z2d(filename)
    stream = generatestream(filename)
    titlestring = 'Stream function at t = ' + str(data['time']) + ' sec.'
    contour(data['r'], data['z'], stream.T, 50, colors='k')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def generatestream(filename):
    data = readzeus.z2d(filename)

    # Generates stream function array of correct size
    stream = data['vr']
    dz = data['z'][10]-data['z'][9]

    for i in range(stream.shape[0]):
        for j in range(stream.shape[1]):
            stream[i,j]=-dz*data['r'][i]*data['vr'][i,j]
            if j != 0:
                stream[i,j]=stream[i,j]+stream[i,j-1]

    return stream


def pltmeasurements(shot_number, calibrationshot_number, calibrationrpm):
    calparams = gencalibrate(calibrationshot_number, calibrationrpm)
    filename = "/home/aroach/MRI/LDV_Data/" + shot_number + "/profile.csv"
    data = scipy.io.array_import.read_array(filename,separator=',',comment="\"")
    x = numpy.zeros(data[:,0].shape[0],dtype=float32)

    for i in range(x.shape[0]):
        x[i] = 0.1*(-data[:,0][i])+calparams['offset']

    y = 100*data[:,3]*calparams['mratio'] #Put data in cm/sec instead of m/sec
    flucy = y*data[:,4]

    errorbar(x, y, flucy, label="LDV Measurements")

def gencalibrate(shot_number, rpm):
    filename = "/home/aroach/MRI/LDV_Data/" + shot_number + "/profile.csv"
    data = scipy.io.array_import.read_array(filename,separator=',',comment="\"")
    r = .1*data[:,0]
    vtheta = 100*data[:,3]
    vthetaerr = data[:,3]*data[:,4]

    for i in range(vtheta.shape[0]-1,-1,-1):
        if(isnan(vtheta[i]) or (vtheta[i] == 0.0)):
            vtheta = delete(vtheta, i)
            r = delete(r, i)
            vthetaerr = delete(vthetaerr, i)

    fitfunc = lambda p, x: p[0] + p[1]*x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    pinit = [140, -5]

    out = optimize.leastsq(errfunc, pinit, args=(r, vtheta, vthetaerr),
                           full_output=1)

    pfinal = out[0]
    m = pfinal[1]
    b = pfinal[0]

    print "m = " + str(m)
    print "b = " + str(b)

    offset = -b/m
    expectedm = 2*3.14159*rpm/60
    mratio = -m/expectedm

    print "Giving values of:"
    print "Offset = " + str(offset)
    print "Mratio = " + str(mratio) + "    (Should be near 1.0)"

    fit = numpy.zeros(r.shape[0],dtype=float32)

    for i in range(r.shape[0]):
        fit[i] = m*r[i]+b

#    clf()
#    plot(r, fit)
#    errorbar(r, vtheta, vthetaerr)
#    show()
#    time.sleep(2)
#    clf()

    valuearray = {'offset': offset, 'mratio': mratio}
    return valuearray
