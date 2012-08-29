from pylab import *
import numpy
import sys
sys.path.append('/home/aroach/MRI/pythonscripts')
import readzeus
import glob
import math

newrin=7.06
newrout=20.30

def avgr(minfile, maxfile, zindex, component, nobc=0):
    allfiles = glob.glob('*.h5')
    files=[]
    for file in allfiles:
        if (file >= minfile) and (file <= maxfile):
           files.append(file)
    numfiles=len(files)

    if(nobc == 0):
        data=readzeus.zmp(files[0])
    else:
        data=readzeus.zmpnobc(files[0])

    valuearray = numpy.zeros((numfiles,data[component].shape[0]))

    i=0
    for file in files:
        if (nobc == 0):
            data=readzeus.zmp(file)
        else:
            data=readzeus.zmpnobc(file)
            
        valuearray[i,]=data[component][:,zindex]
        i=i+1
        
    valuemean=numpy.mean(valuearray,0)
    valuestd=numpy.std(valuearray,0)

    valuereturn = {'r': data['r'], 'z': data['z'][zindex], 'mean': valuemean, 'std': valuestd}
    return valuereturn


def pltavgr(minfile, maxfile, zindex, component, nobc=0):
    data = avgr(minfile, maxfile, zindex, component, nobc=nobc)
    titlestring = '<' + component + '>'
    labelstring = 'z=' + str(data['z'])
    errorbar(data['r'], data['mean'], yerr=data['std'], label=labelstring)
    title(titlestring)
    xlabel("r [cm]")
    ylabel(component + " [cm/sec]")
    legend()
    show()


def plttimeseries(minfile, maxfile, zindex, rindex, component, nobc=0):
    allfiles = glob.glob('*.h5')
    files=[]
    for file in allfiles:
        if (file >= minfile) and (file <= maxfile):
           files.append(file)
    numfiles=len(files)

    if(nobc == 0):
        data=readzeus.zmp(files[0])
    else:
        data=readzeus.zmpnobc(files[0])

    valuearray = numpy.zeros(numfiles)
    time = numpy.zeros(numfiles)

    i=0
    for file in files:
        if(nobc == 0):
            data=readzeus.zmp(file)
        else:
            data=readzeus.zmpnobc(file)
            
        valuearray[i]=data[component][rindex,zindex]
        time[i]=data['time']
        i=i+1

    print(valuearray)
    print(time)
    titlestring = 'Time series of ' + component
    labelstring = 'z=' + str(data['z'][zindex]) + ', r=' + str(data['r'][rindex])
    plot(time, valuearray, label=labelstring)
    title(titlestring)
    xlabel("time [s]")
    ylabel(component + " [cm/sec]")
    legend()
    show()
                                                               
def save_omega_at_height(filename, zindex, savefilename):
    data = readzeus.zmp(filename)
    omega = data['vtheta'][:,zindex]/data['r']
    numpy.savetxt(savefilename, c_[data['r'], omega], fmt="%12.6G")

def pltatheight(filename, component, zindex, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = component + ' at t = ' + str(data['time']) + ' sec.'
    labelstring = 'z = ' + str(data['z'][zindex])
    plot(data['r'], data[component][:,zindex], label=labelstring)
    title(titlestring)
    xlabel("r [cm]")
    legend()
    show()


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


def pltcouette(minfile, maxfile, zindex, omegain, omegaout, nobc=0):
    avgprofile=avgr(minfile, maxfile, zindex, 'vtheta', nobc=nobc)
    couette = gencouette(avgprofile['r'], omegain, omegaout, newrin, newrout)

    clf()
    titlestring = '<V_theta>, ideal Couette profile at z = ' + str(avgprofile['z']) + 'cm'
    plot(couette['r'], couette['vtheta'], label='Ideal Couette')
    errorbar(avgprofile['r'], avgprofile['mean'], yerr=avgprofile['std'], label='Avg. Profile')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("V_theta [cm/sec]")
    legend()
    show()


def pltvthetacontour(filename, nobc=0):
    clf()
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)
        
    titlestring = 'Azimuthal velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vtheta'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")
    show()


def pltheightavgprofile(file, omegain, omegaout, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(file)
    else:
        data = readzeus.zmpnobc(file)
        
    avgprofile = numpy.mean(data['vtheta'],1)
    stdprofile = numpy.std(data['vtheta'],1)

    couette = gencouette(data['r'], omegain, omegaout, newrin, newrout)

    titlestring = 'Z-averaged V_theta'
    plot(couette['r'], couette['vtheta'], label='Ideal Couette')
    errorbar(data['r'], avgprofile, yerr=stdprofile, label='Avg. Profile')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("V_theta [cm/sec]")
    legend()


def pltvthetacontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)
        
    titlestring = 'Azimuthal velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vtheta'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltvrcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)


    titlestring = 'Radial velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vr'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltvzcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = 'Axial velocity [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['vz'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")



def pltbthetacontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = 'Azimuthal magnetic field at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['btheta'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltbrcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = 'Radial magnetic field at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['br'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltbzcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = 'Axial magnetic field [cm/sec] at t = ' + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], data['bz'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltlcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    titlestring = 'Angular momentum [cm^2/sec] at t = ' + str(data['time']) + 'sec.'
    contourf(data['r'], data['z'], data['l'].T, 75)
    colorbar()
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")


def pltstreamcontour(filename, nobc=0, levels=50):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    stream = generatestream(filename, nobc=nobc)

    titlestring = 'Stream function at t = ' + str(data['time']) + ' sec.'
    contour(data['r'], data['z'], stream.T, levels, colors='k')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")

def pltfluxcontour(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    flux = generateflux(filename, nobc=nobc)
    titlestring = 'Flux function at t = ' + str(data['time']) + ' sec.'
    contour(data['r'], data['z'], flux.T, 50, colors='k')
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")

def pltjr(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    jr = generatejr(filename, nobc=nobc)
    titlestring = "$j_r$ at t = " + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], jr.T, 50)
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")
    colorbar()


def pltjz(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    jz = generatejz(filename, nobc=nobc)
    titlestring = "$j_z$ at t = " + str(data['time']) + ' sec.'
    contourf(data['r'], data['z'], jz.T, 50)
    title(titlestring)
    xlabel("r [cm]")
    ylabel("z [cm]")
    colorbar()


def pltj2d(filename, nobc=0):
    if (nobc == 0):
        data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    jr = generatejr(filename, nobc=nobc)
    jz = generatejz(filename, nobc=nobc)
    quiver(data['r'], data['z'], jr.T, jz.T)

def pltjstream(filename, nobc=0, color=1):
    if (nobc == 0):
         data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    jr = generatejr(filename, nobc=nobc)
    jz = generatejz(filename, nobc=nobc)
    jstream = generatejstream(data['r'], data['z'], jr, jz)
    if (color==1):
        contourf(data['r'], data['z'], jstream.T, 50)

    contour(data['r'], data['z'], jstream.T, 50, colors='k')
    xlabel("r [cm]")
    ylabel("z [cm]")

def generatestream(filename, nobc=0):
    if (nobc == 0):
         data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    # Generates stream function array of correct size
    stream = data['vr']
    dz = data['z'][10]-data['z'][9]

    for i in range(stream.shape[0]):
        for j in range(stream.shape[1]):
            stream[i,j]=dz*data['r'][i]*data['vr'][i,j]
            if j != 0:
                stream[i,j]=stream[i,j]+stream[i,j-1]

    return stream


def generatejstream(r, z, jr, jz):
    # Generates stream function array of correct size
    stream = numpy.zeros(jr.shape)
    dz = z[10]-z[9]

    for i in range(stream.shape[0]):
        for j in range(stream.shape[1]):
            stream[i,j]=-dz*r[i]*jr[i,j]
            if j != 0:
                stream[i,j]=stream[i,j]+stream[i,j-1]

    return stream

def generateflux(filename, nobc=0):
    if (nobc == 0):
         data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    # Generates flux function array of correct size
    flux = data['br']
    dz = data['z'][10]-data['z'][9]

    for i in range(flux.shape[0]):
        for j in range(flux.shape[1]):
            flux[i,j]=dz*data['r'][i]*data['br'][i,j]
            if j != 0:
                flux[i,j]=flux[i,j]+flux[i,j-1]

    return flux

def generatejr(filename, nobc=0):
    if (nobc == 0):
         data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    jr = numpy.zeros(data['btheta'].shape)

    dz = data['z'][10]-data['z'][9]

    #Keep in mind that 'i' is the r coordinate, and 'j' is the z coordinate
    #We're going to calculate jr in amps/cm^2. B_theta is in gauss.
    #jr in statamps/cm^2 will require multiplying by c/4pi. The converstion
    #to amps will require dividing by 2.9979e9. So the total conversion
    #factor is 10/(4*pi)
    for i in range(jr.shape[0]):
        for j in range(1,jr.shape[1]-1):
            jr[i,j] = -(data['btheta'][i,j+1]-data['btheta'][i,j-1])/(2*dz)
        jr[i,0] = -(data['btheta'][i,1] - data['btheta'][i,0])/dz

    jr = jr*10.0/(4.0*numpy.pi)
    return jr



def generatejz(filename, nobc=0):
    if (nobc == 0):
         data = readzeus.zmp(filename)
    else:
        data = readzeus.zmpnobc(filename)

    btheta = data['btheta']
    jz = numpy.zeros(btheta.shape)
    

    r = data['r']

    #Keep in mind that 'i' is the r coordinate, and 'j' is the z coordinate
    #We're going to calculate jz in amps/cm^2. B_theta is in gauss.
    #jr in statamps/cm^2 will require multiplying by c/4pi. The converstion
    #to amps will require dividing by 2.9979e9. So the total conversion
    #factor is 10/(4*pi)
    #Jz = (1/r)(d/dr)(rB_theta)
    for i in range(1,jz.shape[0]-1):
        for j in range(jz.shape[1]):
            jz[i,j] = (1.0/r[i])*btheta[i,j] + (btheta[i+1,j] - btheta[i-1,j])/(2*(r[i+1]-r[i-1]))

    #for j in range(jz.shape[1]):
    #    jz[0,j] = (1/r[0])*btheta[0,j] + (btheta[1,j]-btheta[0,j])/(r[1]-r[0])


    jz = jz*10.0/(4.0*numpy.pi)
    return jz

          
