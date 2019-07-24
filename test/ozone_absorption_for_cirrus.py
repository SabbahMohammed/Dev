import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal import savgol_filter
import sys
from scipy.interpolate import interp1d
import os


nml = 0.124
npts = 2**17
fname = "{}nm_npts{}".format(nml,npts)


if os.path.exists(fname):
    raise Exception("the same parameters found in {}".format(fname))
else:
    os.makedirs(fname)



# this function is being used to approx. model infrared basorption spectrum
def gauss(λ,a,b,c):
    return a*np.exp(-(λ-b)**2/(2*c**2))


# load data
data1 = np.loadtxt("ozone 52-131nm.txt")
data2 = np.loadtxt("ozone 110-337.4nm.txt")
data3 = np.loadtxt("ozone 176-845nm.txt")
data4 = np.loadtxt("ozone 213-1100nm.txt")
xray = np.loadtxt("ozone_xray_absorption(10-10000ev).txt")

# concatenate UV and visable data
data= data1[:19,:]
data = np.concatenate((data, data2[:254,:]))
data = np.concatenate((data, data4))


# processing data for x-ray
ev = xray[:,0]
evabs = xray[:,1]
λ = 1240/ev
P = 133.322 #pressure in SI unit (pascal)
V = 1 # assumed to be 1 m^3 does not matter
R = 8.314462 # Gas constant in SI unit
T = 293.15 # Temp in Kelivn
N = P*V/(R*T)*1e-6 # 1e-6 to get it in cm^3
C = N*6.02e23 # number of molecules equal avogadro number times number of moles
σ = -np.log(evabs)/(C) # absorption cross-scetion per molecule
xspl = InterpolatedUnivariateSpline(np.flip(λ, axis=0), np.flip(σ, axis=0))

# estimated IR cross-sections
λ_infra = np.linspace(1e-6,20e-6,4000)
infrared = gauss(λ_infra,20e-24,2.72e-6,0.04e-6)+ gauss(λ_infra,400e-24,3.3e-6,0.04e-6)+ gauss(λ_infra,3e-21,4.75e-6,0.05e-6)+gauss(λ_infra,40e-21,9.6e-6,0.2e-6)
infraspl = InterpolatedUnivariateSpline(λ_infra, infrared)

# get data ready for spline
x = data[:,0]
y = data[:,1]

# spline date
visdataspl = interp1d(x, y)

# create freq grid based on the number of points and the lower nm limit given above
freqdata = np.linspace(0, 2*np.pi*3e8/(nml*1e-9),npts)

# wavelength grid
xdata = 2*np.pi*3e8/freqdata*1e9

# create empty array to carry all absorption cross-section
full_absorption_data = np.array([])

# loop over the freq grid to add the cross-section from different ranges
for i in range(len(xdata)):
    if xdata[i]<=64:
        full_absorption_data = np.append(full_absorption_data, xspl(xdata[i]))
    elif xdata[i]>64 and xdata[i]<1100:
        full_absorption_data = np.append(full_absorption_data, visdataspl(xdata[i]))
    else:
        full_absorption_data = np.append(full_absorption_data, infraspl(xdata[i]*1e-9))


# here, we transfere the absorption cross-section coeffient to the imaginary part(κ(ω)) of index of refraction
imchi = 3e8*full_absorption_data*2.5e21/(2*freqdata)








os.chdir(fname)
np.savetxt("all.txt",data, fmt='%.5e')

#saving absorption cross-section data
np.savetxt("full_absorption_data.txt".format(fname),full_absorption_data[1:], fmt='%.5e')

np.savetxt("imchi.txt".format(fname),imchi[1:], fmt='%.5e')
np.savetxt("freqdata.txt".format(fname),freqdata[1:], fmt='%.5e')

plt.plot(xdata[xdata<1000], full_absorption_data[xdata<1000])
plt.yscale("log")
plt.ylim(1e-25,1e-16)

plt.savefig("image.png")
