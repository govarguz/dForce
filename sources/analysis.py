import math
import numpy as np
import forces



#==============Parameters which are fuction of time=============


def velocity(z,delta_t):

	derivate = np.zeros(len(z))
	for c in np.arange(len(z) - 1):
		derivate[c] = z[c+1]-z[c]
	derivate = derivate/delta_t

	return derivate

def force_ts(option_mask,z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,mutip,musample):
	f = np.zeros(len(z))
	for ii in range(len(z)):
		f[ii] = forces.sumatory(option_mask,z[ii],v[ii],t[ii],f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,mutip,musample) - forces.driven(F0,f0,t[ii])

	return f

def force_dr(F0,f0,t):
	f = np.zeros(len(t))
	for ii in range(len(t)):
		f[ii] = forces.driven(F0,f0,t[ii])

	return f

def identation(z,zc,a0):
	ident = np.zeros(len(z))
	for ii in range(len(z)):
		if (z[ii] + zc) < a0:
			ident[ii] = -(z[ii] + zc - a0)
		else:
			ident[ii] = 0

	return ident



# =============Parameters which are function of zc============

def contact_time(z,zc,a0):
	counter = 0.0
	for ii in range(len(z)):
		if (z[ii] + zc) < a0:
			counter = counter + 1.0
	time = counter/float(len(z))

	return time


def amplitude(z):
	amp = (np.amax(z) - np.amin(z))/2.0

	return amp

def max_force(f):
	max_for = np.amax(np.absolute(f))

	return max_for


def minimum_distance(z,zc,a0):
	minimum = np.amin(z + zc-a0)

	return minimum

def maximum_distance(z,zc):
	maximum = np.amax(z + zc)

	return maximum

def mean_deflection(z):
	maximum = np.amax(z)
	minimum = np.amin(z)
	mean = (maximum+minimum)/2.0

	return mean

def dissipated_energy(v,f):
	total = 0.0
	for ii in range(len(v)):
		total = total - v[ii]*f[ii]

	return total

def virial(z,f):
	total = 0.0
	for ii in range(len(z)-1):
		total = total + z[ii]*f[ii]

	return total


def dissipated_power(v,f):
	maxi = -1000.0
	for ii in range(len(v)-1):
		instant_pow = v[ii]*f[ii]
		if instant_pow > maxi:
			maxi = instant_pow

	return maxi



def mean_deflection(z):
	total = sum(z)
	total = total/float(len(z))

	return total


def phase_shift(z,f):
	f_trans = np.fft.fft(f)
	z_trans = np.fft.fft(z)
	z_amp = np.absolute(z_trans)
	f_amp = np.absolute(f_trans)
	f_ang = np.angle(f_trans)
	z_ang = np.angle(z_trans)

	index_f = np.argmax(f_amp)
	index_z = np.argmax(z_amp)
	shift = -z_ang[index_f] + f_ang[index_f]
	shift = np.arccos(np.cos(shift))*180/math.pi

	return shift


def fourier(z,a):
	z_trans = np.fft.fft(z)
	z_amp = np.absolute(z_trans)
	z_ang = np.angle(z_trans)
	z_ang = z_ang*180/math.pi
	freq = np.fft.fftfreq(len(z), d=a)

	return z_amp, z_ang, freq
