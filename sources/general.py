#!/usr/bin/python
#=======================dFORCE==============================
#The code written below is part of the dFORCE project and has
#been developed at the Forcetool laboratory in Madrid. This 
#software will enable students and researchers from the 
#AFM sphere to easily get an insight into the AFM theory 
#through numerical simulations.
#======================CREDITS==============================
#
# Prof. Ricardo Garcia
# PhD Horacio Vargas
# PhD student Pablo D Garcia
#
#We thank you in advance for sending your feedback and/or 
#suggestions in github.


#======================dForce LICENSE==============================

#Copyright (C) 2014  PhD Horacio V. Guzman and PhD student Pablo D Garcia

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
#@@@@@@@@@@@@@@@@@@@@@2014@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# What we observe is not nature itself, but nature exposed to
# our method of questioning. 
#                                  W. Heisenberg
############################################################
#===========================================================
#
#========================Code lines========================
# IMPORT MODULES
import os
import shutil
import gtk
import gtk.glade
import pygtk # temp
# IMPORT NEWLY VISUALIZATION MODULES
from scipy import *
#import mypypm1
#import mypyeb1
from pylab import *
import numpy as np

import forces
import analysis
import time



# IMPORT NEWLY VISUALIZATION MODULES
from scipy import *
from scipy import integrate

from pylab import *
import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__))




#========FUNCTIONS=================
#|||||||||||||||||||||||||||||||||||||||||||||||



# Function devised to calculate the numerical derivate of a vector (the last position of the output will be always zero)
def diff(vector):

	derivate = np.zeros(len(vector))
	for c in np.arange(len(vector) - 1):
		derivate[c] = vector[c+1]-vector[c]

	return derivate


# PUNTUAL MASS 
def deriv(y,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,option_mask):
	uprime = y[1]
	wprime = - k*y[0]/m -D*y[1]/m + forces.sumatory(option_mask,y[0], y[1],t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc)/m 
	yprime = np.array([uprime,wprime]) 

	return yprime



def runner(b, date,iteration,a):

	print("single")

    # Transofrm the readed data in the rigth units

	f0 = b[0]*1000    #Frequency of maximum response in amplitude tking in account the Q factor (the thermal tunning gives us that). The driving frequency will be this.
	k = b[2]
	Q = b[3]
	if dolly.glade_ui.get_widget("fcorrection_q1").get_active():
		f_ori = f0/math.sqrt(1-1/(2*Q**2))  # Frequency of resonance if there weren't viscosity (in the case Q --> infinity)
	else:
		f_ori = f0
	A01 = b[4]
	Rt = b[5]
	a0 = b[6]
	E =  b[7]
	H = b[8]
	eta = b[9]
	sigmas = b[10]
	sigmat = b[11]
	landadeb = b[12]
	epsilonr = b[13]
	alfa = b[14]
	nc = b[15] 
	m = k/(f_ori*f_ori*4*math.pi*math.pi)
	D = math.sqrt(m*k)/Q

	if dolly.glade_ui.get_widget("copy_fre1").get_active() == 0:
		fd = dolly.glade_ui.get_widget("fd1_sp").get_value()*1000  
		fd = np.array(fd)
	else:
		fd = f0

	F0 = m*A01*math.sqrt((f_ori**2-fd**2)**2 + (fd*f_ori/Q)**2)*4*math.pi*math.pi

	print("************")
	print((f_ori**2-fd**2)**2)
	print((fd*f_ori/Q)**2)
	print(math.sqrt((f_ori**2-fd**2)**2 + (fd*f_ori/Q)**2))
	print(A01*math.sqrt((f_ori**2-fd**2)**2 + (fd*f_ori/Q)**2))
	print(m*A01*math.sqrt((f_ori**2-fd**2)**2 + (fd*f_ori/Q)**2))
	print("-----------")	
	print(m)
	print(f0)
	print(f_ori)
	print(fd)
	print(F0)
	print("************")

	if dolly.glade_ui.get_widget("copy_fre1").get_active() == 0:     #The value of fo is, from here, the value of the driven frequency, not the resonant frequency     
		f0 = fd

	mu = 4  #dummy variable

	tolerance = b[23]

	start=0
	end=b[16]/f0
	numsteps=b[16]*b[17]
	numsteps_ss = b[17]*b[18]
	n_cut = numsteps -numsteps_ss
	t=np.linspace(start,end,numsteps)

	y01=np.array([0.0000,0.000000])   #Initial values

	#ns = math.floor((b[19] - b[20])/b[21])   #number of steps for zc

	#zc = np.linspace(b[20], b[19], ns) 

	zc = np.arange(b[19],b[20],-b[21])
	ns = zc.size

	z_output = np.zeros([ns,12])
	z_output[:,0] = zc*10**9
	t_output = np.zeros([numsteps_ss,8])
	t_output[:,0] = t[n_cut:].copy()*f0


	closer = np.argmin(np.fabs(b[22]*np.ones(len(zc)) - zc))   #Which step in zc is closer to the choosen zc to do ana anlysis

	print(f0)
	print(k)
	print(Q)
	print(A01)
	print(Rt)
	print(a0)
	print(E)
	print(H)
	print(eta)
	print(F0)
	print(float(b[16]))
	print(b[17])	
	print(b[18])


	dolly.abort = 0

	#  Here it starts the bucle for the differents values of zc


	for counter in np.linspace(0, ns-1, ns):


		dolly.glade_ui.get_widget("progressbar1").set_fraction(counter/ns)
		dolly.glade_ui.get_widget("progressbar1").set_text(str(np.around(100*counter/ns,0))+ " %")


		if dolly.abort == 1:
			break 

		option_mask = [dolly.glade_ui.get_widget("dlvo_check").get_active(), dolly.glade_ui.get_widget("vdw_check").get_active(), dolly.glade_ui.get_widget("hertz_check").get_active(), dolly.glade_ui.get_widget("dmt_check").get_active(), dolly.glade_ui.get_widget("tatara_check").get_active(), dolly.glade_ui.get_widget("viscosity_check").get_active()]


		while gtk.events_pending():
			gtk.main_iteration(False)

		print ("zc step: %i" %counter)

		a = t[1] - t[0]
		y, info = integrate.odeint(deriv,y01,t, args = (f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc[counter],eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,option_mask), full_output = 1, rtol = 10**-tolerance, atol = 10**-tolerance)

		stationary = y[n_cut:,0].copy()
		#v = analysis.velocity(stationary,a)
		v = y[n_cut:,1].copy()
		f_ts = analysis.force_ts(option_mask,stationary,v,t[n_cut:],f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc[counter],eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc)


		amplitude = analysis.amplitude(stationary)
		z_output[counter,1] = amplitude*10**9
		phase = analysis.phase_shift(stationary,analysis.force_dr(F0,f0,t[n_cut:]))
		z_output[counter,2] = phase
		minimum_dis = analysis.minimum_distance(stationary,zc[counter],a0)
		z_output[counter,3] = minimum_dis*10**9
		f_max = analysis.max_force(f_ts)
		z_output[counter,5] = f_max*10**12
		dis_en = analysis.dissipated_energy(v,f_ts)
		z_output[counter,6] = a*dis_en*10**20/(b[18])
		vir = analysis.virial(stationary,f_ts)
		z_output[counter,7] = vir*10**20/(b[18]*b[17])
		meand = analysis.mean_deflection(stationary)
		z_output[counter,9] = meand
		t_con = analysis.contact_time(stationary,zc[counter],a0)
		z_output[counter,10] = t_con
		dis_pow = analysis.dissipated_power(v,f_ts)
		z_output[counter,11] = dis_pow*10**16	

		if counter == closer:
			t_output[:,1] = y[n_cut:,0]*10**9
			t_output[:,2] = v*10**6
			t_output[:,3] = f_ts*10**12
			identation = analysis.identation(stationary,zc[counter],a0)
			t_output[:,4] = identation*10**9
			fampt, fphat, freq_bas = analysis.fourier(stationary,a)
			t_output[:,5] = fampt
			t_output[:,6] = fphat	
			t_output[:,7] = freq_bas


		y01=np.array([amplitude*cos(-phase*math.pi/180),amplitude*2*math.pi*f0*cos((-phase + 90)*math.pi/180)])


	name = "dforceproject"
	name  = name + date
	name = name  + "/tdom"
	name = name + iteration
	name = name + ".dfo"

	np.savetxt(name, t_output)

	name = "dforceproject"
	name  = name + date
	name = name  + "/zcdom"
	name = name + iteration
	name = name + ".dfo"

	np.savetxt(name, z_output)



	return t_output



# ==================== Euler-Bernoulii ================
def derivEB(y,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,option_mask):
	f = [y[1], - k[1]*y[0]/m[1] -D[1]*y[1]/m[1] + forces.sumatory(option_mask,y[0] +y[2], y[1] + y[3],t,f0,k[1],Q[1],Rt,a0,E,H,mu,m[1],D[1],F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc)/m[1] ,y[3], - k[2]*y[2]/m[2] -D[2]*y[3]/m[2] + forces.sumatory(option_mask,y[2]+y[0], y[3]+y[1],t,f0,k[2],Q[2],Rt,a0,E,H,mu,m[2],D[2],F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc)/m[2]]

	return f



def runnerEB(b, date,iteration):
	print("....EB....")

    # Transofrm the readed data in the rigth units

	f0 =  np.empty(3)
	f_ori = np.empty(3)
	k =  np.empty(3)
	Q =  np.empty(3)
	A0 =  np.empty(3)
	m =  np.empty(3)
	D =  np.empty(3)
	F0 =  np.empty(3)
	fd = np.empty(3)

	f0[1] = b[0]*1000
	f0[2] = b[1]*1000
	if dolly.glade_ui.get_widget("fcorrection_q1").get_active():
		f_ori[1] = f0[1]/math.sqrt(1-1/(2*Q[1]**2))  # Frequency of resonance if there weren't viscosity (in the case Q --> infinity)
	else:
		f_ori[1] = f0[1]
	if dolly.glade_ui.get_widget("fcorrection_q2").get_active():
		f_ori[2] = f0[2]/math.sqrt(1-1/(2*Q[2]**2))  # Frequency of resonance if there weren't viscosity (in the case Q --> infinity)
	else:
		f_ori[2] = f0[2]
	k[1] = b[2]
	k[2] = b[3]
	Q[1] = b[4]
	Q[2] = b[5]
	A0[1] = b[6]
	A0[2] =  b[7]
	Rt = b[8]
	a0 = b[9]
	E = b[10]
	H = b[11]
	eta = b[12]
	sigmas = b[13]
	sigmat = b[14]
	landadeb = b[15]
	epsilonr = b[16]
	alfa = b[17]
	nc = b[18]
	m[1] = k[1]/(f_ori[1]*f_ori[1]*4*math.pi*math.pi)
	m[2] = k[2]/(f_ori[2]*f_ori[2]*4*math.pi*math.pi)
	D[1] = math.sqrt(m[1]*k[1])/Q[1]
	D[2] = math.sqrt(m[2]*k[2])/Q[2]



	if dolly.glade_ui.get_widget("copy_fre1").get_active() == 0:
		fd[1] = dolly.glade_ui.get_widget("fd1_sp").get_value()*1000  
	else:
		fd[1] = f0[1]
	if dolly.glade_ui.get_widget("copy_fre2").get_active() == 0:
		fd[2] = dolly.glade_ui.get_widget("fd2_sp").get_value()*1000  
	else:
		fd[2] = f0[2]

	F0[1] = m[1]*A0[1]*math.sqrt((f_ori[1]**2-fd[1]**2)**2 + (fd[1]*f_ori[1]/Q[1])**2)*4*math.pi*math.pi
	F0[2] = m[2]*A0[2]*math.sqrt((f_ori[2]**2-fd[2]**2)**2 + (fd[2]*f_ori[2]/Q[2])**2)*4*math.pi*math.pi


	mu = 4  #dummy variable

	tolerance = b[26]


	print(f0)
	print(k)
	print(Q)
	print(A0)
	print(Rt)
	print(a0)
	print(E)
	print(H)
	print(eta)
	print(F0)
	print("************")



	start=0
	end=b[19]/f0[2]
	numsteps=b[19]*b[20]
	numsteps_ss = b[20]*b[21]
	n_cut = numsteps -numsteps_ss
	t=np.linspace(start,end,numsteps)

	y01=np.array([0.0000,0.000000,0.0000,0.00000])    #Initial values

	#ns = math.floor((b[22] - b[23])/b[24])   #number of steps for zc

	#zc = np.linspace(b[23], b[22], ns) 


	zc = np.arange(b[22],b[23],-b[24])
	ns = zc.size

	z_output = np.zeros([ns,17])
	z_output[:,0] = zc*10**9
	t_output = np.zeros([numsteps_ss,21])
	t_output[:,0] = t[n_cut:].copy()*f0[1]


	closer = np.argmin(np.fabs(b[25]*np.ones(len(zc)) - zc))   #Which step in zc is closer to the choosen zc to do ana anlysis

	print(end)
	print(numsteps)
	print(numsteps_ss)
	print(n_cut)
	print(ns)
	print(zc)
	print("************")


	dolly.abort = 0

	#  Here it starts the bucle for the differents values of zc


	for counter in np.linspace(0, ns-1, ns):


		if dolly.abort == 1:
			break 

		dolly.glade_ui.get_widget("progressbar1").set_fraction(counter/ns)
		dolly.glade_ui.get_widget("progressbar1").set_text(str(np.around(100*counter/ns,0))+ " %")

		option_mask = [dolly.glade_ui.get_widget("dlvo_check").get_active(), dolly.glade_ui.get_widget("vdw_check").get_active(), dolly.glade_ui.get_widget("hertz_check").get_active(), dolly.glade_ui.get_widget("dmt_check").get_active(), dolly.glade_ui.get_widget("tatara_check").get_active(), dolly.glade_ui.get_widget("viscosity_check").get_active()]


		while gtk.events_pending():
			gtk.main_iteration(False)


		print ("zc step: %i" %counter)

		a = t[1] - t[0]
		y = integrate.odeint(derivEB,y01,t, args = (f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc[counter],eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,option_mask),rtol = 10**-tolerance, atol = 10**-tolerance)

		stationary1 = y[n_cut:,0].copy()
		stationary2 = y[n_cut:,2].copy()
		stationary_sum = stationary1 + stationary2	


		v1 = analysis.velocity(stationary1,a)
		v2 = analysis.velocity(stationary2,a)
		v_sum = analysis.velocity(stationary_sum,a)

		f_ts = analysis.force_ts(option_mask,stationary_sum,v_sum,t[n_cut:],f0,k[1],Q[1],Rt,a0,E,H,mu,m[1],D[1],F0,zc[counter],eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc)


		amplitude1 = analysis.amplitude(stationary1)
		z_output[counter,1] = amplitude1*10**9
		phase1 = analysis.phase_shift(stationary1,analysis.force_dr(F0[1],f0[1],t[n_cut:]))
		z_output[counter,2] = phase1
		amplitude2 = analysis.amplitude(stationary2)
		z_output[counter,3] = amplitude2*10**9
		phase2 = analysis.phase_shift(stationary2,analysis.force_dr(F0[2],f0[2],t[n_cut:]))
		z_output[counter,4] = phase2
		minimum_dis = analysis.minimum_distance(stationary_sum,zc[counter],a0)
		z_output[counter,5] = minimum_dis*10**9
		maximum_dis = analysis.maximum_distance(stationary_sum,zc[counter])
		z_output[counter,6] = maximum_dis*10**9		
		f_max = analysis.max_force(f_ts)
		z_output[counter,7] = f_max*10**12
		edisp1 = analysis.dissipated_energy(v1,f_ts)
		z_output[counter,8] = a*b[20]*edisp1*10**20/(b[21]*f0[1]/f0[2])
		edisp2 = analysis.dissipated_energy(v2,f_ts)
		z_output[counter,9] = a*b[20]*edisp2*10**20/(b[21]*b[20])
		viri1 = analysis.virial(stationary1,f_ts)
		z_output[counter,10] = viri1*10**20/(b[21]*f0[1]/f0[2])
		viri2 = analysis.virial(stationary2,f_ts)
		z_output[counter,11] = viri2*10**20/(b[21]*b[20])	
		meand = analysis.mean_deflection(stationary_sum)
		z_output[counter,13] = meand
		t_con = analysis.contact_time(stationary_sum,zc[counter],a0)
		z_output[counter,14] = t_con
		pow1 = analysis.dissipated_power(v1,f_ts)
		z_output[counter,15] = pow1*10**16
		pow2 = analysis.dissipated_power(v2,f_ts)
		z_output[counter,16] = pow2*10**16		

		if counter == closer:

			t_output[:,1] = stationary_sum.copy()*10**9
			t_output[:,2] = v_sum.copy()*10**6
			t_output[:,3] = stationary1.copy()*10**9
			t_output[:,4] = v1.copy()*10**6
			t_output[:,5] = stationary2.copy()*10**9
			t_output[:,6] = v2.copy()*10**6
			t_output[:,12] = f_ts.copy()*10**12
			identation = analysis.identation(stationary_sum,zc[counter],a0)
			t_output[:,13] = identation.copy()*10**9
			famp1, fpha1, freq_bas = analysis.fourier(stationary1,a)
			t_output[:,14] = famp1.copy()
			t_output[:,15] = fpha1.copy()	
			famp2, fpha2, freq_bas = analysis.fourier(stationary2,a)
			t_output[:,16] = famp2.copy()
			t_output[:,17] = fpha2.copy()
			fampt, fphat, freq_bas = analysis.fourier(stationary_sum,a)
			t_output[:,18] = fampt.copy()
			t_output[:,19] = fphat.copy()		
			t_output[:,20] = freq_bas.copy()			

		y01=np.array([amplitude1*cos(-phase1*math.pi/180),amplitude1*2*math.pi*f0[1]*cos((-phase1 + 90)*math.pi/180),amplitude2*cos(-phase2*math.pi/180),amplitude2*2*math.pi*f0[2]*cos((-phase2 + 90)*math.pi/180)])


	name = "dforceproject"
	name  = name + date
	name = name  + "/tdom"
	name = name + iteration
	name = name + ".dfo"

	np.savetxt(name, t_output)

	name = "dforceproject"
	name  = name + date
	name = name  + "/zcdom"
	name = name + iteration
	name = name + ".dfo"

	np.savetxt(name, z_output)



	return t_output



#==========Initial problem conditions===============
x=0
t=0
rk4order=2
#==========Initialize output vectors================
#y=zeros(nmax)
#yout=zeros(nmax)
#|||||||||||||||||||||||||||||||||||||||||||||||||||||
class dforce:
# Define modules that will be part of the Dast class
# + Module init initializes the
    def __init__(self):    # defines the module with the input parameter self


        self.glade_ui = gtk.glade.XML('dforceInterfaceV1.glade')    # defines the var self.glade_ui as the same in the est module
        self.glade_ui.signal_autoconnect(self)    # defines the var self.glade_ui as the same in the est module gtk.glade class
        self.window1 = self.glade_ui.get_widget('window1')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.savingd = self.glade_ui.get_widget('savingd')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.uploadingd = self.glade_ui.get_widget('uploadingd')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.credits = self.glade_ui.get_widget('credits')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.vdw_par = self.glade_ui.get_widget('vdw_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	 
        self.dlvo_par = self.glade_ui.get_widget('dlvo_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	        
        self.dmt_par = self.glade_ui.get_widget('dmt_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.tatara_par = self.glade_ui.get_widget('tatara_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	 
        self.hertz_par = self.glade_ui.get_widget('hertz_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	        
        self.jkr_par = self.glade_ui.get_widget('jkr_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.viscosity_par = self.glade_ui.get_widget('viscosity_par')    # defines the var self.glade_ui as the same in the est module gtk.glade class	 
        self.value_table = self.glade_ui.get_widget('value_table')

        # Hide the window when the "x" button is pressed
        def hide_window(window, event):
	        window.hide()
	        return True
	    # Disable the "x" button 
        def no_exit_window(window, event):
	        return True
        self.vdw_par.connect('delete-event', hide_window)
        self.window1.connect('delete-event', no_exit_window)
        self.savingd.connect('delete-event', hide_window)
        self.uploadingd.connect('delete-event', hide_window)
        self.credits.connect('delete-event', no_exit_window)
        self.dlvo_par.connect('delete-event', hide_window)
        self.dmt_par.connect('delete-event', hide_window)
        self.tatara_par.connect('delete-event', hide_window)
        self.hertz_par.connect('delete-event', hide_window)
        self.jkr_par.connect('delete-event', hide_window)
        self.viscosity_par.connect('delete-event', hide_window)
        self.value_table.connect('delete-event', hide_window)


        self.window1.show()

#-----------------------Dissabling buttons in general------------------------------
        #self.glade_ui.get_widget("bmafm_cb").set_sensitive(False)
        #self.glade_ui.get_widget("eb_cb").set_sensitive(False)
        self.glade_ui.get_widget("kc2_sp").set_sensitive(False)
        self.glade_ui.get_widget("f02_sp").set_sensitive(False)
        self.glade_ui.get_widget("q2_sp").set_sensitive(False)
        self.glade_ui.get_widget("jkr_check").set_sensitive(False)
        self.glade_ui.get_widget("zdomload_sp1").set_sensitive(False)
        self.glade_ui.get_widget("tdomload_sp1").set_sensitive(False)
        self.glade_ui.get_widget("Adhesion_check").set_sensitive(False)
        self.glade_ui.get_widget("BECC_check").set_sensitive(False)
        self.glade_ui.get_widget("SLS_check").set_sensitive(False)
        self.glade_ui.get_widget("etst_cb").set_sensitive(False)
        self.glade_ui.get_widget("ptst_cb").set_sensitive(False)

        self.glade_ui.get_widget("pts_x1").set_sensitive(False)
        self.glade_ui.get_widget("pts_y1").set_sensitive(False)
        self.glade_ui.get_widget("pts_x2").set_sensitive(False)
        self.glade_ui.get_widget("pts_y2").set_sensitive(False)
        self.glade_ui.get_widget("pts_x3").set_sensitive(False)
        self.glade_ui.get_widget("pts_y3").set_sensitive(False)
        self.glade_ui.get_widget("ets_x1").set_sensitive(False)
        self.glade_ui.get_widget("ets_y1").set_sensitive(False)
        self.glade_ui.get_widget("ets_x2").set_sensitive(False)
        self.glade_ui.get_widget("ets_y2").set_sensitive(False)
        self.glade_ui.get_widget("ets_x3").set_sensitive(False)
        self.glade_ui.get_widget("ets_y3").set_sensitive(False)                  

        self.glade_ui.get_widget("amps2_cb").set_sensitive(False)
        self.glade_ui.get_widget("fases2_cb").set_sensitive(False)
        self.glade_ui.get_widget("ets2_cb").set_sensitive(False)   
        self.glade_ui.get_widget("pts2_cb").set_sensitive(False)
        self.glade_ui.get_widget("virial2_cb").set_sensitive(False)
        self.glade_ui.get_widget("z1_cb").set_sensitive(False)
        self.glade_ui.get_widget("z2_cb").set_sensitive(False)    
        self.glade_ui.get_widget("z3_cb").set_sensitive(False)
        self.glade_ui.get_widget("z4_cb").set_sensitive(False)  

        self.glade_ui.get_widget("phi2a2_cb").set_sensitive(False)
        self.glade_ui.get_widget("phi2a1_cb").set_sensitive(False)
        self.glade_ui.get_widget("phi1a2_cb").set_sensitive(False)    
        self.glade_ui.get_widget("famp2").set_sensitive(False)
        self.glade_ui.get_widget("fpha2").set_sensitive(False)
        self.glade_ui.get_widget("famp1").set_sensitive(False)
        self.glade_ui.get_widget("fpha1").set_sensitive(False)          
        self.glade_ui.get_widget("fcorrection_q2").set_sensitive(False)
        self.glade_ui.get_widget("a02_sp").set_sensitive(False)
        self.glade_ui.get_widget("fd2_sp").set_sensitive(False)    
        self.glade_ui.get_widget("copy_fre2").set_sensitive(False)



#-----------------------------------------------------------------------------------

	self.flagin=0
	self.flagpm=1
	self.flagfile=0
	self.flagpnam=0
	self.flagaux=0
	self.compareflag = 0     # Flag for compare two different simulations
	self.date = time.strftime("%m,%d-%H,%M")
	self.abort = 0
	#print(self.date)
#========checkbutton get active=================
#|||||||||||||||||||||||||||||||||||||||||||||||
    def tucom(self,cb):				# New
	return self.glade_ui.get_widget(cb).get_active()
#||||||||||||||||||||||||||||||||||||||||||||||||||
#=========Buttons==================================
#||||||||||||||||||||||||||||||||||||||||||||||||||
#=================Quit option========================
    def gtk_main_quit(self, widget):
	gtk.main_quit()
#=================RUN========================
#||||||||||||||||||||||||||||||||||||||||||||||||||||
    def on_run_ex_clicked(self, widget):
	self.glade_ui.get_widget("abort").set_visible(True)
	self.savingd.show()

	while gtk.events_pending():
		gtk.main_iteration(False)   #Refresh the GUI

	#self.flagin=0	# Runs per windows ++
	# recently added code
	if (self.flagaux==0):	# Did you uplaod inputs in Bulk none?
		if (self.glade_ui.get_widget("single_cb").get_active()==1) and ((self.glade_ui.get_widget("bmafm_cb").get_active()==1) or (self.glade_ui.get_widget("pm_cb").get_active()==1)):
			# Point Mass, right?
			self.flagin+=1
			self.flagpm=1
			self.flagfile = 0   # There will not print the load files
			omega0a = self.glade_ui.get_widget("f01_sp").get_value()
			omegaa = self.glade_ui.get_widget("fd1_sp").get_value()
			q= self.glade_ui.get_widget("q1_sp").get_value()	
			kc = self.glade_ui.get_widget("kc1_sp").get_value()
			a0a = self.glade_ui.get_widget("a01_sp").get_value()
			a0 = a0a*1e-9
			emuestraa = self.glade_ui.get_widget("saym_sp").get_value()
			emuestra = emuestraa*1000000.	# in Mpa now
			epuntaa = self.glade_ui.get_widget("ymt_sp").get_value()
			epunta = epuntaa*1000000000.	# in Gpa now
			rada = self.glade_ui.get_widget("tr_sp").get_value()
			rad = rada*1e-9
			rs = self.glade_ui.get_widget("sr_sp").get_value()  #Surface radius
			rs = rs*1e-9
			a00a = self.glade_ui.get_widget("a00_sp").get_value()
			a00 = a00a*1e-9
			hama = self.glade_ui.get_widget("ham_sp").get_value()
			ham = hama*1e-20
			eta = self.glade_ui.get_widget("mvis_sp").get_value()    # Viscosity coeficient
			sigmasa = self.glade_ui.get_widget("sigmas_sp").get_value()  # Surface charge density
			sigmas = sigmasa*10**-3
			sigmata = self.glade_ui.get_widget("sigmat_sp").get_value()  # Tip charge density
			sigmat = sigmata*10**-3
			delta = 0.0
			landadeba = self.glade_ui.get_widget("landadeb_sp").get_value()  # landa debye (nm)
			landadeb = landadeba*10**-9
			epsilon = self.glade_ui.get_widget("epsilon_sp").get_value()   #epsilon
			mutip = self.glade_ui.get_widget("mutip_sp").get_value()     #poisson tip
			musample = self.glade_ui.get_widget("musam_sp").get_value()  #poisson sample
			nper = self.glade_ui.get_widget("nper_sp").get_value()
			npp = self.glade_ui.get_widget("npp_sp").get_value()
			naux=int(nper*npp) 				# need C
			nperfin = self.glade_ui.get_widget("nperfin_sp").get_value()
			dzca = self.glade_ui.get_widget("deltazc_sp").get_value()
			dzc = dzca*1.e-9
			zcmina = self.glade_ui.get_widget("zcmin_sp").get_value()
			zcmin = zcmina*1.e-9
			zcmaxa = self.glade_ui.get_widget("zcmax_sp").get_value()
			zcmax = zcmaxa*1.e-9
			fixzca = self.glade_ui.get_widget("fixedzc_sp").get_value()
			tolerance = self.glade_ui.get_widget("tolerance_sp").get_value()
			if fixzca==0.0:
				fixzc = zcmax
			else:
				fixzc = fixzca*1.e-9
			if emuestra > 0 and epunta > 0:		# NEW
				ebarra=(1./((1.-musample**2)/emuestra+(1.-mutip**2)/epunta))
			else:
				ebarra=0.

			#-----------------tatara----------------------				
			if rs > 0 and self.glade_ui.get_widget("tatara_check").get_active():
				radeff=(1./(1./rs+1./rad))
				alfa = 4.0/3.0*ebarra*math.sqrt(radeff)
				nc = 4*math.pi*rs*rad*emuestra*epunta/(6+mutip+musample-2*mutip**2-2*musample**2)
				rad = radeff
			else:
				alfa = 0
				nc = 0




			#=========================Files handling w Python ===============================
			#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			if self.flagpnam==1 and self.flagin==1:
				self.palo1="dforceproject" + self.date
				os.mkdir(self.palo1)
			elif self.flagpnam==0 and self.flagin==1:
				self.palo1="dforceproject" + self.date
				os.mkdir(self.palo1)
			elif self.flagin>1:
				print "Project folder has been already created"

			date = self.date
			iteration = str(self.flagin)

			#=========================== Saving Sim Inputs ===============================
			#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			tx00=[ '"res freq 1 (kHz)"', '"drive freq 2 (kHz)"', '"k1 (N/m)"', '"Q1 (adim)"', '"A01 (V)"', '"Tip Radius (nm)"', '"a0 (nm)"', '"Reduced YM (Pa)"', '"Hamaker Const(J)"', '"Viscosity (Pa s)"','"sigma surface (adim)"', '"sigma tip(adim)"', '"Debye length"', '"Epsilon"', '"alfa"', '"nc"', '"N periods osc"','"N point period"', '"N per to the ss"', '"zc max (nm)"', '"zc min (nm)"', '"delta zc (nm)"', '"zc fixed td(nm)"','"tolerance"' ]
			tx00=np.array(tx00)
			tx01=[ omega0a, omegaa, kc, q, a0, rad, a00, ebarra, ham, eta, sigmas, sigmat, landadeb, epsilon, alfa, nc, nper, npp, nperfin, zcmax, zcmin, dzc, fixzc,tolerance]	# Converted to an unique array
			tx01=np.array(tx01)
			f1 = open("dforceproject"+ date + "/inputs"+ iteration + ".ini","w+")
			for i in np.ndindex(tx00.shape):
				filerow='%s      %5.5e\n' % (tx00[i],tx01[i])
				f1.write(filerow)
			f1.close()
#	elif: ADD FILE HAS BEEN LOADED IN ADVANCE
#=========================Call fortran functions ===============================
			print naux
			#mypypm1.mypmmoda.mainbim(naux)  # Calling the function mainbim F90 PM
			print "Runnig f90 module PM"


			t_output = runner(tx01,date,iteration,self)

			self.glade_ui.get_widget("abort").set_visible(False)
			self.savingd.hide()


#=========================Files handling w Python ===============================
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

			self.piat=os.getcwd()
			self.naind="inputs"+str(self.flagin)+".ini"
			self.nazcd="zcdom"+str(self.flagin)+".dfo"
			self.natd="tdom"+str(self.flagin)+".dfo"
			#os.rename("inputs.ini",self.naind)
			#os.rename("zcdom.dfo",self.nazcd)
			#os.rename("tdom.dfo",self.natd)
			#os.mkdir(self.palo1)
			#shutil.copy(self.naind,self.piat+"/"+self.palo1)
			#shutil.copy(self.nazcd,self.piat+"/"+self.palo1)
			#shutil.copy(self.natd,self.piat+"/"+self.palo1)
			#os.remove(self.naind)
			#os.remove(self.nazcd)
			#os.remove(self.natd)
			#os.chdir(self.piat+"/"+self.palo1)
#========================Obtaining outputs zc domain ==========================
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
		elif self.glade_ui.get_widget("bmafm_cb").get_active()==1 and self.glade_ui.get_widget("eb_cb").get_active()==1:
			self.flagin+=1
			self.flagpm=0
			# New default variables loading
			#mypyeb1.myebmoda.invols=[1.0, 0.288]
			#mypyeb1.myebmoda.temp=0.0
			#mypyeb1.myebmoda.nombre1=self.glade_ui.get_widget("nombre1_txt").get_value()
			#nomaux2="zcdom"+str(flagin)
			#mypyeb1.myebmoda.nombre1=array('nomaux2')
			#nomaux="tdom"+str(flagin)
			#mypyeb1.myebmoda.nombre2=array('nomaux')
			# Previous vars and parameters loading
			if (self.flagfile==0):
				omega0a = self.glade_ui.get_widget("f01_sp").get_value()
				omega02a = self.glade_ui.get_widget("f02_sp").get_value()
				omegaa = self.glade_ui.get_widget("fd1_sp").get_value()
				q= self.glade_ui.get_widget("q1_sp").get_value()
				q2= self.glade_ui.get_widget("q2_sp").get_value()
				kc= self.glade_ui.get_widget("kc1_sp").get_value()
				kc2= self.glade_ui.get_widget("kc2_sp").get_value()
				a0a = self.glade_ui.get_widget("a01_sp").get_value()
				a0 = a0a*1e-9  #OUT
				a02a = self.glade_ui.get_widget("a02_sp").get_value()
				a02 = a02a*1e-9 #OUT
				emuestraa = self.glade_ui.get_widget("saym_sp").get_value()
				emuestra =emuestraa*1000000.	# in Mpa now
				epuntaa = self.glade_ui.get_widget("ymt_sp").get_value()	
				epunta=epuntaa*1000000000.	# in Gpa now
				rada = self.glade_ui.get_widget("tr_sp").get_value()
				rad = rada*1e-9
				rs = self.glade_ui.get_widget("sr_sp").get_value()  #Surface radius
				rs = rs*1e-9
				a00a = self.glade_ui.get_widget("a00_sp").get_value()
				a00 = a00a*1e-9
				hama = self.glade_ui.get_widget("ham_sp").get_value()
				ham = hama*1e-20
				eta = self.glade_ui.get_widget("mvis_sp").get_value()
				epsilon1 = self.glade_ui.get_widget("lrd_sp").get_value()	
				epsilon2 = self.glade_ui.get_widget("cdis_sp").get_value()	
				#delta = 0.0
				ljmina = self.glade_ui.get_widget("ljdep_sp").get_value()		
				ljmin = ljmina*1e-9
				lengtha = self.glade_ui.get_widget("ljlen_sp").get_value()	
				length = lengtha*1e-9
				mtipa = self.glade_ui.get_widget("mmtip_sp").get_value()
				mtip = mtipa*1e-15
				msample = self.glade_ui.get_widget("mmsam_sp").get_value()
				nper = self.glade_ui.get_widget("nper_sp").get_value()
				npp =self.glade_ui.get_widget("npp_sp").get_value()
				naux=int(nper*npp) 				# need C
				eta = self.glade_ui.get_widget("mvis_sp").get_value()    # Viscosity coeficient
				sigmasa = self.glade_ui.get_widget("sigmas_sp").get_value()  # Surface charge density
				sigmas = sigmasa*10**-3
				sigmata = self.glade_ui.get_widget("sigmat_sp").get_value()  # Tip charge density
				sigmat = sigmata*10**-3
				delta = 0.0
				landadeba = self.glade_ui.get_widget("landadeb_sp").get_value()  # landa debye (nm)
				landadeb = landadeba*10**-9
				epsilon = self.glade_ui.get_widget("epsilon_sp").get_value()   #epsilon
				mutip = self.glade_ui.get_widget("mutip_sp").get_value()     #poisson tip
				musample = self.glade_ui.get_widget("musam_sp").get_value()  #poisson sampl
				nperfin = self.glade_ui.get_widget("nperfin_sp").get_value()
				dzca = self.glade_ui.get_widget("deltazc_sp").get_value()
				dzc= dzca*1.e-9
				zcmina = self.glade_ui.get_widget("zcmin_sp").get_value()
				zcmin = zcmina*1.e-9
				zcmaxa = self.glade_ui.get_widget("zcmax_sp").get_value()
				zcmax = zcmaxa*1.e-9
				fixzca = self.glade_ui.get_widget("fixedzc_sp").get_value()
				tolerance = self.glade_ui.get_widget("tolerance_sp").get_value()
				if fixzca==0.0:
					fixzc = zcmax
				else:
					fixzc = fixzca*1.e-9
				if emuestra > 0 and epunta > 0:		# E effective
					ebarra=(1./((1.-musample**2)/emuestra+(1.-mutip**2)/epunta))
				else:
					ebarra=0.

				#-----------------tatara----------------------				
				if rs > 0 and self.glade_ui.get_widget("tatara_check").get_active():
					radeff=(1./(1./rs+1./rad))
					alfa = 4.0/3.0*ebarra*math.sqrt(radeff)
					nc = 4*math.pi*rs*rad*emuestra*epunta/(6+mutip+musample-2*mutip**2-2*musample**2)
					rad = radeff
				else:
					alfa = 0
					nc = 0


				#=========================Files handling w Python ===============================
				#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
				if self.flagpnam==1 and self.flagin==1:
					self.palo1="dforceproject" + self.date
					os.mkdir(self.palo1)
				elif self.flagpnam==0 and self.flagin==1:
					self.palo1="dforceproject" + self.date
					os.mkdir(self.palo1)
				elif self.flagin>1:
					print "Project folder has been already created"

				date = self.date
				iteration = str(self.flagin)

				#=========================== Saving Sim Inputs ===============================
				#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
				tx00=[ '"res freq 1 (kHz)"', '"res freq 2 (kHz)"', '"k1 (N/m)"', '"k2 (N/m)"', '"Q1 (adim)"', '"Q2 (adim)"', '"A01 (V)"', '"A02 (V)"', '"Tip Radius (nm)"', '"a0 (nm)"', '"Reduced YM (Pa)"', '"Hamaker Const(J)"', '"Viscosity (Pa s)"','"sigma surface (adim)"', '"sigma tip(adim)"', '"LJdepth (nN)"', '"LJlength (nm)"', '"Tip magn (A m2)"', '"Surf magn (A m2)"', '"N periods osc"','"N point period"', '"N per to the ss"', '"zc max (nm)"', '"zc min (nm)"', '"delta zc (nm)"', '"zc fixed td (nm)"' ]
				tx00=np.array(tx00)
				tx01=[ omega0a, omega02a, kc, kc2, q, q2, a0, a02, rad, a00, ebarra, ham, 	eta, sigmas, sigmat, landadeb, epsilon, alfa, nc, nper, npp, nperfin, zcmax, zcmin, dzc,fixzc,tolerance ]	# Converted to an unique array
				tx01=np.array(tx01)
				f1 = open("dforceproject"+ date + "/inputs"+ iteration + ".ini","w+")	#np.savetxt('timedom'+str(zco)+'.out', tx1)
				for i in np.ndindex(tx00.shape):
					filerow='%s      %5.5e\n' % (tx00[i],tx01[i])
					f1.write(filerow)
				f1.close()
			else:
				print "Mete patinya"	# only for EB with inputs
			#============== Define a direct way of obtaining the naux=npp*nper ================
				datain = np.genfromtxt(inputs.ini)	# new
				din = datain[:,1]
				naux=din[20]*din[21]
#=========================Call fortran functions EX ===============================
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||	
			print naux
			#mypyeb1.myebmoda.mainbim(naux)  # Calling the function mainbim F90 EB
			print "Runnig f90 module EB"
			'''
			if self.flagpnam==1 and self.flagin==1:
				self.palo1=self.pnamed
				os.mkdir(self.palo1)
			elif self.flagpnam==0 and self.flagin==1:
				self.palo1="dforceproject" + self.date
				os.mkdir(self.palo1)
			elif self.flagin>1:
				print "Project folder has been already created"
			'''


			date = self.date
			iteration = str(self.flagin)
			t_output = runnerEB(tx01,date,iteration)


			self.glade_ui.get_widget("abort").set_visible(False)
			self.savingd.hide()

#=========================Files handling w Python ===============================
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			self.piat=os.getcwd()
			self.naind="inputs"+str(self.flagin)+".ini"
			self.nazcd="zcdom"+str(self.flagin)+".dfo"
			self.natd="tdom"+str(self.flagin)+".dfo"
			#os.rename("inputs.ini",self.naind)
			#os.rename("zcdom.dfo",self.nazcd)
			#os.rename("tdom.dfo",self.natd)
			#os.mkdir(self.palo1)
			#shutil.copy(self.naind,self.piat+"/"+self.palo1)
			#shutil.copy(self.nazcd,self.piat+"/"+self.palo1)
			#shutil.copy(self.natd,self.piat+"/"+self.palo1)
			#os.remove(self.naind)
			#os.remove(self.nazcd)
			#os.remove(self.natd)
			#os.chdir(self.piat+"/"+self.palo1)
		else:
			#self.flagin=0
			print "Mist hyper Mist"

#------------------------------------------------------------------------------------------------------------------------------------	
    


    def on_plot_ex_clicked(self, widget):

	palette=['b.--','r.--']
	matplotlib.rcParams.update({'font.size': 16,'figure.figsize' : [11, 10]}) # http://matplotlib.org/users/customizing.html
	
	for comp_counter in range(self.compareflag+1):
		if self.flagpm==1:
			#========Reading files for point-mass outputs===========
			#||||||||||||||||||||||||||||||||||||||||||
			if self.flagfile == 0:                 # Is there any file load?
				os.chdir(self.piat+"/"+self.palo1)	
				tdom = np.genfromtxt(self.natd)
				zcdom = np.genfromtxt(self.nazcd)
			else:
				directionstdom = [self.tdomload, self.tdomload1]
				directionszdom = [self.zdomload, self.zdomload1]
				tdom = np.genfromtxt(directionstdom[comp_counter])
				zcdom = np.genfromtxt(directionszdom[comp_counter])				
			tt = tdom[:,0]
			zt = tdom[:,1]
			vt = tdom[:,2]
			forcet = tdom[:,3]
			ident = tdom[:,4]
			famptex = tdom[:,5] 
			fphatex = tdom[:,6] 	
			freq_basex = tdom[:,7] 
			zc = zcdom[:,0]
			amp1d = zcdom[:,1]
			fase1d = zcdom[:,2]
			dmind = zcdom[:,3]
			dmaxd = zcdom[:,4] # new
			fmaxd = zcdom[:,5]
			ets1d = zcdom[:,6]
			vts1d = zcdom[:,7]
			fmedd = zcdom[:,8] # new
			defld = zcdom[:,9]
			tcd = zcdom[:,10]# new
			pts1d = zcdom[:,11] # new		
		else:
			os.chdir(self.piat+"/"+self.palo1)	# new trying
			#========Reading files for Euler Bernoulli outputs===========
			#||||||||||||||||||||||||||||||||||||||||||
			tdom = np.genfromtxt(self.natd)
			tt = tdom[:,0]
			zt = tdom[:,1]
			vt = tdom[:,2]
			z1 = tdom[:,3]
			v1 = tdom[:,4]
			z2 = tdom[:,5]
			v2 = tdom[:,6]
			z3 = tdom[:,7]
			v3 = tdom[:,8]
			z4 = tdom[:,9]
			v4 = tdom[:,10]
			zcref = tdom[:,11]
			forcet = tdom[:,12]
			ident = tdom[:,13]
			famp1ex = tdom[:,14]
			fpha1ex =tdom[:,15]
			famp2ex = tdom[:,16]
			fpha2ex =tdom[:,17]
			famptex = tdom[:,18]
			fphatex =tdom[:,19]
			freq_basex = tdom[:,20]
			# Bimodal data
			print "Creating vectors from files for Bimodal ARRAYS ZC distance"
			#========Material1=================
			zcdom = np.genfromtxt(self.nazcd)
			zc = zcdom[:,0]
			amp1d = zcdom[:,1]
			fase1d = zcdom[:,2]
			amp2d = zcdom[:,3]
			fase2d = zcdom[:,4]
			dmind = zcdom[:,5]
			dmaxd = zcdom[:,6]
			fmaxd = zcdom[:,7]
			ets1d = zcdom[:,8]
			ets2d = zcdom[:,9]
			pts1d = zcdom[:,15]	# new
			pts2d = zcdom[:,16]	# new
			vts1d = zcdom[:,10]
			vts2d = zcdom[:,11]
			fmedd = zcdom[:,12]
			defld = zcdom[:,13]
			tcd = zcdom[:,14]
		#==================================================
		# End of getting arrays from simulated DATA
		#======================PLOT========================
		if self.flagin!=0 or self.flagfile == 1:

			#========Building the big Switch===========
			#||||||||||||||||||||||||||||||||||||||||||
			def amp1():
			#========Amp 1 vs. zc=================
				figure(1)
				grid(False)
				#title('$A_1$($z_c$)')
				ylabel('$A_1 \; (nm)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,amp1d,palette[comp_counter])
				savefig('amp1Vzc'+str(self.flagin)+'.png')
			def pha1():
				#========Phase 1 vs. zc=================
				figure(2)
				grid(False)
				#title('$\phi_1$($z_c$)')
				ylabel('$\phi_1 \; (deg)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,fase1d,palette[comp_counter])
				savefig('phi1Vzc'+str(self.flagin)+'.png')
			def den1():
				#========Diss Energy 1st vs. zc=================
				figure(3)
				grid(False)
				#title('First disp. energy vs. Average distance')
				ylabel('$E_{ts1} * 10^{-20} \; (J)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,ets1d,palette[comp_counter])
				savefig('ets1Vzc'+str(self.flagin)+'.png')
			def dpo1():
				#========Diss Power 1st vs. zc=================
				figure(4)
				grid(False)
				#title('First disp. power vs. Average distance')
				ylabel('$P_{ts1} * 10^{-20} \; (W)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,pts1d,palette[comp_counter])
				savefig('pts1Vzc'+str(self.flagin)+'.png')
			def vir1():
				#========Virial 1st vs. zc=================
				figure(5)
				grid(False)
				#title('Virial first vs. Average distance')
				ylabel('$V_{ts1} * 10^{-20} \; (J)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,vts1d,palette[comp_counter])
				savefig('vts1Vzc'+str(self.flagin)+'.png')
			def def1():
				#========Deflection vs. zc=================
				figure(6)
				grid(False)
				#title('Deflection vs. Average distance')
				ylabel('$Defl \; (nm)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,defld,palette[comp_counter])
				savefig('deflVzc'+str(self.flagin)+'.png')
			def dmin1():
				#========Dmin vs. zc=================
				figure(7)
				grid(False)
				#title('Minimum Distance vs. Average distance')
				ylabel('$dmin \; (nm)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,dmind,palette[comp_counter])
				savefig('dminVzc'+str(self.flagin)+'.png')
			def maf():
				#========Fmax vs. zc=================
				figure(8)
				grid(False)
				#title('Maximum Force vs. Average distance')
				ylabel('$F_{max} \; (nN)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,fmaxd,palette[comp_counter])
				savefig('fmaxVzc'+str(self.flagin)+'.png')
			def tco():
				#========Contact time vs. zc=================
				figure(9)
				grid(False)
				#title('Contact Time (tc/T) vs. Average distance (nm)')
				ylabel('$Contact \;Time \; (t/T)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,tcd,palette[comp_counter])
				savefig('ctimeVzc'+str(self.flagin)+'.png')
			def amp2():
					#========Amp 2 vs. zc=================
				figure(10)
				grid(False)			
				#title('$A_2$($z_c$)')
				ylabel('$A_2 \; (nm)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,amp2d,palette[comp_counter])
				savefig('amp2Vzc'+str(self.flagin)+'.png')
			def pha2():
				#========Phase 2 vs. zc=================
				figure(11)
				grid(False)	
				#title('$\phi_2$($z_c$)')
				ylabel('$\phi_2 \; (deg)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,fase2d,palette[comp_counter])
				savefig('phi2Vzc'+str(self.flagin)+'.png')
			def den2():
				#========Diss Energy 2nd vs. zc=================
				figure(12)
				grid(False)
				#title('Second disp. energy vs. Average distance')
				ylabel('$E_{ts2} * 10^{-20} \; (J)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,ets2d,palette[comp_counter])
				savefig('ets2Vzc'+str(self.flagin)+'.png')
			def dpo2():
				#========Diss Power 1st vs. zc=================
				figure(13)
				grid(False)
				#title('Second disp. power vs. Average distance')
				ylabel('$P_{ts2} * 10^{-20} \; (W)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,pts2d,palette[comp_counter])
				savefig('pts2Vzc'+str(self.flagin)+'.png')
			def vir2():
				#========Virial 2nd vs. zc=================
				figure(14)
				grid(False)		
				#title('Virial second vs. Average distance ')
				ylabel('$V_{ts2} * 10^{-20} \; (J)$',weight='bold',size='x-large')
				xlabel('$z_c\; (nm)$',weight='bold',size='x-large')
				plot(zc,vts2d,palette[comp_counter])
				savefig('vts2Vzc'+str(self.flagin)+'.png')
				#----------------------------------------------------------------
				#===================Time Domain==================================
				#----------------------------------------------------------------
			def ztt():
				#========Instantaneous position vs. time=================
				figure(15)
				grid(False)
				#title("$z_T$($t/T$)")
				ylabel('$z_T \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,zt,palette[comp_counter])
				savefig('totalposTime'+str(self.flagin)+'.png')
			def vtt():
				#========Instantaneous velocity vs. time=================
				figure(16)
				grid(False)
				#title("$v_T$($t/T$)")
				ylabel('$v_T \; (nm/\mu s)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,vt,palette[comp_counter])
				savefig('totalveloTime'+str(self.flagin)+'.png')
			def ftt():
				#========Instantaneous force vs. time=================
				figure(17)
				grid(False)
				#title('$f_T$($t/T$)')
				ylabel('$f_T \; (pN)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,forcet,palette[comp_counter])
				savefig('totalforceTime'+str(self.flagin)+'.png')
			def zt1():
				#========Instantaneous position 1 vs. time=================
				figure(18)
				grid(False)
				#title("$z_1$($t/T$)")
				ylabel('$z_1 \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,z1,palette[comp_counter])
				savefig('FirstposTime'+str(self.flagin)+'.png')
			def zt2():
				#========Instantaneous position 2 vs. time=================
				figure(19)
				grid(False)
				#title("$z_2$($t/T$)")
				ylabel('$z_2 \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,z2,palette[comp_counter])
				savefig('SecondposTime'+str(self.flagin)+'.png')
			def zt3():
				#========Instantaneous position 3 vs. time=================
				figure(20)
				grid(False)
				#title("$z_3$($t/T$)")
				ylabel('$z_3 \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,zt,palette[comp_counter])
				savefig('ThirdposTime'+str(self.flagin)+'.png')
			def zt4():
				#========Instantaneous position 4 vs. time=================
				figure(21)
				grid(False)
				#title("$z_4$($t/T$)")
				ylabel('$z_4 \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,zt,palette[comp_counter])
				savefig('FourthposTime'+str(self.flagin)+'.png')
			def p1a1():
				#========Phase 1 vs. Amp 21=================
				figure(22)
				grid(False)
				#title("$\phi_1$($A_1$)")
				ylabel('$\phi_1 \; (deg)$',weight='bold',size='x-large')
				xlabel('$A_1\; (nm)$',weight='bold',size='x-large')
				plot(amp1d,fase1d,palette[comp_counter])
				savefig('phi1Vamp1'+str(self.flagin)+'.png')
			def p2a2():
				#========Phase 2 vs. Amp 2=================
				figure(23)
				grid(False)
				#title("$\phi_2$($A_2$)")
				ylabel('$\phi_2 \; (deg)$',weight='bold',size='x-large')
				xlabel('$A_2\; (nm)$',weight='bold',size='x-large')
				plot(amp2d,fase2d,palette[comp_counter])
				savefig('phi2Vamp2'+str(self.flagin)+'.png')
			def p2a1():
				#========Phase 2 vs. Amp 1=================
				figure(24)
				grid(False)
				#title("$\phi_2$($A_1$)")
				ylabel('$\phi_2 \; (deg)$',weight='bold',size='x-large')
				xlabel('$A_1\; (nm)$',weight='bold',size='x-large')
				plot(amp1d,fase2d,palette[comp_counter])
				savefig('phi2Vamp1'+str(self.flagin)+'.png')
			def p1a2():
				#========Phase 1 vs. Amp 2=================
				figure(25)
				grid(False)
				#title("$\phi_1$($A_2$)")
				ylabel('$\phi_1 \; (deg)$',weight='bold',size='x-large')
				xlabel('$A_2\; (nm)$',weight='bold',size='x-large')
				plot(amp2d,fase1d,palette[comp_counter])
				savefig('phi1Vamp2'+str(self.flagin)+'.png')
			def famp1graph():
				#========Fourier amplitude 1=================
				figure(26)
				grid(False)
				#title("Amplitude 1")
				ylabel('$A_1\; (logarithmic  units)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],np.log(famp1ex[:freq_basex.size/2]),palette[comp_counter])
				savefig('famp1'+str(self.flagin)+'.png')
			def fpha1graph():
				#========Fourier phase 1=================
				figure(27)
				grid(False)
				#title("Phase 1")
				ylabel('$\phi_1 \; (deg)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],fpha1ex[:freq_basex.size/2],palette[comp_counter])
				savefig('fpha1'+str(self.flagin)+'.png')

			def famp2graph():
				#========Fourier amplitude 2=================
				figure(28)
				grid(False)
				#title("Amplitude 2")
				ylabel('$A_2\; (logarithmic  units)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],np.log(famp2ex[:freq_basex.size/2]),palette[comp_counter])
				savefig('famp2'+str(self.flagin)+'.png')
			def fpha2graph():
				#========Fourier phase 2=================
				figure(29)
				grid(False)
				#title("Phase 2")
				ylabel('$\phi_2 \; (deg)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],fpha2ex[:freq_basex.size/2],palette[comp_counter])
				savefig('fpha2'+str(self.flagin)+'.png')

			def famptgraph():
				#========Fourier amplitude total=================
				figure(30)
				grid(False)
				#title("Amplitude total")
				ylabel('$A \; (logarithmic \, units)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],np.log(famptex[:freq_basex.size/2]),palette[comp_counter])
				savefig('fampT'+str(self.flagin)+'.png')
			def fphatgraph():
				#========Fourier phase total=================
				figure(31)
				grid(False)
				#title("Phase total")
				ylabel('$\phi \; (deg)$',weight='bold',size='x-large')
				xlabel('$Freq  \; (Hz)$',weight='bold',size='x-large')
				plot(freq_basex[:freq_basex.size/2],fphatex[:freq_basex.size/2],palette[comp_counter])
				savefig('fphaT'+str(self.flagin)+'.png')
			def iden():
				#========Fourier phase total=================
				figure(32)
				grid(False)
				#title("Identation")
				ylabel('$Identation \; (nm)$',weight='bold',size='x-large')
				xlabel('$(t/T)$',weight='bold',size='x-large')
				plot(tt,ident,palette[comp_counter])
				savefig('iden'+str(self.flagin)+'.png')




			#========Custom figures===========
			def multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,counter):

				#deleting axix before move them
				def make_patch_spines_invisible(ax):
					ax.set_frame_on(True)
					ax.patch.set_visible(False)
					for sp in ax.spines.itervalues():
					    sp.set_visible(False)


				fig = plt.figure(100+counter)
				grid(False)
				title("Custom graph")
				ax = fig.add_subplot(111)

				fig.subplots_adjust(right=0.8, left = 0.2)

				x_axis = np.arange(len(curves_mixed))
				counter = 0	
				for ii in range(len(curves_mixed)):
					if opt_mixed_x[ii] == 1:
						x_axis = curves_mixed[ii]
						ax.set_xlabel(names_mixed[ii])

				for ii in range(len(curves_mixed)):				
					if opt_mixed_y[ii] == 1:
						counter = counter + 1
						y_axis = curves_mixed[ii]
						if counter == 1:
							lns1 = ax.plot(x_axis, y_axis, 'k-', label = names_mixed[ii])
							ax.set_ylabel(names_mixed[ii])
							ax.yaxis.label.set_color('k')
							ax.tick_params(axis='y', colors = 'k')
						elif counter == 2:
							ax2 = ax.twinx()
							lns2 = ax2.plot(x_axis, y_axis, 'r-', label = names_mixed[ii])
							ax2.set_ylabel(names_mixed[ii])
							ax2.yaxis.label.set_color('r')
							ax2.tick_params(axis='y', colors = 'r')
							ax2.spines['right'].set_color('r')
						elif counter == 3:
							ax3 = ax.twinx()
							lns2 = ax3.plot(x_axis, y_axis, 'g-', label = names_mixed[ii])
							ax3.spines["right"].set_position(("axes", 1.2))
							make_patch_spines_invisible(ax3)
							ax3.spines["right"].set_visible(True)
							ax3.set_ylabel(names_mixed[ii])
							ax3.yaxis.label.set_color('g')
							ax3.tick_params(axis='y', colors = 'g')
							ax3.spines['right'].set_color('g')
						elif counter == 4:
							ax4 = ax.twinx()
							lns2 = ax4.plot(x_axis, y_axis, 'b-', label = names_mixed[ii])
							ax4.spines["right"].set_position(("axes", -0.3))
							make_patch_spines_invisible(ax4)
							ax4.spines["right"].set_visible(True)
							ax4.set_ylabel(names_mixed[ii])
							ax4.yaxis.label.set_color('b')
							ax4.tick_params(axis='y', colors = 'b')
							ax4.spines['right'].set_color('b')



			#Create a Tuple for the command self.glade_ui.get_widget
			curves = [amp1, pha1, den1, dpo1, vir1, def1, dmin1, maf, tco, amp2, pha2, den2, dpo2, vir2, ztt, vtt, ftt, zt1, zt2, zt3, zt4, p1a1, p2a2, p2a1, p1a2,famp1graph,fpha1graph,famp2graph,fpha2graph,famptgraph,fphatgraph,iden]
			opt1=[self.tucom("ampst_cb"),self.tucom("fasest_cb"),self.tucom("etst_cb"),self.tucom("ptst_cb"),self.tucom("virialt_cb"),self.tucom("mdef_cb"),self.tucom("mdis_cb"),self.tucom("maforce_cb"),self.tucom("ctime_cb"),self.tucom("amps2_cb"),self.tucom("fases2_cb"),self.tucom("ets2_cb"),self.tucom("pts2_cb"),self.tucom("virial2_cb"),self.tucom("zt_cb"),self.tucom("vt_cb"),self.tucom("ft_cb"),self.tucom("z1_cb"),self.tucom("z2_cb"),self.tucom("z3_cb"),self.tucom("z4_cb"),self.tucom("phi1a1_cb"),self.tucom("phi2a2_cb"),self.tucom("phi2a1_cb"),self.tucom("phi1a2_cb"),self.tucom("famp1"),self.tucom("fpha1"),self.tucom("famp2"),self.tucom("fpha2"),self.tucom("fampt"),self.tucom("fphat"),self.tucom("it_cb")]
			


			#===========Draw custom figures=================
			if self.flagpm == 1:
				for ii in [1,2,3]:
					end = str(ii)
					opt_mixed_x = [self.tucom("amp_x" + end),self.tucom("phase_x"+end),self.tucom("zc_x"+end),self.tucom("pts_x" + end),self.tucom("ets_x"+end),self.tucom("dmin_x"+end),self.tucom("fmax_x" + end),self.tucom("tcd_x" + end),self.tucom("vts_x" + end)]
					opt_mixed_y = [self.tucom("amp_y" + end),self.tucom("phase_y"+end),self.tucom("zc_y"+end),self.tucom("pts_y" + end),self.tucom("ets_y"+end),self.tucom("dmin_y"+end),self.tucom("fmax_y" + end),self.tucom("tcd_y" + end),self.tucom("vts_y" + end)]
					curves_mixed = [amp1d, fase1d, zc, pts1d, ets1d, dmind,fmaxd,tcd,vts1d]
					names_mixed = ["amp (nm)","$\phi\; (deg)$","zc (nm)","dissipated power","$Energy * 10^{-20} \; (J)$","minimum distance (nm)","maximum force (pN)","contact time (t/T)","virial"]
					if sum(opt_mixed_x) > 0 :	
						multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,ii)

			else:

				opt_mixed_x = [self.tucom("amp_x1"),self.tucom("phase_x1"),self.tucom("zc_x1"),self.tucom("pts_x1"),self.tucom("ets_x1"),self.tucom("dmin_x1"),self.tucom("fmax_x1"),self.tucom("tcd_x1"),self.tucom("vts_x1" ),self.tucom("amp_x2"),self.tucom("phase_x2"),self.tucom("zc_x2"),self.tucom("pts_x2"),self.tucom("ets_x2"),self.tucom("dmin_x2"),self.tucom("fmax_x2"),self.tucom("tcd_x2"),self.tucom("vts_x2" )]
				opt_mixed_y = [self.tucom("amp_y1"),self.tucom("phase_y1"),self.tucom("zc_y1"),self.tucom("pts_y1"),self.tucom("ets_y1"),self.tucom("dmin_y1"),self.tucom("fmax_y1"),self.tucom("tcd_y1"),self.tucom("vts_y1" ),self.tucom("amp_y2"),self.tucom("phase_y2"),self.tucom("zc_y2"),self.tucom("pts_y2"),self.tucom("ets_y2"),self.tucom("dmin_y2"),self.tucom("fmax_y2"),self.tucom("tcd_y2"),self.tucom("vts_y2" )]
				curves_mixed = [amp1d, fase1d, zc, pts1d, ets1d, dmind, fmaxd, tcd, vts1d, amp2d, fase2d, zc, pts2d, ets2d, dmind, fmaxd, tcd, vts2d]
				names_mixed = ["amp 1  (nm)","$\phi_1 \; (deg)$","zc (nm)","$Power_1 * 10^{-16} \; (W)$","$Energy_1 * 10^{-20} \; (J)$","dmin (nm)","maximum force (pN)","contact time (t/T)","$Virial_1 * 10^{-20} \; (J)$","amp 2 (nm)","$\phi_2 \; (deg)$"," (nm)","$Power_2 * 10^{-16} \; (W)$","$Energy_2 * 10^{-20} \; (J)$","dmin (nm)","maximum force (pN)","contact time (t/T)","$Virial_2 * 10^{-20} \; (J)$"]
				if sum(opt_mixed_x) > 0 :	
					multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,4)


			#============Bias figure=====================
			if self.tucom("bias") == 1:
				if self.flagpm == 1:
					opt_mixed_x = [1,0,0,0,0]
					opt_mixed_y = [0,self.tucom("zt_cb"),self.tucom("vt_cb"),self.tucom("ft_cb"),self.tucom("it_cb")]
					curves_mixed = [tt,zt, vt,forcet,ident]
					names_mixed = ["time (t/T)","z (nm)","$v_T \; (nm/\mu s)$","forc (pN)","identation (nm)"]
					if sum(opt_mixed_x) > 0 :	
						multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,2)
				else:
					opt_mixed_x = [1,0,0,0,0,0,0]
					opt_mixed_y = [0,self.tucom("zt_cb"),self.tucom("vt_cb"),self.tucom("ft_cb"),self.tucom("it_cb"),self.tucom("z1_cb"),self.tucom("z2_cb")]
					curves_mixed = [tt,zt, vt,forcet,ident,z1,z2]
					names_mixed = ["time (t/T)","z (nm)","$v_T \; (nm/\mu s)$","force (pN)","identation (nm)","z1 (nm)","z2 (nm)"]
					if sum(opt_mixed_x) > 0 :	
						multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,3)








			#||||||||||||||||||||||||||||||||||||||||||
			#=========Big switch ready to use==============
			print "This is www.forceforfuture.es"
			for ii in range(len(opt1)):
				if opt1[ii]==1:
					curves[ii]()	#calling ccurves.py and saving plots
				else:
					print "You did not visualize", ii+1	#-----Message
			#self.flagin+=1
			if self.flagfile == 0:
				os.chdir(self.piat)
		else:
			print "You have nothing to visualize"	#-----Message
	show()		# intented for plotting the bounch of curves...



#--------------------------------End of run button----------------------------------------------------------------------------------------------------
 

#--------------------Load button-------------------------------------------------

    def on_load_clicked(self,widget):
        	self.uploadingd.show()


    def on_compare_sp_clicked(self,widget):
    	if (self.tucom("compare_sp")):
    		self.glade_ui.get_widget("zdomload_sp1").set_sensitive(True)
    		self.glade_ui.get_widget("tdomload_sp1").set_sensitive(True)
    		self.compareflag = 1
    	else:
    		self.glade_ui.get_widget("zdomload_sp1").set_sensitive(False)
    		self.glade_ui.get_widget("tdomload_sp1").set_sensitive(False)
    		self.compareflag = 0


    def on_notfile_clicked(self, widget):
        	self.uploadingd.hide()


    def on_okfile_clicked(self, widget):
        	self.tdomload=self.glade_ui.get_widget("tdomload_sp").get_filename()
        	print(self.tdomload)
        	self.zdomload=self.glade_ui.get_widget("zdomload_sp").get_filename()
        	print(self.zdomload)

        	self.tdomload1=self.glade_ui.get_widget("tdomload_sp1").get_filename()
        	print(self.tdomload1)
        	self.zdomload1=self.glade_ui.get_widget("zdomload_sp1").get_filename()
        	print(self.zdomload1)
        	self.flagfile=1

        	self.uploadingd.hide()


#--------------Show value table button-------------------------------------------------

    def on_show_value_clicked(self,widget):
    	self.value_table.show()

    def on_value_table_exit_clicked(self,widget):
    	self.value_table.hide()

#---------------------Forces buttons----------------------------------------------------------

# The main values are referenced in:
# Sample radious = sr_sp
# Sample Young = saym_sp
# Sample mu = musam_sp
# Sample relative epsilon = epsilon_sp
# Sample sigma = sigmas_sp
# Tip sigma = sigmat_sp
# Sample Deby lenght = landadeb_sp
# Sample Hamaker = ham_sp
# Sample atomic radious = a00_sp
# Sample viscosity = mvis_sp

# Tip radius radious = tr_sp
# Tip Young's modulus = ymt_sp
# Tip poisson = mutip_sp



    def on_vdw_check_clicked(self,widget):
    	if (self.tucom("vdw_check")):
    		self.glade_ui.get_widget("dlvo_check").set_state(False)	

        	local = self.glade_ui.get_widget("ham_sp").get_value()  
        	self.glade_ui.get_widget("ham_vdw").set_value(local)

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_vdw").set_value(local)


        	self.vdw_par.show()

    def on_vdw_exit_clicked(self,widget):
    	self.vdw_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("ham_vdw").get_value()  
    	self.glade_ui.get_widget("ham_sp").set_value(local)

    	local = self.glade_ui.get_widget("tr_vdw").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)




    def on_dlvo_check_clicked(self,widget):
    	if (self.tucom("dlvo_check")):
    		self.glade_ui.get_widget("vdw_check").set_state(False)	

        	local = self.glade_ui.get_widget("epsilon_sp").get_value()  
        	self.glade_ui.get_widget("epsilon_dlvo").set_value(local)

        	local = self.glade_ui.get_widget("landadeb_sp").get_value()  
        	self.glade_ui.get_widget("landadeb_dlvo").set_value(local)

        	local = self.glade_ui.get_widget("sigmas_sp").get_value()  
        	self.glade_ui.get_widget("sigmas_dlvo").set_value(local)

        	local = self.glade_ui.get_widget("sigmat_sp").get_value()  
        	self.glade_ui.get_widget("sigmat_dlvo").set_value(local) 

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_dlvo").set_value(local)         	


        	self.dlvo_par.show()

    def on_dlvo_exit_clicked(self,widget):
    	self.dlvo_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("epsilon_dlvo").get_value()  
    	self.glade_ui.get_widget("epsilon_sp").set_value(local)

    	local = self.glade_ui.get_widget("landadeb_dlvo").get_value()  
    	self.glade_ui.get_widget("landadeb_sp").set_value(local)

    	local = self.glade_ui.get_widget("sigmas_dlvo").get_value()  
    	self.glade_ui.get_widget("sigmas_sp").set_value(local)

    	local = self.glade_ui.get_widget("sigmat_dlvo").get_value()  
    	self.glade_ui.get_widget("sigmat_sp").set_value(local)  

    	local = self.glade_ui.get_widget("tr_dlvo").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)  




    def on_dmt_check_clicked(self,widget):
    	if (self.tucom("dmt_check")):
    		self.glade_ui.get_widget("tatara_check").set_state(False)
    		self.glade_ui.get_widget("hertz_check").set_state(False)
    		self.glade_ui.get_widget("jkr_check").set_state(False)

        	local = self.glade_ui.get_widget("saym_sp").get_value()  
        	self.glade_ui.get_widget("saym_dmt").set_value(local)

        	local = self.glade_ui.get_widget("a00_sp").get_value()  
        	self.glade_ui.get_widget("a00_dmt").set_value(local)

        	local = self.glade_ui.get_widget("ham_sp").get_value()  
        	self.glade_ui.get_widget("ham_dmt").set_value(local)

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_dmt").set_value(local)  

        	local = self.glade_ui.get_widget("ymt_sp").get_value()  
        	self.glade_ui.get_widget("ymt_dmt").set_value(local) 

        	local = self.glade_ui.get_widget("mutip_sp").get_value()  
        	self.glade_ui.get_widget("mutip_dmt").set_value(local)

        	local = self.glade_ui.get_widget("musam_sp").get_value()  
        	self.glade_ui.get_widget("musam_dmt").set_value(local)



        	self.dmt_par.show()

    def on_dmt_exit_clicked(self,widget):
    	self.dmt_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("saym_dmt").get_value()  
    	self.glade_ui.get_widget("saym_sp").set_value(local)

    	local = self.glade_ui.get_widget("a00_dmt").get_value()  
    	self.glade_ui.get_widget("a00_sp").set_value(local)

    	local = self.glade_ui.get_widget("ham_dmt").get_value()  
    	self.glade_ui.get_widget("ham_sp").set_value(local)

    	local = self.glade_ui.get_widget("tr_dmt").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)  

    	local = self.glade_ui.get_widget("ymt_dmt").get_value()  
    	self.glade_ui.get_widget("ymt_sp").set_value(local) 

    	local = self.glade_ui.get_widget("mutip_dmt").get_value()  
    	self.glade_ui.get_widget("mutip_sp").set_value(local)

    	local = self.glade_ui.get_widget("musam_dmt").get_value()  
    	self.glade_ui.get_widget("musam_sp").set_value(local)



    def on_hertz_check_clicked(self,widget):
    	if (self.tucom("hertz_check")):
    		self.glade_ui.get_widget("tatara_check").set_state(False)
    		self.glade_ui.get_widget("dmt_check").set_state(False)
    		self.glade_ui.get_widget("jkr_check").set_state(False)

        	local = self.glade_ui.get_widget("saym_sp").get_value()  
        	self.glade_ui.get_widget("saym_hertz").set_value(local)

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_hertz").set_value(local)

        	local = self.glade_ui.get_widget("ymt_sp").get_value()  
        	self.glade_ui.get_widget("ymt_hertz").set_value(local)

        	local = self.glade_ui.get_widget("mutip_sp").get_value()  
        	self.glade_ui.get_widget("mutip_hertz").set_value(local)

        	local = self.glade_ui.get_widget("musam_sp").get_value()  
        	self.glade_ui.get_widget("musam_hertz").set_value(local) 

        	self.hertz_par.show()

    def on_hertz_exit_clicked(self,widget):

    	self.hertz_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("saym_hertz").get_value()  
    	self.glade_ui.get_widget("saym_sp").set_value(local)

    	local = self.glade_ui.get_widget("tr_hertz").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)

    	local = self.glade_ui.get_widget("ymt_hertz").get_value()  
    	self.glade_ui.get_widget("ymt_sp").set_value(local)

    	local = self.glade_ui.get_widget("mutip_hertz").get_value()  
    	self.glade_ui.get_widget("mutip_sp").set_value(local)

    	local = self.glade_ui.get_widget("musam_hertz").get_value()  
    	self.glade_ui.get_widget("musam_sp").set_value(local)


    def on_jkr_check_clicked(self,widget):
    	if (self.tucom("jkr_check")):
    		self.glade_ui.get_widget("tatara_check").set_state(False)
    		self.glade_ui.get_widget("hertz_check").set_state(False)
    		self.glade_ui.get_widget("dmt_check").set_state(False)
        	self.jkr_par.show()

    def on_jkr_exit_clicked(self,widget):
    	self.jkr_par.hide()




    def on_tatara_check_clicked(self,widget):
    	if (self.tucom("tatara_check")):
    		self.glade_ui.get_widget("dmt_check").set_state(False)
    		self.glade_ui.get_widget("hertz_check").set_state(False)
    		self.glade_ui.get_widget("jkr_check").set_state(False)

    		self.glade_ui.get_widget("viscosity_check").set_state(False)
    		self.glade_ui.get_widget("viscosity_check").set_sensitive(False)     		

        	local = self.glade_ui.get_widget("ymt_sp").get_value()  
        	self.glade_ui.get_widget("ymt_tat").set_value(local)

        	local = self.glade_ui.get_widget("saym_sp").get_value()  
        	self.glade_ui.get_widget("saym_tat").set_value(local)

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_tat").set_value(local)

        	local = self.glade_ui.get_widget("sr_sp").get_value()  
        	self.glade_ui.get_widget("sr_tat").set_value(local)  

        	local = self.glade_ui.get_widget("mutip_sp").get_value()  
        	self.glade_ui.get_widget("mutip_tat").set_value(local)

        	local = self.glade_ui.get_widget("musam_sp").get_value()  
        	self.glade_ui.get_widget("musam_tat").set_value(local) 


        	self.tatara_par.show()

        else:
    		self.glade_ui.get_widget("viscosity_check").set_sensitive(True)   


    def on_tatara_exit_clicked(self,widget):
    	self.tatara_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("ymt_tat").get_value()  
    	self.glade_ui.get_widget("ymt_sp").set_value(local)

    	local = self.glade_ui.get_widget("saym_tat").get_value()  
    	self.glade_ui.get_widget("saym_sp").set_value(local)

    	local = self.glade_ui.get_widget("tr_tat").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)

    	local = self.glade_ui.get_widget("sr_tat").get_value()  
    	self.glade_ui.get_widget("sr_sp").set_value(local)  

    	local = self.glade_ui.get_widget("mutip_tat").get_value()  
    	self.glade_ui.get_widget("mutip_sp").set_value(local)

    	local = self.glade_ui.get_widget("musam_tat").get_value()  
    	self.glade_ui.get_widget("musam_sp").set_value(local) 


    def on_viscosity_check_clicked(self,widget):
    	if (self.tucom("viscosity_check")):

        	local = self.glade_ui.get_widget("tr_sp").get_value()  
        	self.glade_ui.get_widget("tr_vis").set_value(local)

        	local = self.glade_ui.get_widget("mvis_sp").get_value()  
        	self.glade_ui.get_widget("mvis_vis").set_value(local)

        	self.glade_ui.get_widget("pts_x1").set_sensitive(True)
        	self.glade_ui.get_widget("pts_y1").set_sensitive(True)
        	self.glade_ui.get_widget("pts_x2").set_sensitive(True)
        	self.glade_ui.get_widget("pts_y2").set_sensitive(True)
        	self.glade_ui.get_widget("pts_x3").set_sensitive(True)
        	self.glade_ui.get_widget("pts_y3").set_sensitive(True)
        	self.glade_ui.get_widget("ets_x1").set_sensitive(True)
        	self.glade_ui.get_widget("ets_y1").set_sensitive(True)
        	self.glade_ui.get_widget("ets_x2").set_sensitive(True)
        	self.glade_ui.get_widget("ets_y2").set_sensitive(True)
        	self.glade_ui.get_widget("ets_x3").set_sensitive(True)
        	self.glade_ui.get_widget("ets_y3").set_sensitive(True)   



        	self.glade_ui.get_widget("etst_cb").set_sensitive(True)
        	self.glade_ui.get_widget("ptst_cb").set_sensitive(True)
        	if self.tucom("eb_cb"):
	        	self.glade_ui.get_widget("ets2_cb").set_sensitive(True)
	        	self.glade_ui.get_widget("pts2_cb").set_sensitive(True)       		

        	self.viscosity_par.show()

        else:

        	self.glade_ui.get_widget("pts_x1").set_sensitive(False)
        	self.glade_ui.get_widget("pts_y1").set_sensitive(False)
        	self.glade_ui.get_widget("pts_x2").set_sensitive(False)
        	self.glade_ui.get_widget("pts_y2").set_sensitive(False)
        	self.glade_ui.get_widget("pts_x3").set_sensitive(False)
        	self.glade_ui.get_widget("pts_y3").set_sensitive(False)
        	self.glade_ui.get_widget("ets_x1").set_sensitive(False)
        	self.glade_ui.get_widget("ets_y1").set_sensitive(False)
        	self.glade_ui.get_widget("ets_x2").set_sensitive(False)
        	self.glade_ui.get_widget("ets_y2").set_sensitive(False)
        	self.glade_ui.get_widget("ets_x3").set_sensitive(False)
        	self.glade_ui.get_widget("ets_y3").set_sensitive(False) 


        	self.glade_ui.get_widget("etst_cb").set_sensitive(False)
        	self.glade_ui.get_widget("ptst_cb").set_sensitive(False)
        	if self.tucom("eb_cb"):
	        	self.glade_ui.get_widget("ets2_cb").set_sensitive(False)
	        	self.glade_ui.get_widget("pts2_cb").set_sensitive(False)           	

    def on_viscosity_exit_clicked(self,widget):
    	self.viscosity_par.hide()
    	while gtk.events_pending():
	   		gtk.main_iteration(False)   #Refresh the GUI
    	local = self.glade_ui.get_widget("tr_vis").get_value()  
    	self.glade_ui.get_widget("tr_sp").set_value(local)

    	local = self.glade_ui.get_widget("mvis_vis").get_value()  
    	self.glade_ui.get_widget("mvis_sp").set_value(local)





#------------------------------Credits button-----------------------------------------------------   
    def on_about_clicked(self, widget):
        self.credits.show()
        self.window1.hide()

    def on_exit_credits_clicked(self, widget):
        self.credits.hide()
        self.window1.show()
       
#--------------------------------------------------------------------------------------------------


#------------------------------Abort button-----------------------------------------------------   
    def on_abort_clicked(self, widget):
        self.abort = 1
       
#--------------------------------------------------------------------------------------------------


    def on_save_ex_clicked(self, widget):
        print "Your simulation data and curves have been saved please register to www.dforce.net"
        self.savingd.show()

	#||||||||||||||||||||||||||||||||||||||||||
	#=========OK and NOK options==============
	def on_oksave_clicked(self, widget):		#not working
		self.pnamed=self.glade_ui.get_widget("pname_tin").get_text()
		self.pownerd=self.glade_ui.get_widget("powner_tin").get_text()
		self.flagpnam=1
		print "The project:",self.pnamed,"owned by:",self.pownerd
	def on_notsave_clicked(self, widget):
		print "No name..."
#------------------------------------------------------------------------------------------------------------------------------------
    def on_inputs_ex_clicked(self, widget):
	print "Expected to save the file in the directory name as the main.py is running"
        self.uploadingd.show()
	#self.infname=self.glade_ui.get_widget("inpfname_tin").get_text()
	# Pending on Window menu for saving a file
	#||||||||||||||||||||||||||||||||||||||||||
	#=========OK and NOK options==============

	def on_stopulfile_clicked(self, widget):
		self.flagfile=0
		print "Back to gui parameter inputs"

	#||||||||||||||||||||||||||||||||||||||||||||||||||
	#========= CheckButton PM ==================================

	#|||||||||   New Checkbutton  CALLBACK |||||||||||||||||||||||||
#    def on_pm_cb_clicked(self,widget):
#	self.glade_ui.get_widget("eb_cb").set_sensitive(False)
		
#    def on_pm_cb_leave(self,widget):
#	self.glade_ui.get_widget("eb_cb").set_sensitive(True)
#    def on_pm_cb_released(self,widget):
#	self.glade_ui.get_widget("eb_cb").set_sensitive(True)

	#|||||||||   New Checkbutton  CALLBACK |||||||||||||||||||||||||
#    def on_eb_cb_clicked(self,widget):
#	self.glade_ui.get_widget("pm_cb").set_sensitive(False)
		
#    def on_eb_cb_leave(self,widget):
#	self.glade_ui.get_widget("pm_cb").set_sensitive(True)
#    def on_eb_cb_released(self,widget):
#	self.glade_ui.get_widget("pm_cb").set_sensitive(True)
	#|||||||||   New Checkbutton  CALLBACK |||||||||||||||||||||||||
#    def on_single_cb_clicked(self,widget):
#	self.glade_ui.get_widget("bmafm_cb").set_sensitive(False)
#	self.glade_ui.get_widget("eb_cb").set_sensitive(False)

		
#    def on_single_cb_leave(self,widget):
#	self.glade_ui.get_widget("bmafm_cb").set_sensitive(True)
#    def on_single_cb_released(self,widget):
#	self.glade_ui.get_widget("bmafm_cb").set_sensitive(True)
	#||||||||| New Checkbutton  CALLBACK Crossed butts |||||||||||||||||||

    def on_bmafm_cb_clicked(self,widget):
    	if self.tucom("bmafm_cb"):

    		self.glade_ui.get_widget("single_cb").set_sensitive(False)
    		self.glade_ui.get_widget("pm_cb").set_sensitive(False)
    		self.glade_ui.get_widget("single_cb").set_state(False)
    		self.glade_ui.get_widget("pm_cb").set_state(False)

    		self.glade_ui.get_widget("kc2_sp").set_sensitive(True)
    		self.glade_ui.get_widget("f02_sp").set_sensitive(True)
    		self.glade_ui.get_widget("q2_sp").set_sensitive(True)
    		self.glade_ui.get_widget("fcorrection_q2").set_sensitive(True)
    		self.glade_ui.get_widget("a02_sp").set_sensitive(True)
    		self.glade_ui.get_widget("fd2_sp").set_sensitive(True)    
    		self.glade_ui.get_widget("copy_fre2").set_sensitive(True)

    		self.glade_ui.get_widget("amps2_cb").set_sensitive(True)
    		self.glade_ui.get_widget("fases2_cb").set_sensitive(True)
    		self.glade_ui.get_widget("virial2_cb").set_sensitive(True)
    		self.glade_ui.get_widget("z1_cb").set_sensitive(True)     		
    		self.glade_ui.get_widget("z2_cb").set_sensitive(True)    
    		self.glade_ui.get_widget("phi2a2_cb").set_sensitive(True)
    		self.glade_ui.get_widget("phi2a1_cb").set_sensitive(True)
    		self.glade_ui.get_widget("phi1a2_cb").set_sensitive(True)    
    		self.glade_ui.get_widget("famp2").set_sensitive(True)
    		self.glade_ui.get_widget("fpha2").set_sensitive(True)
    		self.glade_ui.get_widget("famp1").set_sensitive(True)
    		self.glade_ui.get_widget("fpha1").set_sensitive(True)
    	else:
    		self.glade_ui.get_widget("single_cb").set_sensitive(True)
    		self.glade_ui.get_widget("pm_cb").set_sensitive(True)

    		self.glade_ui.get_widget("kc2_sp").set_sensitive(False)
    		self.glade_ui.get_widget("f02_sp").set_sensitive(False)
    		self.glade_ui.get_widget("q2_sp").set_sensitive(False)
    		self.glade_ui.get_widget("fcorrection_q2").set_sensitive(False)
    		self.glade_ui.get_widget("a02_sp").set_sensitive(False)
    		self.glade_ui.get_widget("fd2_sp").set_sensitive(False)    
    		self.glade_ui.get_widget("copy_fre2").set_sensitive(False)

    		self.glade_ui.get_widget("amps2_cb").set_sensitive(False)
    		self.glade_ui.get_widget("fases2_cb").set_sensitive(False)
    		self.glade_ui.get_widget("virial2_cb").set_sensitive(False)
    		self.glade_ui.get_widget("z1_cb").set_sensitive(False) 
    		self.glade_ui.get_widget("z2_cb").set_sensitive(False)    
    		self.glade_ui.get_widget("phi2a2_cb").set_sensitive(False)
    		self.glade_ui.get_widget("phi2a1_cb").set_sensitive(False)
    		self.glade_ui.get_widget("phi1a2_cb").set_sensitive(False)    
    		self.glade_ui.get_widget("famp2").set_sensitive(False)
    		self.glade_ui.get_widget("fpha2").set_sensitive(False)
    		self.glade_ui.get_widget("famp1").set_sensitive(False)
    		self.glade_ui.get_widget("fpha1").set_sensitive(False)     		  
    		self.glade_ui.get_widget("ets2_cb").set_sensitive(False)
    		self.glade_ui.get_widget("pts2_cb").set_sensitive(False)  



#    def on_bmafm_cb_leave(self,widget):
#	self.glade_ui.get_widget("single_cb").set_sensitive(True)
#	self.glade_ui.get_widget("pm_cb").set_sensitive(True)
#    def on_bmafm_cb_released(self,widget):
#	self.glade_ui.get_widget("single_cb").set_sensitive(True)
#	self.glade_ui.get_widget("pm_cb").set_sensitive(True)

	#|||||||||   New Checkbutton  CALLBACK Crossed buttons||||||||||||||||||||



#------------------------------------------------------------------------------------------------------------------------------------
#============================run main=========================================
#------------------------------------------------------------------------------------------------------------------------------------	
if __name__ == '__main__':
    dolly = dforce()
    gtk.main()
