
# coding: utf-8

# In[1]:

from tkinter import *
from tkFileDialog   import askopenfilename
from tkFileDialog import asksaveasfilename
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.optimize import nnls
from scipy.linalg import toeplitz
from scipy.special import kv, iv, gamma
from scipy.integrate import quad
from cmath import *
from math import ceil,floor
import mpmath as mp
mp.dps = 25; mp.pretty = True
global F
global R
F = 96485.0
R = 8.314
def cot(phi):
    return 1.0/tan(phi)
def csc(phi):
    return 1.0/sin(phi)
def coth(x):
    return 1/tanh(x)

#===========================================================================================
#===========================================================================================
#===========================================================================================
#IMPORT ALL WRITTEN MODULES
#===========================================================================================
#===========================================================================================
#===========================================================================================

from Koutecky_Levich             import *
from Tafel                       import *
from Randles_Rev                 import *
from Randles_Irr                 import *
from Cottrell                    import *
from Impedance                   import *
from DRT_Transformation          import *
from AC_Cyclic_Voltammetry       import *
from Cyclic_Voltammetry          import *
from RIS_Cyclic_Voltammetry      import *
from ACCV_Superimposer           import *
from CV_Superimposer             import *
from RISCV_Superimposer          import *
from LSA_Cyclic_Voltammetry      import *
from LSACV_Superimposer          import *
from STC_Cyclic_Voltammetry      import *
from STCCV_Superimposer          import *
from NST_Chronoamperometry       import *
from NSTCA_Superimposer          import *


#==========================================================================================
#Define Data selection-functions
#==========================================================================================

def CV_Data():
    Get_CV_Data()
def CA_Data():
    Get_NSTCA_Data() 
def ACCV_Data():
    Get_ACCV_Data() 
def LSACV_Data():
    Get_LSACV_Data()
def STCCV_Data():
    Get_STCCV_Data()
def ImpFit_Data():
    Get_ImpFit_Data()
def ImpDRT_Data():
    Get_ImpDRT_Data()
    
def NotNow():
    Fenster = Toplevel()                                                         
    Fenster.title("Function will follow soon")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "This function is not available\nright now. Tim is doing his best\nto make it soon.")

    
    
#===========================================================================================
#===========================================================================================
#===========================================================================================
#Build main Window
#===========================================================================================
#===========================================================================================
#===========================================================================================

root = Tk()
root.title("Polarographica")                         
root.geometry("800x500")

menu = Menu(root)                                               
root.config(menu=menu)

#=================================================================================
#=================================================================================

filemenu = Menu(menu)                                           
menu.add_cascade(label="File", menu=filemenu)           
filemenu.add_command(label="Open CV/LSV for fitting", command=CV_Data)
filemenu.add_command(label="Open AC-CV/AC-LSV for FFT-Analysis", command=ACCV_Data)
filemenu.add_command(label="Open n-Step CA for fitting", command=CA_Data)
filemenu.add_command(label="Open LSA-CV for fitting", command=LSACV_Data)
filemenu.add_command(label="Open STC-CV for fitting", command=STCCV_Data)
filemenu.add_command(label="Open Impedance/Nyquist for fitting", command=ImpFit_Data)
filemenu.add_command(label="Open Impedance/Nyquist for DRT analysis", command=ImpDRT_Data)

#=================================================================================
#=================================================================================
#=================================================================================
#Classical Evaluation techniques
#=================================================================================
#=================================================================================
#=================================================================================

SemiInfDiff   = Menu(menu)
menu.add_cascade(label="Classical Evaluations", menu=SemiInfDiff)
SemiInfDiff.add_command(label="Randles Sevcik Reversible", command=Get_RSREV_Data)
SemiInfDiff.add_command(label="Randles Sevcik Irreversible", command=Get_RSIRR_Data)
SemiInfDiff.add_command(label="Koutecky-Levich Analysis", command=Get_KL_Data)
SemiInfDiff.add_command(label="Tafel and Shape-Analysis", command=Get_Tafel_Data)
SemiInfDiff.add_command(label="Cottrell-Analysis", command=Get_Cottrell_Data)

#=================================================================================
#=================================================================================
#=================================================================================
#Voltamperometric Simulation-techniques
#=================================================================================
#=================================================================================
#=================================================================================

SimVoltamperometric = Menu(menu)
menu.add_cascade(label="Voltamperometric Simulators", menu=SimVoltamperometric)

#----------------------------
#Cyclic Voltammetry simulators
#----------------------------
CVSimulators   = Menu(SimVoltamperometric)
CVSimulators.add_command(label="Planar Semi-Infinite Diffusion CV Simulator",               command=Semi_Inf_Planar)
CVSimulators.add_command(label="Planar Finite Reflective Diffusion CV Simulator",           command=Finit_Planar)
CVSimulators.add_command(label="Planar Finite Transmissive Diffusion CV Simulator",         command=Finit_Planar_Trans)
CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Simulator", command=Semi_Inf_Zyl_Ext)
CVSimulators.add_command(label="External Cylindrical Finite Diffusion CV Simulator",        command=Finit_Zyl_Ext)
CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion CV Simulator",        command=Finit_Zyl_Int)
CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion CV Simulator",   command=Semi_Inf_Sphere_Ext)
CVSimulators.add_command(label="Internal Spherical Finite Diffusion CV Simulator",          command=Finit_Sphere_Int)
CVSimulators.add_command(label="--------------------------------------------------------------------")
CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion CV Simulator",               command=Statistical_Finit_Planar)
CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Ext)
CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Int)
CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Simulator",   command=Statistical_Finit_Sphere_Int)
CVSimulators.add_command(label="--------------------------------------------------------------------")
CVSimulators.add_command(label = 'Superimpose different CVs', command=CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Cyclic Voltammetry', menu = CVSimulators)
#----------------------------
#AC-Cyclic Voltammetry Simulators
#----------------------------
ACCVSimulators   = Menu(SimVoltamperometric)
ACCVSimulators.add_command(label="Planar Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Planar)
ACCVSimulators.add_command(label="Planar Finite Reflective Diffusion AC-CV Simulator", command=AC_Finit_Planar)
ACCVSimulators.add_command(label="Planar Finite Transmissive Diffusion AC-CV Simulator", command=AC_Finit_Transm_Planar)
ACCVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Zyl_Ext)
ACCVSimulators.add_command(label="External Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Finit_Zyl_Ext)
ACCVSimulators.add_command(label="Internal Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Finit_Zyl_Int)
ACCVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Sphere_Ext)
ACCVSimulators.add_command(label="Internal Spherical Finite Diffusion AC-CV Simulator", command=AC_Finit_Sphere_Int)
ACCVSimulators.add_command(label="--------------------------------------------------------------------")
ACCVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Planar)
ACCVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Zyl_Ext)
ACCVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Zyl_Int)
ACCVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Sphere_Int)
ACCVSimulators.add_command(label="--------------------------------------------------------------------")
ACCVSimulators.add_command(label = 'Superimpose different AC-CVs', command=ACCV_SIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'AC-Cyclic Voltammetry', menu = ACCVSimulators)
#----------------------------
#Large Sine Amplitude Cyclic Voltammetry simulators
#----------------------------
LSA_CVSimulators   = Menu(SimVoltamperometric)
LSA_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion LSA-CV Simulator",               command=LSA_Semi_Inf_Planar)
LSA_CVSimulators.add_command(label="Planar Finite Reflective Diffusion LSA-CV Simulator",           command=LSA_Finit_Planar)
LSA_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion LSA-CV Simulator",         command=LSA_Finit_Planar_Trans)
LSA_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion LSA-CV Simulator", command=LSA_Semi_Inf_Zyl_Ext)
LSA_CVSimulators.add_command(label="External Cylindrical Finite Diffusion LSA-CV Simulator",        command=LSA_Finit_Zyl_Ext)
LSA_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion LSA-CV Simulator",        command=LSA_Finit_Zyl_Int)
LSA_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion LSA-CV Simulator",   command=LSA_Semi_Inf_Sphere_Ext)
LSA_CVSimulators.add_command(label="Internal Spherical Finite Diffusion LSA-CV Simulator",          command=LSA_Finit_Sphere_Int)
LSA_CVSimulators.add_command(label="--------------------------------------------------------------------")
LSA_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion LSA-CV Simulator",               command=LSA_Statistical_Finit_Planar)
LSA_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion LSA-CV Simulator", command=LSA_Statistical_Finit_Zyl_Ext)
LSA_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion LSA-CV Simulator", command=LSA_Statistical_Finit_Zyl_Int)
LSA_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion LSA-CV Simulator",   command=LSA_Statistical_Finit_Sphere_Int)
LSA_CVSimulators.add_command(label="--------------------------------------------------------------------")
LSA_CVSimulators.add_command(label = 'Superimpose different LSA-CVs', command=LSA_CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Large Sine Amplitude Voltammetry', menu = LSA_CVSimulators)
#----------------------------
#n-Step Chronoamperometry Simulators
#----------------------------
CASimulators   = Menu(SimVoltamperometric)
CASimulators.add_command(label="Planar Semi-Infinite Diffusion n-Step CA Simulator",               command=NST_Semi_Inf_Planar)
CASimulators.add_command(label="Planar Finite Reflective Diffusion n-Step CA Simulator",           command=NST_Finit_Planar)
CASimulators.add_command(label="Planar Finite Transmissive Diffusion n-Step CA Simulator",         command=NST_Finit_Planar_Trans)
CASimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion n-Step CA Simulator", command=NST_Semi_Inf_Zyl_Ext)
CASimulators.add_command(label="External Cylindrical Finite Diffusion n-Step CA Simulator",        command=NST_Finit_Zyl_Ext)
CASimulators.add_command(label="Internal Cylindrical Finite Diffusion n-Step CA Simulator",        command=NST_Semi_Inf_Sphere_Ext)
CASimulators.add_command(label="External Spherical Semi-Infinite Diffusion n-Step CA Simulator",   command=NST_Semi_Inf_Sphere_Ext)
CASimulators.add_command(label="Internal Spherical Finite Diffusion n-Step CA Simulator",          command=NST_Finit_Sphere_Int)
CASimulators.add_command(label="--------------------------------------------------------------------")
CASimulators.add_command(label="Statistically Weighted Planar Finite Diffusion n-Step CA Simulator",                command=NST_Statistical_Finit_Planar)
CASimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion n-Step CA Simulator",  command=NST_Statistical_Finit_Zyl_Ext)
CASimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion n-Step CA Simulator",  command=NST_Statistical_Finit_Zyl_Int)
CASimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion n-Step CA Simulator",    command=NST_Statistical_Finit_Sphere_Int)
CASimulators.add_command(label="--------------------------------------------------------------------")
CASimulators.add_command(label = 'Superimpose different LSA-CVs', command=NST_CASIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'n-Step Chronoamperometry', menu = CASimulators)
#----------------------------
#Cyclic Staicase Voltammetry simulators
#----------------------------
STC_CVSimulators   = Menu(SimVoltamperometric)
STC_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion Staircase-CV Simulator",               command=STC_Semi_Inf_Planar)
STC_CVSimulators.add_command(label="Planar Finite Reflective Diffusion Staircase-CV Simulator",           command=STC_Finit_Planar)
STC_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion Staircase-CV Simulator",         command=STC_Finit_Planar_Trans)
STC_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion Staircase-CV Simulator", command=STC_Semi_Inf_Zyl_Ext)
STC_CVSimulators.add_command(label="External Cylindrical Finite Diffusion Staircase-CV Simulator",        command=STC_Finit_Zyl_Ext)
STC_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion Staircase-CV Simulator",        command=STC_Finit_Zyl_Int)
STC_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion Staircase-CV Simulator",   command=STC_Semi_Inf_Sphere_Ext)
STC_CVSimulators.add_command(label="Internal Spherical Finite Diffusion Staircase-CV Simulator",          command=STC_Finit_Sphere_Int)
STC_CVSimulators.add_command(label="--------------------------------------------------------------------")
STC_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion Staircase-CV Simulator",               command=STC_Statistical_Finit_Planar)
STC_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion Staircase-CV Simulator", command=STC_Statistical_Finit_Zyl_Ext)
STC_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion Staircase-CV Simulator", command=STC_Statistical_Finit_Zyl_Int)
STC_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion Staircase-CV Simulator",   command=STC_Statistical_Finit_Sphere_Int)
STC_CVSimulators.add_command(label="--------------------------------------------------------------------")
STC_CVSimulators.add_command(label = 'Superimpose different Staircase-CVs', command=STC_CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Staircase Cyclic Voltammetry', menu = STC_CVSimulators)


#----------------------------
#RIS_Cyclic Voltammetry simulators
#----------------------------
RIS_CVSimulators   = Menu(SimVoltamperometric)
RIS_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion RIS-CV Simulator",               command=RIS_Semi_Inf_Planar)
RIS_CVSimulators.add_command(label="Planar Finite Reflective Diffusion RIS-CV Simulator",           command=RIS_Finit_Planar)
RIS_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion RIS-CV Simulator",         command=RIS_Finit_Planar_Trans)
RIS_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion RIS-CV Simulator", command=RIS_Semi_Inf_Zyl_Ext)
RIS_CVSimulators.add_command(label="External Cylindrical Finite Diffusion RIS-CV Simulator",        command=RIS_Finit_Zyl_Ext)
RIS_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion RIS-CV Simulator",        command=RIS_Finit_Zyl_Int)
RIS_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion RIS-CV Simulator",   command=RIS_Semi_Inf_Sphere_Ext)
RIS_CVSimulators.add_command(label="Internal Spherical Finite Diffusion RIS-CV Simulator",          command=RIS_Finit_Sphere_Int)
RIS_CVSimulators.add_command(label="--------------------------------------------------------------------")
RIS_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion RIS-CV Simulator",               command=RIS_Statistical_Finit_Planar)
RIS_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion RIS-CV Simulator", command=RIS_Statistical_Finit_Zyl_Ext)
RIS_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion RIS-CV Simulator", command=RIS_Statistical_Finit_Zyl_Int)
RIS_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion RIS-CV Simulator",   command=RIS_Statistical_Finit_Sphere_Int)
RIS_CVSimulators.add_command(label="--------------------------------------------------------------------")
RIS_CVSimulators.add_command(label = 'Superimpose different RIS-CVs', command=RISCV_SIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Random-input Signal Voltammetry', menu = RIS_CVSimulators)



#=================================================================================
#=================================================================================
#=================================================================================
#Voltamperometric Evaluation-techniques
#=================================================================================
#=================================================================================
#=================================================================================
EvalVoltamperometric = Menu(menu)
menu.add_cascade(label="Voltamperometric Evaluation", menu=EvalVoltamperometric)
#----------------------------
#Cyclic Voltammetry fitters
#----------------------------
CVFitters   = Menu(menu)
CVFitters.add_command(label="Planar Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Planar_FITTER)
CVFitters.add_command(label="Planar Finite Reflective Diffusion CV Fitter", command=Finit_Planar_FITTER)
CVFitters.add_command(label="Planar Finite Transmissive Diffusion CV Fitter",         command=Finit_Planar_Trans_FITTER)
CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Zyl_Ext_FITTER)
CVFitters.add_command(label="External Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Ext_FITTER)
CVFitters.add_command(label="Internal Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Int_FITTER)
CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Sphere_Ext_FITTER)
CVFitters.add_command(label="Internal Spherical Finite Diffusion CV Fitter", command=Finit_Sphere_Int_FITTER)
CVFitters.add_command(label="--------------------------------------------------------------------")
CVFitters.add_command(label="Statistically Weighted Planar Finite Reflective Diffusion CV Fitter", command=Statistical_Finit_Planar_FITTER)
CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Ext_FITTER)
CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Int_FITTER)
CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Fitter", command=Statistical_Finit_Sphere_Int_FITTER)
CVFitters.add_command(label="--------------------------------------------------------------------")
CVFitters.add_command(label = 'Superimpose different CVs for fitting', command=CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Cyclic Voltammetry Fitters', menu = CVFitters)
#----------------------------
#ACCV-Fourier-Transformer
#----------------------------
EvalVoltamperometric.add_command(label = 'FFT-ACCV-Analysis', command = FFTACCV_Experimental)
#----------------------------
#n-Step Chronoamperometry fitters
#----------------------------
CAFitters   = Menu(menu)
CAFitters.add_command(label="Planar Semi-Infinite Diffusion n-Step CA Fitter",               command=NST_Semi_Inf_Planar_FITTER)
CAFitters.add_command(label="Planar Finite Diffusion n-Step CA Fitter",                      command=NST_Finit_Planar_FITTER)
CAFitters.add_command(label="Planar Finite Transmissive Diffusion n-Step CA Fitter",         command=NST_Finit_Planar_Trans_FITTER)
CAFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion n-Step CA Fitter", command=NST_Semi_Inf_Zyl_Ext_FITTER)
CAFitters.add_command(label="External Cylindrical Finite Diffusion n-Step CA Fitter",        command=NST_Finit_Zyl_Ext_FITTER)
CAFitters.add_command(label="Internal Cylindrical Finite Diffusion n-Step CA Fitter",        command=NST_Finit_Zyl_Int_FITTER)
CAFitters.add_command(label="External Spherical Semi-Infinite Diffusion n-Step CA Fitter",   command=NST_Semi_Inf_Sphere_Ext_FITTER)
CAFitters.add_command(label="Internal Spherical Finite Diffusion n-Step CA Fitter",          command=NST_Finit_Sphere_Int_FITTER)
CAFitters.add_command(label="--------------------------------------------------------------------")
CAFitters.add_command(label="Statistically Weighted Planar Finite Diffusion n-Step CA Fitter",               command=NST_Statistical_Finit_Planar_FITTER)
CAFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion n-Step CA Fitter", command=NST_Statistical_Finit_Zyl_Ext_FITTER)
CAFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion n-Step CA Fitter", command=NST_Statistical_Finit_Zyl_Int_FITTER)
CAFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion n-Step CA Fitter",   command=NST_Statistical_Finit_Sphere_Int_FITTER)
CAFitters.add_command(label="--------------------------------------------------------------------")
CAFitters.add_command(label = 'Superimpose different n-Step CAs for fitting', command=NST_CAFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'n-Step Chronoamperometry Fitters', menu = CAFitters)
#----------------------------
#Large Sine Amplitude Cyclic Voltammetry simulators
#----------------------------
LSA_CVFitters   = Menu(menu)
LSA_CVFitters.add_command(label="Planar Semi-Infinite Diffusion LSA-CV Fitter",               command=LSA_Semi_Inf_Planar_FITTER)
LSA_CVFitters.add_command(label="Planar Finite Reflective Diffusion LSA-CV Fitter",           command=LSA_Finit_Planar_FITTER)
LSA_CVFitters.add_command(label="Planar Finite Transmissive Diffusion LSA-CV Fitter",         command=LSA_Finit_Planar_Trans_FITTER)
LSA_CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion LSA-CV Fitter", command=LSA_Semi_Inf_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="External Cylindrical Finite Diffusion LSA-CV Fitter",        command=LSA_Finit_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="Internal Cylindrical Finite Diffusion LSA-CV Fitter",        command=LSA_Finit_Zyl_Int_FITTER)
LSA_CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion LSA-CV Fitter",   command=LSA_Semi_Inf_Sphere_Ext_FITTER)
LSA_CVFitters.add_command(label="Internal Spherical Finite Diffusion LSA-CV Fitter",          command=LSA_Finit_Sphere_Int_FITTER)
LSA_CVFitters.add_command(label="--------------------------------------------------------------------")
LSA_CVFitters.add_command(label="Statistically Weighted Planar Finite Diffusion LSA-CV Fitter",               command=LSA_Statistical_Finit_Planar_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion LSA-CV Fitter", command=LSA_Statistical_Finit_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion LSA-CV Fitter", command=LSA_Statistical_Finit_Zyl_Int_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion LSA-CV Fitter",   command=LSA_Statistical_Finit_Sphere_Int_FITTER)
LSA_CVFitters.add_command(label="--------------------------------------------------------------------")
LSA_CVFitters.add_command(label = 'Superimpose different LSA-CVs for fitting', command=LSA_CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Large Sine Amplitude Voltammetry Fitters', menu = LSA_CVFitters)
4
#----------------------------
#Staircase Cyclic Voltammetry simulators
#----------------------------
STC_CVFitters   = Menu(menu)
STC_CVFitters.add_command(label="Planar Semi-Infinite Diffusion STC-CV Fitter",               command=STC_Semi_Inf_Planar_FITTER)
STC_CVFitters.add_command(label="Planar Finite Reflective Diffusion STC-CV Fitter",           command=STC_Finit_Planar_FITTER)
STC_CVFitters.add_command(label="Planar Finite Transmissive Diffusion STC-CV Fitter",         command=STC_Finit_Planar_Trans_FITTER)
STC_CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion STC-CV Fitter", command=STC_Semi_Inf_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="External Cylindrical Finite Diffusion STC-CV Fitter",        command=STC_Finit_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="Internal Cylindrical Finite Diffusion STC-CV Fitter",        command=STC_Finit_Zyl_Int_FITTER)
STC_CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion STC-CV Fitter",   command=STC_Semi_Inf_Sphere_Ext_FITTER)
STC_CVFitters.add_command(label="Internal Spherical Finite Diffusion STC-CV Fitter",          command=STC_Finit_Sphere_Int_FITTER)
STC_CVFitters.add_command(label="--------------------------------------------------------------------")
STC_CVFitters.add_command(label="Statistically Weighted Planar Finite Diffusion STC-CV Fitter",               command=STC_Statistical_Finit_Planar_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion STC-CV Fitter", command=STC_Statistical_Finit_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion STC-CV Fitter", command=STC_Statistical_Finit_Zyl_Int_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion STC-CV Fitter",   command=STC_Statistical_Finit_Sphere_Int_FITTER)
STC_CVFitters.add_command(label="--------------------------------------------------------------------")
STC_CVFitters.add_command(label = 'Superimpose different STC-CVs for fitting', command=STC_CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Staircase Cyclic Voltammetry Fitters', menu = STC_CVFitters)





#=================================================================================
#=================================================================================
#=================================================================================
#Electrrochemical Impedance Spectroscopy
#=================================================================================
#=================================================================================
#=================================================================================

Impedance_Techniques = Menu(menu)
menu.add_cascade(label="Impedance Techniques", menu=Impedance_Techniques)

PEIS_simulators   = Menu(menu)
PEIS_simulators.add_command(label="Planar Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Plan_Imp_Sim)
PEIS_simulators.add_command(label="Planar Finite Transmissive Diffusion Impedance Simulator", command=Fin_Plan_Trans_Imp_Sim)
PEIS_simulators.add_command(label="Planar Finite Reflective Diffusion Impedance Simulator", command=Fin_Plan_Ref_Imp_Sim)
PEIS_simulators.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Zyl_Imp_Sim)
PEIS_simulators.add_command(label="Spherical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Sphe_Imp_Sim)
PEIS_simulators.add_command(label="Spherical Internal Finite Diffusion Impedance Simulator", command=Finit_Int_Sphe_Imp_Sim)
Impedance_Techniques.add_cascade(label = 'Potentiostatic EIS Simulators', menu = PEIS_simulators)

PEIS_fitters   = Menu(menu)
PEIS_fitters.add_command(label="Planar Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Plan_Imp_Fit)
PEIS_fitters.add_command(label="Planar Finite Transmissive Diffusion Impedance Fitter", command=Fin_Plan_Trans_Imp_Fit)
PEIS_fitters.add_command(label="Planar Finite Reflective Diffusion Impedance Fitter", command=Fin_Plan_Ref_Imp_Fit)
PEIS_fitters.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Zyl_Imp_Fit)
PEIS_fitters.add_command(label="Spherical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Sphe_Imp_Fit)
PEIS_fitters.add_command(label="Spherical Internal Finite Diffusion Impedance Fitter", command=Finit_Int_Sphe_Imp_Fit)
Impedance_Techniques.add_cascade(label = 'Potentiostatic EIS Fitters', menu = PEIS_fitters)

DRT_Transformation  = Menu(menu)
DRT_Transformation.add_command(label="DRT-Tools-NNLS-DRT", command=DRT_Tools_NNLS_DRT)
DRT_Transformation.add_command(label="JSTT-NNLS-DRT", command=JSTTNNLS_DRT)
Impedance_Techniques.add_cascade(label = 'DRT-Analysis', menu = DRT_Transformation)



#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================

mainloop()

