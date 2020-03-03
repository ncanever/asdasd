
# coding: utf-8

# In[ ]:

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
global File_Was_Loaded


F = 96485.0
R = 8.314
File_Was_Loaded = 0

def cot(phi):
    return 1.0/tan(phi)

def csc(phi):
    return 1.0/sin(phi)

def coth(x):
    return 1/tanh(x)


def No_Loaded_File_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("No loaded file found")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! You loaded no file that\nshould be analyzed. Go to\nopen file first and select\na file that should be analyzed.")


# In[ ]:

def Open_ImpDRT_File():
    root = Toplevel()
    root.title("Your Impedance Data")
    global File_Was_Loaded
    File_Was_Loaded = 1
    global ExpData   #needed later to differentiate experimental and simulated data of ACCV
    ExpData      = 1
    Delim        = Delimiter
    global data
    
    
    if Delim == 1:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t')

        
    if Delim == 0:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='')
    
    
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    
    f  = Figure(figsize=(9, 4), dpi=100)        
    b1 = f.add_subplot(121)
    b2 = f.add_subplot(122)   
        
    
    #============================================================================
    #define entry window
    #============================================================================
    #============================================================================
    
    global Potenzial
    global Strom
    global Frequency_Array
    global Time
    
   
    #==================================================================================
    Frequency         =  data[FromRowxx:ToRowxx:Readeveryxx,0:1:1] 
    Potenzial         =  data[FromRowxx:ToRowxx:Readeveryxx,1:2:1] * UmrPot            
    Strom             =  data[FromRowxx:ToRowxx:Readeveryxx,2:3:1] * UmRStrom
        
    #==================================================
    #turn frequencies alsways from high to low
        
    if Frequency[1] > Frequency[0]:            
        Frequency     = Frequency[::-1]
        Potenzial     = Potenzial[::-1]             
        Strom         = Strom[::-1]
   
    #==================================================================================                 
    
    if  Potenzial[0,0] > Potenzial[1,0]:
        Potenzial  = Potenzial[::-1]
        Strom      = Strom[::-1]
            
            
    global Weite
    Weite = Potenzial.shape[0]
    global NumMeas
    NumMeas = Potenzial.shape[1]
    
    for i in range(Weite):
        Stromarrays             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))                          
        Potenzialarrays         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:])) 
        LenPotArrays            = len(Potenzialarrays)
        Frequency_Array         = Frequency  
        
       
    #=========================================================
    #final plotting
    #=========================================================
    
    b1.set_title("Nyquist-Plot")
    b1.plot(Potenzial, -Strom, linestyle='-',marker='',color='k')
    b1.set_xlabel('Z_real / Ohm', fontsize=12)
    b1.set_ylabel('-Z_imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            

    b2.set_title("Bode-Plot")    
    b2.plot(np.log10(Frequency), -Strom, linestyle='-',marker='',color='b', label = '-Z.imag')
    b2.plot(np.log10(Frequency), Potenzial, linestyle='-',marker='',color='r',label = 'Z.real')
    b2.set_xlabel('log10(freq /s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm  ', fontsize=12)
    b2.plot()
    b2.legend()

    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
    #=========================================================  

         
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)





# In[ ]:

def Get_ImpDRT_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Data")                         
    Fenster.geometry("400x270")                                            

    FromRowxxx_label      = Label(Fenster,text="Read from row")
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)                                               
    FromRowxxx_Eingabe.grid(row=0, column=1)

    ToRowxxx_label       = Label(Fenster,text="Read to row")
    ToRowxxx_label.grid(row=1, column=0)
    ToRowxxx_Eingabe     = Entry(Fenster)
    ToRowxxx_Eingabe.grid(row=1, column=1)
    
    Readeveryxxx_label   = Label(Fenster,text="Read every")
    Readeveryxxx_label.grid(row=2, column=0)
    Readeveryxxx_Eingabe = Entry(Fenster)
    Readeveryxxx_Eingabe.grid(row=2, column=1)
    
    
    UmrPot_Label = Label(Fenster,text="Z_real Factor to be Ohm")
    UmrPot_Label.grid(row=3, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.grid(row=3, column=1)
        
    UmRStrom_Label = Label(Fenster,text="Z_imag Factor to be Ohm")
    UmRStrom_Label.grid(row=4, column=0)
    UmRStrom_Eingabe = Entry(Fenster)
    UmRStrom_Eingabe.grid(row=4, column=1)
    
    HasToBe_Label  = Label(Fenster,text="Data Order")
    HasToBe_Label.grid(row=6, column=0)
    HasToBe_Label1  = Label(Fenster,text="Freq___Z.real___Z.Imag")
    HasToBe_Label1.grid(row=6, column=1)
        
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)



    def Accept_ImpDRT():
              
        global FromRowxx
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        global ToRowxx
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        global Readeveryxx
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        global UmRStrom#
        UmRStrom      = (float(UmRStrom_Eingabe.get()))
        global UmrPot#
        UmrPot        = (float(UmrPot_Eingabe.get()))
        
        global DesiReactxx
        global FirstComesxx

        global Delimiter
        Delimiter         = var5.get()
    
    
    def Next_ImpDRT():
        Open_ImpDRT_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=Accept_ImpDRT)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=Next_ImpDRT)
    Next.grid(row=9, column=0)  

# In[ ]:

#============================================================================================================   
#                 -------> JSTT_NNLS_DRT <-----------
#============================================================================================================   
def JSTTNNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        JSTTNNLS_DRT_Entry_Window()
        
    
      
def JSTTNNLS_DRT_Entry_Window():         
    Fenster = Toplevel()  
    Fenster.geometry("350x250")
    Fenster.title("JSTT-NNLS-DRT")
    
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*")
    Resol_Label.grid(row=0, column=0)
    Resol_Eingabe = Entry(Fenster)                                               
    Resol_Eingabe.grid(row=0, column=1)
        
    reg_par_Label = Label(Fenster, text="Tikhonov regul. parameter*")
    reg_par_Label.grid(row=1, column=0)
    reg_par_Eingabe = Entry(Fenster)                                               
    reg_par_Eingabe.grid(row=1, column=1)
    
    Gauss_decay_Label = Label(Fenster,text="Gaussian-decay-factor*")
    Gauss_decay_Label.grid(row=2, column=0)
    Gauss_decay_Eingabe = Entry(Fenster)
    Gauss_decay_Eingabe.grid(row=2, column=1)
    
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Use imaginary part", variable=var1).grid(row=3, column=0, sticky=W)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Use real part", variable=var2).grid(row=4, column=0, sticky=W)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Save_Data", variable=var3).grid(row=5, column=0, sticky=W)
    
    def AcceptParams():
        global Resol
        global reg_par
        global Gauss_decay
        global Imag_User
        global Real_User
        global As_txt_saver
    

        
        Resol             = (float(Resol_Eingabe.get()))
        reg_par           = (float(reg_par_Eingabe.get()))
        Gauss_decay       = (float(Gauss_decay_Eingabe.get()))
       
        Imag_User         = var1.get()
        Real_User         = var2.get()
        As_txt_saver      = var3.get()
        
        
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=6,column=0)
    button=Button(Fenster,text="Next",command=JSTTNNLS_DRT_Transformation).grid(row=7,column=0)
    
    

    
#========================================================================================================
#========================================================================================================
def JSTTNNLS_DRT_Transformation():
    
    root = Toplevel()
    root.title("DRT-Function")
    
    
    global Z_imag_Meas
    global Z_real_Meas
    global frequencies
    Z_imag_Meas     = np.squeeze(Strom)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Potenzial)    #deshalb hier die Umbenennung
    frequencies     = Frequency_Array
    
    frequencies     = Frequency_Array  #data[:,0]
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    frequencies     = frequencies[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    
    
    
    #==============================================
    #Interpolate higher resolution
    #==============================================
    
    Interpolation_Imag = InterpolatedUnivariateSpline(np.log10(frequencies[::-1]),Z.imag[::-1])
    Interpolation_Real = InterpolatedUnivariateSpline(np.log10(frequencies[::-1]),Z.real[::-1])
    Resolutor          = Resol     #calculate n-times resolution
    Freq_Resol         = np.logspace(np.log10(frequencies[-1]),np.log10(frequencies[0]),Resolutor*len(frequencies))
    Z_Resol_Imag       = Interpolation_Imag(np.log10(Freq_Resol))
    Z_Resol_Real       = Interpolation_Real(np.log10(Freq_Resol))
    frequencies        = Freq_Resol[::-1]
    Z                  = Z_Resol_Real[::-1] + 1j*Z_Resol_Imag[::-1]
    #==============================================
    
    def make_tau(f_max,f_min,M):
        tau = np.zeros(M)
        tau[0] = 1/(2*np.pi*f_max) #first value in array tau is tau_min
        tau[M-1] = 1/(2*np.pi*f_min) #last value is tau_max
        
        for k in range(1,M-1):
            tau[k] =10**(np.log10(tau[0])+(k)/float(M-1)*np.log10(tau[-1]/tau[0]))
        return tau
    
    
    
    global tau
    tau =make_tau(10**8,10**(-8),len(frequencies))

    
    #==========================================================================================

    def build_matrix(freq,tau):
        '''
        set up the data matrix X needed to compute vector w
        which fullfills Xw=y with y being the experimental Z-data
        '''    
        X = np.zeros((len(freq),len(tau)),dtype=np.complex)
        for m in range(len(freq)):
            for n in range(len(tau)):
                X[m][n] = 1/(1+2j*np.pi*freq[m]*tau[n])
        return X

    X = build_matrix(frequencies,tau)
    
    def compute_DRT(X,data,damp):
        
        A = np.vstack((X,damp*np.ones(X.shape)))
        c = np.squeeze(np.concatenate([data,np.zeros(data.shape)]))

        b = nnls(A,c)[0]
        return b

    
    if Imag_User == 1:
        Im_DRT = compute_DRT(-X.imag,np.abs(Z.imag),reg_par)
    
    if Real_User == 1:
        Re_DRT = compute_DRT(X.real,np.abs((Z.real-Z.real[0])),reg_par) 
    
    
    #==========================================================================================
    #Gaussige Verschmierung
    #==========================================================================================
    Knechtarray_Re = np.zeros((len(tau),len(tau)))
    Knechtarray_Im = np.zeros((len(tau),len(tau)))
    
    for i in range(len(tau)):
        for j in range(len(tau)):
            epsilon = Gauss_decay
            if Real_User == 1:
                Knechtarray_Re[i,j] = Re_DRT[i]*np.exp(-epsilon*(np.abs(np.log10(tau[i])-np.log10(tau[j]))**2)) 
            if Imag_User == 1:
                Knechtarray_Im[i,j] = Im_DRT[i]*np.exp(-epsilon*(np.abs(np.log10(tau[i])-np.log10(tau[j]))**2)) 
    
    global broadened_Re
    global broadened_Im
    broadened_Re = np.zeros(len(tau))
    broadened_Im = np.zeros(len(tau))
    
    for i in range(len(tau)):
        if Real_User == 1:
            broadened_Re[i] = np.sum(Knechtarray_Re[::,i])
        if Imag_User == 1:
            broadened_Im[i] = np.sum(Knechtarray_Im[::,i])

    
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    if Imag_User == 1:
        b.plot(np.log10(tau[::]),broadened_Im,color='b', linestyle ='-',label = 'Im_DRT')
    if Real_User == 1:
        b.plot(np.log10(tau[::]),broadened_Re,color='r', linestyle ='-', label = 'Re_DRT')
    b.legend()    
    b.set_xlabel('log10(tau / s)', fontsize=12)
    b.set_ylabel('gamma(tau)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if As_txt_saver == 1:
        JSTTNNLS_DRT_as_Txt_Saver()
    

def JSTTNNLS_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))

    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        fi.write("Resol")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(reg_par))
        fi.write("\n")
        fi.write("Gauss_Decay")
        fi.write("\t")
        fi.write(str(Gauss_decay))
        fi.write("\n")
        
        fi.write("\n")
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\n")
        
        for i in range(len(tau)):
            fi.write(str(np.asscalar(np.log10(tau[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(broadened_Im[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(broadened_Re[i])))
                fi.write("\t")
            fi.write("\n")
    
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_Meas [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_Meas [Ohm]")
        fi.write("\t")
        fi.write("\n")
            
        for i in range(len(Z_real_Meas)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\t")
            fi.write("\n")
    
    root.destroy()


# In[ ]:

#============================================================================================================   
#                 -------> DRT_TOOLS_NNLS_DRT <-----------
#============================================================================================================ 
def DRT_Tools_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        DRT_Tools_NNLS_DRT_Entry_Window()

#============================================================================================================   
#    -------> define all DRT Tools functions in python 2.7 language except quadprog <-----------
#       as there is not quadprog in python this transfom uses nnls algorthm in the end
#============================================================================================================ 

        
def quad_format_combined(A_re,A_im,b_re,b_im,M_re,M_im,damp): 
    H = 2*(0.5*(np.dot(A_re.T,A_re)+np.dot(A_im.T,A_im))+damp*M_re)
    c = -2*0.5*(np.dot(b_im.T,A_im)+np.dot(b_re.T,A_re))
    return H , c

def quad_format(A,b,M,damp):
    H = 2*(np.dot(A.T,A)+damp*M)
    c = -2*np.dot(b.T,A)
    return H , c    

##############################################################################
#map array to gamma
##############################################################################

def map_array_to_gamma(freq_map, freq_coll, x, epsilon, rbf_type):
    if rbf_type == "gaussian":
        rbf = lambda y,y0: np.exp(-(epsilon*(y-y0))**2)
    elif rbf_type == "C0_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))
    elif rbf_type == "C2_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))*(1+abs(epsilon*(y-y0)))
    elif rbf_type == "C4_matern":
        rbf = lambda y,y0: 1/3.0*np.exp(-abs(epsilon*(y-y0)))*(3+3*abs(epsilon*(y-y0))+abs(epsilon*(y-y0))**2)
    elif rbf_type == "C6_matern":
        rbf = lambda y,y0: 1/15.0*np.exp(-abs(epsilon*(y-y0)))*(15+15*abs(epsilon*(y-y0))+6*abs(epsilon*(y-y0))**2+abs(epsilon*(y-y0))**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda y,y0: 1.0/(1+(epsilon*(y-y0))**2)
    else:
        print('ERROR - Unexpected RBF input at map_array_to_gamma')
    y0  = -np.log(freq_coll)
    out_gamma = np.zeros(len(freq_map))   
    
    
    for iter_freq_map in range(len(freq_map)):
        freq_map_loc = freq_map[iter_freq_map]
        y = -np.log(freq_map_loc)
        rbf_temp = rbf(y,y0)
        out_gamma[iter_freq_map] = np.dot(x.T,rbf_temp)
    
    return out_gamma

##############################################################################
#inner products
##############################################################################
def inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)
    if rbf_type == "gaussian":
        out_IP = epsilon**3*(3-6*a**2+a**4)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon**3*(1+abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon**3/6.0*(3+3*abs(a)-6*abs(a)**2+abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon**3/30.0*(45+45*abs(a)-15*abs(a)**3-5*abs(a)**4+abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon**3/140.0*(2835+2835*abs(a)+630*abs(a)**2-315*abs(a)**3-210*abs(a)**4-42*abs(a)**5+abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 48*(16+5*a**2*(-8 + a**2))*np.pi*epsilon**3/((4 + a**2)**5)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf_2')

def inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)

    if rbf_type == "gaussian":
        out_IP = -epsilon*(-1+a**2)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon*(1-abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon/6.0*(3+3*abs(a)-abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon/30.0*(105+105*abs(a)+30*abs(a)**2-5*abs(a)**3-5*abs(a)**4-abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon/140.0*(10395 +10395*abs(a)+3780*abs(a)**2+315*abs(a)**3-210*abs(a)**4-84*abs(a)**5-14*abs(a)**6-abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 4*epsilon*(4-3*a**2)*np.pi/((4+a**2)**3)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf.')
    
#############################################################################
#g_i
#############################################################################
def g_i(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)
    else:
        print('ERROR - Unexpected RBF input at g_i')
    
    integrand_g_i = lambda x: 1.0/(1+alpha**2 *np.exp(2*x))*rbf(x)
    
    out_val  = quad(integrand_g_i, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]#,'RelTol',1E-9,'AbsTol',1e-9);       
    return out_val 


#############################################################################
#g_ii
#############################################################################

def g_ii(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)

    else:
        print('ERROR - Unexpected RBF input at g_ii')
    
    integrand_g_ii = lambda x: alpha/(1.0/np.exp(x)+alpha**2 *np.exp(x))*rbf(x)

    out_val  = quad(integrand_g_ii, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]       
    return out_val 
    
#############################################################################
#compute_L_re
#############################################################################
def compute_L_re(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq+1))
    
    for p in range(N_freq-1):
        delta_loc          = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p+1]  = -1/delta_loc
        out_L_temp[p,p+2]  =  1/delta_loc
        
    out_L[::,1::] = out_L_temp   
    return out_L

#############################################################################
#compute_L_im
############################################################################
def compute_L_im(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq))
    
    for p in range(N_freq-1):
        delta_loc        = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p]  = -1/delta_loc
        out_L_temp[p,p+1]=  1/delta_loc
        
    out_L[::,2::] = out_L_temp   
    return out_L

#############################################################################
#assemble A
############################################################################

def assemble_A_im(freq, epsilon, rbf_type,L=0):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_im_temp = np.zeros(len(freq))
    out_A_im      = np.zeros((len(freq), len(freq)+2))
    if std_freq/mean_freq < 1:  
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_ii(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type)

        out_A_im_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_im_temp[iter_freq_n, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type) 
    out_A_im[:, 2::] = out_A_im_temp
    if L==1:
        out_A_im[:,0] = -2*np.pi*(freq[:])         
    return out_A_im


def assemble_A_re(freq, epsilon, rbf_type):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_re_temp = np.zeros(len(freq))
    out_A_re      = np.zeros((len(freq), len(freq)+2))
    
    if std_freq/mean_freq < 1:  #(error in frequency difference <1% make sure that the terms are evenly distributed)
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_i(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)

        out_A_re_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_re_temp[iter_freq_n, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)
    
    out_A_re[:, 2::] = out_A_re_temp
    out_A_re[:,1] = 1
    return out_A_re
    
##############################################################################
#assemble M
##############################################################################
def assemble_M_re(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_re_temp = np.zeros(len(freq))
    out_M_re      = np.zeros((len(freq)+2, len(freq)+2))
    
    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_re[2::, 2::] = out_M_re_temp
    
    return out_M_re


def assemble_M_im(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_im_temp = np.zeros(len(freq))
    out_M_im      = np.zeros((len(freq)+2, len(freq)+2))

    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
            
            out_M_im_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_im_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_im[2::, 2::] = out_M_im_temp
    
    return out_M_im

##############################################################################
#----------------------------------------------------------------------------#
##############################################################################
    
def compute_A_im(freq,L):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_im = np.zeros((N_freq,N_freq+2))
    if L == 1:
        out_A_im[::,0] = -2*np.pi*freq[::]
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    
    return out_A_im

def compute_A_re(freq):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_re = np.zeros((N_freq,N_freq+2))
    out_A_re[::,1] = 1
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    return out_A_re

#=====================================================================================================================
#translation of all DRT functions from matlab finished. The main-window will be defined as separate function
#=====================================================================================================================
        
    
def DRT_Tools_NNLS_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("350x450")
    Fenster.title("DRT_Tools_NNLS_DRT")
    
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*")
    Resol_Label.grid(row=0, column=0)
    Resol_Eingabe = Entry(Fenster)                                               
    Resol_Eingabe.grid(row=0, column=1)
        
    reg_par_Label = Label(Fenster, text="Tikhonov regul. parameter*")
    reg_par_Label.grid(row=1, column=0)
    reg_par_Eingabe = Entry(Fenster)                                               
    reg_par_Eingabe.grid(row=1, column=1)
    
    epsilon_Label = Label(Fenster,text="epsilon*")
    epsilon_Label.grid(row=2, column=0)
    epsilon_Eingabe = Entry(Fenster)
    epsilon_Eingabe.grid(row=2, column=1)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Use imaginary part", variable=var1).grid(row=3, column=0, sticky=W)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Use real part", variable=var2).grid(row=4, column=0, sticky=W)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Combined re_im", variable=var3).grid(row=5, column=0, sticky=W)
    
    Sep_Label = Label(Fenster,text="Choose ONE type of radial basis function*")
    Sep_Label.grid(row=6, column=0)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="gaussian", variable=var4).grid(row=7, column=0, sticky=W)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="C2_Matern", variable=var5).grid(row=8, column=0, sticky=W)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="C4_Matern", variable=var6).grid(row=9, column=0, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="C6_Matern", variable=var7).grid(row=10, column=0, sticky=W)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Inv._quadratic", variable=var8).grid(row=11, column=0, sticky=W)
    
    Sep2_Label = Label(Fenster,text="Choose ONE type of derivative order*")
    Sep2_Label.grid(row=12, column=0)
    
    var9 = IntVar()
    Checkbutton(Fenster, text="1st_order", variable=var9).grid(row=13, column=0, sticky=W)
    
    var10 = IntVar()
    Checkbutton(Fenster, text="2nd_order", variable=var10).grid(row=14, column=0, sticky=W)
    
    var11 = IntVar()
    Checkbutton(Fenster, text="Save_Data", variable=var11).grid(row=15, column=0, sticky=W)
    
    
    def AcceptParams():
        global Resol
        global Tikh_Par
        global eps
        global Imag_User
        global Real_User
        global Comb_User
        global deru
        global RBF_Type
        global As_txt_saver
    

        Resol             = (float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        eps               = (float(epsilon_Eingabe.get()))
       
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        gaussrbf          = var4.get()
        c2rbf             = var5.get()
        c4rbf             = var6.get()
        c6rbf             = var7.get()
        invqurbf          = var8.get()
        deru1             = var9.get()
        deru2             = var10.get()
        
        #rbf_type
        if gaussrbf == 1:
            RBF_Type = 'gaussian'
        if c2rbf == 1:
            RBF_Type = 'C2_matern'
        if c4rbf == 1:
            RBF_Type = 'C4_matern'
        if c6rbf == 1:
            RBF_Type = 'C6_matern'
        if invqurbf  == 1:
            RBF_Type = 'inverse_quadratic'
          
        #der_used
        if deru1  == 1:
            deru = '1st-order'
        if deru2  == 1:
            deru = '2nd-order'
            
    
        As_txt_saver      = var11.get()
        
        
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=16,column=0)
    button=Button(Fenster,text="Next",command=DRT_Tools_NNLS_DRT_Transformation).grid(row=17,column=0)

#====================================================================================================================    
#here is the GUI translation in python 2.7
#====================================================================================================================    


def DRT_Tools_NNLS_DRT_Transformation():
    ##############################################################################
    
    root = Toplevel()
    root.title("DRT_Tools_DRT-Function")
    
    
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    Z_imag_Meas     = np.squeeze(Strom)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Potenzial)    #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array)
                       
    
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
  
    
    ##############################################################################
    #initialise parameters and bounds
    epsilon       = eps
    rbf_type      = RBF_Type
    der_used      = deru
    damp          = Tikh_Par

    taumax        = ceil(np.max(np.log10(1./freq)))+1
    taumin        = floor(np.min(np.log10(1./freq)))-1
    freq_out      = np.logspace(-taumin,-taumax,Resol*len(freq))

    ##############################################################################
    A_im          = assemble_A_im(freq,epsilon,rbf_type)
    M_im          = assemble_M_im(freq,epsilon,rbf_type,der_used)
    H_im, f_im    = quad_format(A_im,-Z.imag,M_im,damp)

    A_re          = assemble_A_re(freq,epsilon,rbf_type)
    M_re          = assemble_M_re(freq,epsilon,rbf_type,der_used)
    H_re, f_re    = quad_format(A_re,(Z.real),M_re,damp)

    H_co, f_co    = quad_format_combined(A_re, A_im,(Z.real), -Z.imag, M_re, M_im,damp)

    x_imag_2      = nnls(H_im,np.abs(f_im))[0]
    x_real_2      = nnls(H_re,np.abs(f_re))[0]
    x_comb_2      = nnls(H_co,np.abs(f_co))[0]
    
    global gamma_imag_fine
    global gamma_real_fine
    global gamma_comb_fine
    global TAUARRAY
    
    TAUARRAY         = 1./freq_out 

    gamma_imag_fine  = map_array_to_gamma(freq_out,freq,x_imag_2[2:],epsilon,rbf_type)
    gamma_real_fine  = map_array_to_gamma(freq_out,freq,x_real_2[2:],epsilon,rbf_type)
    gamma_comb_fine  = map_array_to_gamma(freq_out,freq,x_comb_2[2:],epsilon,rbf_type)
    
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    if Imag_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_imag_fine,color='b', linestyle ='-',label = 'Im_DRT')
    if Real_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_real_fine,color='r', linestyle ='-', label = 'Re_DRT')
    if Comb_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_comb_fine,color='k', linestyle ='-', label = 'Comb_DRT')
    
    b.legend()    
    b.set_xlabel('log10(tau / s)', fontsize=12)
    b.set_ylabel('gamma(tau)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if As_txt_saver == 1:
        DRT_Tools_NNLS_DRT_as_Txt_Saver()
        
        
def DRT_Tools_NNLS_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))

    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        fi.write("DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("epsilon")
        fi.write("\t")
        fi.write(str(eps))
        fi.write("\n")
        fi.write("RBF-Type:")
        fi.write("\t")
        if RBF_Type == 'gaussian':
            fi.write("gaussian")
        if RBF_Type == 'C2_matern':
            fi.write("C2_matern")
        if RBF_Type == 'C4_matern':
            fi.write("C4_matern")
        if RBF_Type == 'C6_matern':
            fi.write("C6_matern")
        if RBF_Type == 'inverse_quadratic':
            fi.write("inverse_quadratic")
        fi.write("\n")
        fi.write("Derivative-Order:")
        fi.write("\t")
        if deru == '1st-order':
            fi.write("1st-order")
        if deru == '2nd-order':
            fi.write("2nd-order")
    
        
        fi.write("\n")
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")
        
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_fine[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_fine[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_fine[i])))
                fi.write("\t")
            fi.write("\n")
    
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_Meas [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_Meas [Ohm]")
        fi.write("\t")
        fi.write("\n")
            
        for i in range(len(Z_real_Meas)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\t")
            fi.write("\n")
    
    root.destroy()

