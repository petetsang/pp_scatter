#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.pyplot import figure

figure(figsize=(8, 6), dpi=80)

# pylab has sines and cosines
from pylab import *

#######################################################################################
#############################
#######################################################################################
# Plot data from files
import pylab
list_of_files = [
#     ('data/24gev.dat','23.5 GeV'),
#    ('data/31gev.dat','30.7 GeV'),
#    ('data/45gev.dat','44.7 GeV'),
#    ('data/53gev.dat','52.8 GeV'),
#    ('data/63gev.dat','62.5 GeV'),
#     ('data/7tev.dat','7TeV'),
#    ('data/8tev.dat','8TeV'),
#    ('data/13tev.dat','13TeV'),
#    ('data/pip.dat','pip'),
#     ('data/isr_totalcrosssection.dat','isr_total_crosssection'),
#      ('data/lhc_totalcrosssection.dat','lhc_total_crosssection'),
      ('data/total_cross_section_data.dat','total_cross_section')
]

datalist1 = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files ]

for data, label in datalist1:
     pylab.plot(data[:,0], data[:,1],linestyle=' ',  color='black', marker='o', markersize=4  # marker =r'$\cdot$', markersize=15
     
)
     pylab.errorbar(data[:,0], data[:,1], yerr=data[:,2],linestyle=' ', capsize=2, marker='', c='b',ecolor='black',markersize=10
     )

#################################################################################
######### sigma_elastic_total data
#################################################################################
list_of_files = [
    #     ('data/24gev.dat','23.5 GeV'),
    #    ('data/31gev.dat','30.7 GeV'),
    #    ('data/45gev.dat','44.7 GeV'),
    #    ('data/53gev.dat','52.8 GeV'),
    #    ('data/63gev.dat','62.5 GeV'),
    #     ('data/7tev.dat','7TeV'),
    #    ('data/8tev.dat','8TeV'),
    #    ('data/13tev.dat','13TeV'),
    #    ('data/pip.dat','pip'),
    #     ('data/isr_totalcrosssection.dat','isr_total_crosssection'),
    #      ('data/lhc_totalcrosssection.dat','lhc_total_crosssection'),
    #      ('data/total_cross_section_data.dat','total_cross_section'),
          ('data/total_elastic_cross_section.dat','total_elastic_cross_section')
    ]

datalist = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files ]

for data, label in datalist:
         pylab.plot(data[:,0], data[:,1],linestyle=' ',  color='black', marker='*', markersize=4  # marker =r'$\cdot$', markersize=15
         
    )
         pylab.errorbar(data[:,0], data[:,1], yerr=data[:,2],linestyle=' ', capsize=2, marker=' ', c='b',ecolor='black',markersize=10
         )
     
     
#################choose your energy
E0=20.0
E1=23.5
E2=30.7
E3=44.7
E4=52.8
E5=62.5
E6=7000.
E7=8000.0
E8=13000.0
#Choose an energy to calculate
E=E6

#######################################################################################
# define variables m is in units of GeV**2



# t_min=0.005
# K=.389
# g=7.0
# p=-0.1083*1
# lambda1=0.2010855
# kappa=-0.0040392
# m=0.154574
# mbar=0.41328
# N=9
# s = E**2
# Nem=0.0003
# Nem_p= Nem * 10

t_min=0.005*.7
K=.389
g=7.0
p=-0.02*4.1
lambda1=0.4*1.0
kappa=-0.0099*1
m=0.181*1.8
mbar=0.41*.9
N=9
Nem = 0.0003
Nem_p = Nem*10


s = E**2

#######################################################################################
############## DELTA Q ENERGY DEPENDENCE ON VERTICES
#######################################################################################


#delta_q = (lambda1/m)*(1/(math.log(sqrt(s)/6*m)))
#delta_q= (lambda1/m)*(6*m/sqrt(s))**p
#delta_q = (lambda1/m)*(log(sqrt(s)/(6*m)))   #foissart bound related or something
#delta_q = (lambda1/m)*(log(sqrt(s)/(6*m))**p)*1


#######################################################################################
############## ENERGY DEPENDENCE OF DELTA_Q(S)
#######################################################################################


def delta_q(z):
    return (lambda1/m)*(6*m/  (sqrt(z)) )**(-1*p)

def delta_q_plot(z): 
    return (lambda1/m)*(6*m/z)**(1*p)

#######################################################################################
############## TOTAL CROSS-SECTION
#######################################################################################


def term1(z):
    return ( 1/(1*sqrt(2))  ) * g *  delta_q_plot(z)**2

def term2(z):
    return (g**2*kappa**2/(16*pi)) *    (mbar**2) *delta_q_plot(z)**4

def term3(z):
    return (g**3*m**4*kappa/(4*sqrt(2)*pi))   *(mbar**2/(m**2+(mbar**2/2))**2)*delta_q_plot(z)**4

def sigmatotal(x):
    return( K*N*(term1(x) - 0*term2(x) - 0*term3(x)))

def sigmatotal_1bundle(x):
    return ( K *N*(term1(x)) )

def sigmatotal_1bundle_em(z):
    return (
     K*N * ( (2*Nem_p/t_min) + ((1/(1*sqrt(2)))*g * delta_q_plot(z)**2 )
        )
)



##################################################################################
####### sigma_elastic_total
##################################################################################

EulerGamma= .577216



################## Sept 17th, 2024
############### without E&M
def analytic_elastic_sigma(y):
    return(
        N*K/(4*pi)*(
            0*(Nem**2+Nem_p**2)/t_min
            - 0*(1/sqrt(2))*(Nem_p-Nem)*g*delta_q_plot(y)**2*(EulerGamma+log(t_min/(4*m**2)) )
            + 0*Nem_p*g**2*delta_q_plot(y)**2*kappa            
            + (1/2)*g**2*delta_q_plot(y)**4*m**2
            + (1/sqrt(2))*g**3*delta_q_plot(y)**4*kappa*((m**4*mbar**2)/(m**2+(mbar**2)/2)**2)
            + (1/8)*g**4*delta_q_plot(y)**4*kappa**2*mbar**2
            
            )
        
        )



######### with E&M
def analytic_elastic_sigma_em(y):
    return(
        N*K*(
            (Nem**2+Nem_p**2)/t_min
            - (1/sqrt(2))*(Nem_p-Nem)*g*delta_q_plot(y)**2*(EulerGamma+log(t_min/(4*m**2)))
            + Nem_p*g**2*delta_q_plot(y)**2*kappa            
            - 0*(1/(8*pi))*Nem*g**2*delta_q_plot(y)**4*m**2*(EulerGamma + log(t_min/(8*m**2)))
            + (1/2)*g**2*delta_q_plot(y)**4*m**2
            + (1/sqrt(2))*g**3*delta_q_plot(y)**4*kappa*((m**4*mbar**2)/(m**2+(mbar**2)/2)**2)
            + (1/8)*g**4*delta_q_plot(y)**4*kappa**2*mbar**2
            - 0*(1/(6*pi*sqrt(2)))*g**3*delta_q_plot(y)**6*m**4
            + 0*(1/(8*pi)**2)*g**4*delta_q_plot(y)**8*m**6
            
            )
          )
        






#######################################################################################
########### PLOTS
#######################################################################################

t=np.arange(1400,14000,10)
#### sigmatotalold works better, July 31, 2023
y1= sigmatotal_1bundle(t)
#y1= sigmatotal1(t)
######### old sigmatotal_em
y2= sigmatotal_1bundle_em(t)
#y2= sigmatotal_em1(t)
y3= analytic_elastic_sigma(t)
#y4= analytic_elastic_sigma(t)
y4= analytic_elastic_sigma_em(t)
#y5= analytic_elastic_sigma2(t)


plot1, = plt.semilogx(t, y1, color='g',linewidth=2, linestyle="dashed",label=str(E) + r"$GeV^{\ 2}$")
plot2, = plt.semilogx(t, y2, color='r', label=str(E2))
plot3, = plt.semilogx(t, y3, color='b', label=str(E3))
plot4, = plt.semilogx(t, y4, color='m', )
#plot5, = plt.semilogx(t, y5, color='y', )

plt.xlim(10E0, 10e5)
plt.ylim(0,140)

#plt.title(r'Elastic pp-scattering differential cross-section at ' + str(E) + r' $GeV^{\ 2}$')
plt.xlabel(r'$\sqrt{s}\  [GeV]$',size=17)
plt.ylabel(r'$\sigma_{total} (green), \sigma_{elastic} (red)[mb]$',size=17)

legend1= plt.legend(['$\sigma_{total}$ data','$\sigma_{elastic}$ data'
    , '$\sigma_{total}$'
    , '$\sigma_{total}\ with\  EM$' 
    , '$\sigma_{elastic}$'
    , '$\sigma_{elastic} \ with \ EM$'

   

    # 'LHC region, '  
    # , "\n" r"$\sigma_{total}$" 
    # + "\n" " in green dash, " 
    # + "\n" + r"$\sigma_{elastic}$ in red" 
    # + r"$\sqrt{s}=$"+ str(E) + r" $GeV^{}$"
    
  # , 'p = ' + str(p)[:8] 
  
  # + "\n" "g=" + str(g)[:3] 
  # + "\n" "m=" + str(m)[:6] 
  # + "\n" r"$\bar{m}=$" + str(mbar)[:6] 
  # + "\n" r"$\kappa=$" + str(kappa)[:9] 
  # + "\n" r"$\lambda=$" + str(lambda1)[:6]
  # + "\n" "N= " + str(N_)[:3]
  # + "\n" "N'=" + str(N_p)[:3]
#            plot1,
#                'g=' + str(g)
#                  + '\n' + 'p=' + str(p)
#                  + '\n' + 'k=' + str(k) + 'millibarns'
#                  + '\n' + r'$m_{ext}=$' + str(m_ext)
#                  + '\n' + r'$m_{int}=$' + str(m_int) 
#                  + '\n' + r'$\beta_{ext}=$' + str(beta_ext)
#                  + '\n' + r'$\beta_{int}=$' + str(beta_int)
#                  + '\n' + r' $\bar{k}$=' + str(lambda_kappa_bar),
             #plot2
             #plot3,
             #plot4,
             #plot5,
             #plot6,
             #plot7, 'g=' + str(g) + ', p=' + str(p)  + '\n k=' + str(k) + ' millibarns' + '\n m='  + str(m) + r'$GeV^{\ 2}$' + '\n' + r' $\beta_{ext}=$' + str(beta) + '\n' + r' $\lambda\bar{k}$=' + str(lambda_kappa_bar)
 #             ],fontsize=9,loc='center left', bbox_to_anchor=(1, 0.5))
               ],loc=2,fontsize=9)
#legend2 = plt.legend(['test'])

plt.gca().add_artist(legend1)

plt.grid(True)

plt.show(block=True)