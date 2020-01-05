# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 18:52:36 2017

@author: BenjaSmea
"""
import numpy as np

def PressureDrop(pipelengtheqv,viscosity,flowrate,head,pipeid,rr,density):
    velocity = flowrate/(np.pi*(pipeid/2)**2)
    reynum = (pipeid * velocity) / viscosity
    if reynum < 2300:
        fricfac = 64/reynum
    else:
        fricfac = 1
        f = 0.000001
        for x in range(0, 1000):
            #Colebrook equation as a function of friction factor
            fx = 1/np.sqrt(f)+2*np.log10((rr/(3.7*pipeid*1000))+(2.51/(reynum*np.sqrt(f))))
            #First derivative of colebrook equation
            fxprime = -0.5*(f**(-3/2))*(1+((2*2.51)/(np.log10((rr/(3.7*pipeid*1000))+2.51/(reynum*np.sqrt(f)))*reynum)))
            if fx >= -0.0001 and fx <= 0.0001:
                fricfac = f
                break
            #Newton Raphson iterative solution formula
            f = f - (fx/fxprime) 
            if x > 999:
                print("Solution invalid equation did not converge")
    Dpfric = fricfac * (pipelengtheqv/pipeid)*((density*velocity**2)/2)*1*10**-5
    Dphead = density*9.81*head*1*10**-5
    Dptotal = Dpfric+Dphead
    return Dptotal, velocity, reynum

fluidkey = str(input("Choose fluid Name (e,g, HFO20, Diesel25, Water): ")) 

viscosity = ""
density = ""
pdtotal = 0
with open("fluidlist.txt") as fluidfile:
    for line in fluidfile:
        if line.split(',')[0] == fluidkey:
            viscosity = float(line.split(',')[3])
            density = float(line.split(',')[1])
            print('Found fluid: ' + fluidkey)
            
pipeno= int(input("No of pipes?: "))   

for x in range(0, pipeno):
    flowrate= float(input("Choose pipe flow rate in l/hr for pipe "  + str(x+1) + ": "))*float(2.778*10**-7)
    pipekey = str(input("Choose Pipe Size for pipe " + str(x+1) + ": (e.g. DN80, DN100): ")) 
    print('Calculating ' + str(x+1) + ' of ' + str(pipeno) + ' pipes in piping system')
    with open("pipelist.txt") as pipefile:
        for line in pipefile:
            if line.split(',')[0] == pipekey:
                pipeid = float(line.split(',')[1])/1000
                print('Found pipe: ' + pipekey)
                
    pipelen = float(input("Pipe length in m for pipe "  + str(x+1) + "?: ")) 

    #fluidlist = np.atleast_1d(np.genfromtxt("fluidlist.txt", delimiter=",", dtype=None,comments='#'))
    #pipelist = np.atleast_1d(np.genfromtxt("pipelist.txt", delimiter=",", dtype=None,comments='#'))
    #print('\n \nList of pipe sizes in file (1-n)\n\n' + str(pipelist))
    #print('\n \nList of fluid types in file (1-n) \n\n' + str(fluidlist))
      
    pipelengtheqv = pipelen
    head = 0
    rr=0.046
    pd, velocity, reynum  = PressureDrop(pipelengtheqv,viscosity,flowrate,head,pipeid,rr,density)
    pdtotal += pd
    print('Pipe ID: ' + str(pipeid) + 'm')
    print('Velocity: ' + str(velocity) + 'm/s')
    print('Flowrate used for pipe '  + str(x+1) + ' is : ' + str(flowrate) + 'm^3/s')
    print('Pressure drop for pipe '  + str(x+1) + " is :" + str(pd) + ' Bar\n')
    
print('Viscosity used is: ' + str(viscosity) + 'm^2/s')
   
print('Total Pressure drop is :' + str(pdtotal) + ' Bar')

input("Press Enter to quit")