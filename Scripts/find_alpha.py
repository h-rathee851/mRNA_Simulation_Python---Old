
# coding: utf-8

# In[15]:
"""
find_alpha.py

"""
import numpy as np
from scipy.optimize import brentq
import math
from helper_methods_cy import *
from simulation_up_cy import *


def rho_alpha(alpha,sim):

    """
    Returns the value of calculated total density for a given alpha after running the simulation.

    Inputs:
           alpha <float>: Inisitiation Rate
           sim <Object Simulation> : An instance of the class Simulation.
           
    Outputs:
            density <float>: Total density of the simulation
    
    """
      
    
    #Set alpha  
    sim.ChangeAlpha(alpha)

    sim.Reach_steady()
    sim.Run_Sim()
    
    return sim.getTot_den()


def func_for_brent(alpha,exp_den,*args):
    """
    Returns the value of calculated total density - the experimental density to apply brents alrogithm.

    Inputs:
           alpha <float>: Inisitiation Rate
           exp_den <float> : Experimentally calculated Density
           *args : Additional arguments for rho_alpha. <Object Simulation>
           
    
    """
    
    return rho_alpha(alpha,*args) - exp_den

def plot_rho(low_bound,upper_bound,sim, rho_exp):

    """
    Plot the calculated density against alpha.

    Inputs:
           low_bound <float>: Lower bound for alpha
           upper_bound <float> : Upper bound for alpha
           sim <Object Simulation> : An instance of the class Simulation.
           rho_exp <float> : Experimentally calculated Densit

    
    """
    # Make empty lists
    alpha_ = []
    totden_ = []
    n = 0
    #Specify number of steps
    steps = 5  
    
    for i in np.linspace(low_bound,upper_bound,steps):

        #Add values to lists
        totden_ += [func_for_brent(i,rho_exp,sim)]
        alpha_ += [i]
        n += 1
        print(str(n*100/steps)+'% done')
        
        
    #Plot values    
    plt.plot(alpha_,totden_)
    plt.show()
    
def calcAlpha(args,brac):

    """
    Calculate alpha for a using brents for a given bracket

    Inputs:
           args : Additional arguments for func_for_brent. exp_den <float> and sim <Object Simulation>
           brac list<float> : Bounds for alpha

    Outputs:
           alpha <float> : Inisitiation rate calculated by brents
           *aux : Other output of brents function.

    
    """ 
    
    alpha, *aux = brentq(func_for_brent,*brac,args = args, full_output= True)
    
    return alpha,aux

def calcAvgAndErr_rho(alpha,sim):

    """
    Calculate mean value of and error in calculated rho_ for an alpha

    Inputs:
           alpha <float> : Inisitiation Rate
           sim <Object Simulation> : Upper bound for alpha

    Outputs:
           mean <float> : Mean value of rho_
           error <float> : Error in mean
           rel_error <float> : Relative error

    
    """ 

    #Create list and add values to it
    rho_ = []
    
    for i in range(5):
        
        rho_.append(rho_alpha(alpha,sim))

    
    rho_ = np.array(rho_)
    mean = np.mean(rho_) #Mean
    error = np.std(rho_)/math.sqrt(5) #Error in mean
    rel_error = error*100/mean #Relative error
    
    return mean,error,rel_error

def calcAvgAndErr_alpha(args,brac):

    """
    Calculate mean value of and error in calculated alpha

    Inputs:
           args : Additional arguments for calcAlpha. exp_den <float> and sim <Object Simulation> 
           brac list<float> : Bounds for alpha

    Outputs:
           mean <float> : Mean value of alpha
           error <float> : Error in mean
           rel_error <float> : Relative error

    
    """ 
    alpha_ = []
    
    for i in range(5):
        
        a = calcAlpha(args,brac)
        alpha_.append(a[0])
        
    
    alpha_ = np.array(alpha_)
    mean = np.mean(alpha_)
    error = np.std(alpha_)/math.sqrt(5)
    rel_error = error*100/mean
    
    return mean,error,rel_error
