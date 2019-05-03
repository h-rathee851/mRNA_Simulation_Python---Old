
# coding: utf-8

# In[4]:

import numpy as np
import math
import random

def Extract_hopping_rates(file,alpha,beta):
    
    
    """
    Extract hopping rates from a given file and add initiation and exit rate
    
    Inputs:
            file <file>: File to read the hopping rates from
            alpha <float>: Inisization rate
            beta <float>: Exit rate
   
    Output:
            w_ <np.array>: Array containing all hopping rates
    
    """
    
    # Read hopping rates and append to list
    lines = file.readlines()
    file.close()
    
    w_ = []
    for line in lines:
        
        codon,w = line.split(' ',1)
        w = float(w)
        w_.append(w)
        
    w_ = [alpha] + w_ + [beta]
    
    return np.array(w_)

def Get_tau(a):
    
    """
    Generate a random time interval according to the distribution a exp(-a*tau). Where a = sum of available hopping rates
    
    Inputs:
            w_ <np.array> : Array containing available hopping rates
            
    outputs:
            tau <float>: Random time interval
    """
    np.random.seed(2019)
    r = random.random()
    #a = np.sum(w_)
    return (1.0/a) * math.log(1.0/r)

def Get_mu(w_,a):
    
    """
    Generate a random hopping event according to the distribution w_i/a. Where a = sum of available hopping rates
    and w_i are individual hopping rates.
    
    Inputs:
            w_ <np.array> : Array containing available hopping rates
            
    outputs:
            tau <float>: Random index of selected hopping
    
    """
    np.random.seed(2019)
    r = random.random()
    w_ = np.cumsum(w_)
    
    for i in range(0,len(w_)):
        
        if w_[i] >= r*a :
            
            return i
             
        
def Find_pos_jumps(rib_pos_,L):
    
    """
    Find ribosimes that can move in the next step of the simulation
    
    Inputs:
            rib_pos_ <np.array> : Array containing current ribosome position
            L <int>: Ribosome size
            
    outputs:
            JUMP_POS_IND_ <np.array>: Indices of ribosomes that can move in the next step of the simulation
    
    """

    if len(rib_pos_) == 0: #Initiate 
        
        return np.array([0])
    
    jump_pos_ind_ = []
    
    if rib_pos_[0] > L+1: # If first L-sites clear, a new ribosime can join
        
        jump_pos_ind_ += [0]
        
    for i in range(0,len(rib_pos_)-1): # If distence beteen ribosomes > L, it can move
        
        if rib_pos_[i+1] - rib_pos_[i] > L:
            
            jump_pos_ind_ += [rib_pos_[i]]
        
    jump_pos_ind_ += [rib_pos_[-1]] # Last ribosome can always move
    
    return np.array(jump_pos_ind_)


# In[ ]:



