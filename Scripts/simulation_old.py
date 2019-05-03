
# coding: utf-8



import numpy as np
import math
import random
from helper_methods import *
import matplotlib.pyplot as plt




# B = 18.03 , alpha = 0.1, start position 2 E = 0.001, rel_error < 1%



class Simulation:
    
    def __init__(self,L, rel_tol = 1, sweeps = 0):
        
        self.L = L
        self.sweeps = sweeps
        self.rel_tol = rel_tol
        
        
    def ExtractRates(self,file_name,alpha,beta):
    
        self.w_ = Extract_hopping_rates(open(file_name, 'r'), alpha,beta)
        self.LEN = len(self.w_) - 1
        self.RUNS = self.sweeps * self.LEN
        
    def ChangeAlpha(self,alpha):
        
        self.w_[0] = alpha
        
        return
    def Reach_steady(self):
        
        self.rib_pos_ = np.array([], dtype=int)
        t = 0
        
        no_rib_enter = 0
        no_rib_exit = 0
        
        runs = 0
        
        
        entry_current_new = 0
        
        while True:
            
            t, enter, exit, move =  self.__Sim_run(t)
            
            runs += 1
            
            if enter:
                no_rib_enter += 1
            elif exit:
                no_rib_exit += 1
              
            if runs == 1:
                
                entry_current_new = no_rib_enter/t
           
            
            if self.RUNS == 0 and ((runs/self.LEN)%500) == 0:
               
                entry_current_old = entry_current_new
                entry_current_new = no_rib_enter/t
                   
            
                if self.__checkSteadyState(entry_current_old,entry_current_new):
                
                    break
                
            elif self.RUNS !=0 and runs >=self.RUNS:
                
                break
            
            
           
        if self.RUNS == 0:
            
            self.RUNS = runs
            
            
    def Run_Sim(self):
        
        self.TIME = 0
        self.current = 0
        self.rho_ = np.zeros(self.LEN)
        
        no_rib_enter = 0
    
        
        for i in range(self.RUNS):
            
            self.TIME,enter, *aux = self.__Sim_run(self.TIME,True)
            
            if enter:
                
                no_rib_enter += 1
    
    
        self.current = no_rib_enter/self.TIME
    
        self.rho_ = self.rho_/self.TIME
        self.rho_ = self.rho_[1:]
        

    def __Sim_run(self, t = 0, steady_reached = False):
        
        enter = False
        exit = False
        move = False
        
        pos_w_ = []
    
        pos_jump_ = Find_pos_jumps(self.rib_pos_,self.L) 
        
    
        for i in pos_jump_:
        
            pos_w_ += [self.w_[i]]
        
        tau = Get_tau(pos_w_)
        
        t += tau
    
        index = Get_mu(pos_w_)
        
        if steady_reached:
            
            self.__calcDensity(tau)
    
        if pos_jump_[index] == 0:
        
            self.rib_pos_ = np.insert(self.rib_pos_,0,2)
            enter = True
        
        elif pos_jump_[index] == self.LEN:
        
            self.rib_pos_ = np.delete(self.rib_pos_, -1)
            exit = True
        
        else:
        
            rib_pos_index = np.where(self.rib_pos_ == pos_jump_[index])
        
            rib_pos_index = rib_pos_index[0][0]
        
            self.rib_pos_[rib_pos_index] += int(1)
        
            move = True
    
 
        return t, enter, exit, move



    def __checkSteadyState(self,entry_current_old,entry_current_new):
        
        if abs((entry_current_old - entry_current_new)*100/entry_current_old) < self.rel_tol:
            
            return True
        
        else: 
            
            return False
        
  
    def __calcDensity(self,tau):
        
        for i in self.rib_pos_:
            
            self.rho_[i-1] += tau
        
        return
    
    def getTot_den(self):
        
        return np.sum(self.rho_)/len(self.rho_)

