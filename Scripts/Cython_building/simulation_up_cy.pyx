# cython: profile=True

"""
simulation_up_cy.pyx

"""



#Imports
import numpy as np
import array
from cpython cimport array
cimport numpy as np
import math
import random
from helper_methods_cy import *
import matplotlib.pyplot as plt
from cpython cimport bool
from libc.math cimport log
cimport cython
# ctypedef np.int_t DTYPE_t
@cython.boundscheck(False) # turn off bounds-checking for entire function


cdef class Simulation:
    
    """
    Class containing the main methods to run the simulation
    """

    # Class variables
    cdef public int L, sweeps, RUNS, LEN
    cdef public double rel_tol, TIME, current, a
    cdef array.array  pos_w_,rib_pos_, pos_jump_, w_
 #   cdef public np.ndarray pos_w_,rib_pos_, pos_jump_
    cdef public np.ndarray rho_
  #  cdef dict __dict__
        
    def __init__(self,int L, double  rel_tol = 0.1, int sweeps = 0):
        
        """
        Class initiator
        
        Inputs:
                L <int>: Size of ribosome
                rel_tol <float>: Maximum Relative current difference to determine steady state
                sweeps <int>: Number of sweeps to the the simulation for. (Both to reach
                                steady state and after)
                                
        """
        self.L = L
        self.sweeps = sweeps
        self.rel_tol = rel_tol
        
        return 
        
    def ExtractRates(self,str file_name,float alpha,float beta):
        
         
        """
        Extract hopping rates from the file.
        
        Inputs:
                file_name <string>: Path/Name of file containing the rates
                alpha <float>: Inisitiation Rate
                beta <float>: Exit rate
                                
        """
     #   assert self.w_.dtype == np.double_t
        self.w_ = array.array('d',(Extract_hopping_rates(open(file_name, 'r'), alpha,beta)))
        self.LEN = len(self.w_) - 1
        self.RUNS = self.sweeps * self.LEN
        
    def ChangeAlpha(self,float alpha):
        """
        Change the inisitiation rate
        
        Inputs:
                alpha <float>: Inisitiation Rate
                                
        """
        
        self.w_[0] = alpha
        return
    
    def Reach_steady(self):
        
        """
        Reach the steady state for the simulation
                                
        """

        # Initialize lists needed for simulation
        self.pos_jump_ = array.array('i',[0])# List of ribisome positions that can possible jump in next step
        self.pos_w_ = array.array('d',[self.w_[0]]) # List of hopping rates for ribosomes that can jump
        self.rib_pos_ = array.array('i',[]) # List of positions of ribosomes
        self.a = self.w_[0]
        # Assert their data types
    #    assert self.pos_jump_.dtype == np.int 
    #    assert self.rib_pos_.dtype == np.int
    #    assert self.pos_w_.dtype == np.double
        
        
        #cdef np.ndarray[np.int_t, ndim=1] pos_jump_cy = np.zeros(1, dtype = int)
     #  cdef np.ndarray[np.float_t, ndim=1] rho_
       # cdef np.ndarray[np.int_t, ndim=1] rib_pos_cy = np.array([], dtype = int)
       # aux = [pos_jump_cy,pos_w_cy,rib_pos_cy]

        aux = [self.pos_jump_,self.pos_w_,self.rib_pos_,self.w_]

        # Initialize variables of simulation
        cdef float t = 0
        cdef int no_rib_enter = 0
        cdef int no_rib_exit = 0
        cdef int runs = 0
        cdef float entry_current_new = -1
        cdef float entry_current_old
        cdef bool enter = False
        cdef bool do_exit = False
        cdef bool move = False
        
        while True:
            
            #Run simulation step
            t, enter, do_exit, move =  self.__Sim_run(t)
            
            runs += 1
            
            #Update number of ribosomes entered and exited
            if enter:
                no_rib_enter += 1
            elif do_exit:
                no_rib_exit += 1
            
         
            
            if self.RUNS == 0 and ((runs/self.LEN)%500) == 0:
               
                #Calculate currents
                entry_current_old = entry_current_new
                entry_current_new = no_rib_enter/t

            #    print("Here"+ str(entry_current_old)+str(entry_current_new))
                # Check if Steady state reached if runs not specified
                if self.__checkSteadyState(entry_current_old,entry_current_new):
                    
                
                    break
            
            # If runs != 0 then stop if number of runs reached
            elif self.RUNS !=0 and runs >=self.RUNS:
                
                break
            
            
        # If runs not specified, set runs of the sumulation as the runs it took to reach steady state   
        if self.RUNS == 0:
            
            self.RUNS = runs

            
    def Run_Sim(self):
        
        """
        Run simulation after reaching steady state
                                
        """
        #Run simulation
        

     #   cdef np.ndarray[np.int_t, ndim=1] pos_jump_cy = self.pos_jump_
    #    cdef np.ndarray[np.float_t, ndim=1] test_ = self.pos_w_
     #   cdef np.ndarray[np.int_t, ndim=1] rib_pos_cy = self.rib_pos_
     #   aux = [pos_jump_cy,pos_w_cy,rib_pos_cy]


        aux = [self.pos_jump_,self.pos_w_,self.rib_pos_,self.w_]

        #Initialize simulation variables
        self.TIME = 0
        self.current = 0
        cdef int LEN = self.LEN
        cdef bool enter = False
        cdef bool do_exit = False
        cdef bool move = False
        #cdef np.ndarray[np.double_t, ndim=1] rho_
        self.rho_ = np.zeros(self.LEN)
        assert self.rho_.dtype == np.double 
        
        cdef int no_rib_enter = 0
        cdef int i
        # Run simulation
        for i in range(self.RUNS):
            
            self.TIME, enter, do_exit, move = self.__Sim_run(self.TIME,steady_reached = True)
            
            if enter:
                
                no_rib_enter += 1
            
        
        # Calculate current
        self.current = no_rib_enter/self.TIME
        # Calculate density profile
        self.rho_ = self.rho_/self.TIME
        self.rho_ = self.rho_[1:]


    cdef __Sim_run(self, double t, bool steady_reached = False):

        """
        Simulation step
        
        Inputs:
               t <double> : current real time in simulation
               aux list<np.ndarray> : List of numpy arrays used for simulation
               steady_reached <bool> : Specify if steady state reached

        Outputs:
                t <double> : Real time after simulation step
                enter <bool> : Specify whether a new ribosome entered in the step
                do_exit <bool> : Specify whether a ribosome exited during the step
                move <bool> : Specify whether a ribosome moved during the step
                aux list<np.ndarray> : List of numpy arrays updated during the step

        """
        
        cdef bool enter = False
        cdef bool do_exit = False
        cdef bool move = False
        
        
        # Get time for this step
        cdef double tau = Get_tau(self.a)

        # Update simulation time
        t += tau
        #Get step index for the step
        cdef Py_ssize_t index = Get_mu(self.pos_w_,self.a)
        cdef Py_ssize_t rib_pos_index

        # If steady state reached, calculate density
        if steady_reached:
            self.__calcDensity(tau,self.rho_)
    
        # Inisitiation step
        if self.pos_jump_[index] == 0:

            #Update values
            self.rib_pos_.insert(0,2)
            self.a -= self.w_[0]

            # Check if the new entered ribosime can move in the next step and update values
            if len(self.rib_pos_) > 1:
                if self.rib_pos_[1] > self.L+2: 
                    self.pos_jump_[index] = 2
                    self.pos_w_[index] = self.w_[2]
                    self.a += self.w_[2]
                else:
                    del self.pos_jump_[index]
                    del self.pos_w_[index]
            else:
                self.pos_jump_[index] = 2
                self.pos_w_[index] = self.w_[2]
                self.a += self.w_[2]
                
            enter = True
            
        # Exit step
        elif self.pos_jump_[index] == self.LEN:
            #Update values
            del self.rib_pos_[-1]
            self.a -= self.w_[-1]

            # Check if second last robisome can move in next step and update values
            if len(self.rib_pos_) >= 1 and self.rib_pos_[-1] == self.LEN-self.L:
                self.pos_jump_[index] = self.rib_pos_[-1]
                self.pos_w_[index] = self.w_[self.rib_pos_[-1]]
                self.a += self.w_[self.rib_pos_[-1]]
            
            else:
                del self.pos_jump_[index]
                del self.pos_w_[index]
                
            do_exit = True

        # Move step    
        else:
            #Update values
            rib_pos_index = self.rib_pos_.index(self.pos_jump_[index])
            self.a -= self.pos_w_[index]
            self.rib_pos_[rib_pos_index] += int(1)
            
            #Check if moved ribosome can move in next step and update values
            if len(self.rib_pos_) > rib_pos_index+1:
                if self.rib_pos_[rib_pos_index+1] - self.rib_pos_[rib_pos_index] > self.L:
                    self.pos_jump_[index] = self.rib_pos_[rib_pos_index]
                    self.pos_w_[index] = self.w_[self.rib_pos_[rib_pos_index]]
                    self.a += self.w_[self.rib_pos_[rib_pos_index]]
                else:
                    del self.pos_jump_[index]
                    del self.pos_w_[index]
            else:
                
                self.pos_jump_[index] = self.rib_pos_[rib_pos_index]
                self.pos_w_[index] = self.w_[self.rib_pos_[rib_pos_index]]
                self.a += self.w_[self.rib_pos_[rib_pos_index]]

            # Check if ribosome behing the moved ribosome can move in next step and update values
            if rib_pos_index != 0 and self.rib_pos_[rib_pos_index] - self.rib_pos_[rib_pos_index-1] == self.L+1 :
                
                self.pos_jump_.insert(index,self.rib_pos_[rib_pos_index-1]) 
                self.pos_w_.insert(index,self.w_[self.rib_pos_[rib_pos_index-1]])
                self.a += self.w_[self.rib_pos_[rib_pos_index-1]]
            
                
            elif rib_pos_index == 0 and self.rib_pos_[rib_pos_index] == self.L+2:
                self.pos_jump_.insert(0,0) 
                self.pos_w_.insert(0,self.w_[0])
                self.a += self.w_[0]
        
            move = True
    
 
        return t, enter, do_exit, move



    cdef bool __checkSteadyState(self,float entry_current_old,float entry_current_new):

        """
         Check if steady state reached
    
         Inputs:
              entry_current_old <float> : Old entry current
              entry_current_new <float> : New entry current
         Output
             <bool>

        """
        
        #If entry_current of ribosomes changed by less than relative tolerance, return true else return false
        if abs((entry_current_old - entry_current_new)*100/entry_current_old) < self.rel_tol:
            
            return True
        
        else: 
            
            return False
        
  
    cdef void __calcDensity(self,double tau ,np.ndarray[np.double_t, ndim=1] rho_ ):

        """
        Update density profile after a simulation step

        """
        cdef Py_ssize_t i
        for i in self.rib_pos_:
            
            rho_[i-1] += tau
        
        return
    
    def getTot_den(self):

        """
        Calculate total density given the density profile

        """
        
        return np.sum(self.rho_)/len(self.rho_)



