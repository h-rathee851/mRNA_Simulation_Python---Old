# cython: profile=True

"""
simulation_up_cy.py

"""



#Imports
import numpy as np
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
    cdef public np.ndarray pos_w_,rib_pos_, pos_jump_
    cdef public np.ndarray w_, rho_
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
        self.w_ = Extract_hopping_rates(open(file_name, 'r'), alpha,beta)
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
        self.pos_jump_ = np.zeros(1,dtype=int)# List of ribisome positions that can possible jump in next step
        self.pos_w_ = np.array([self.w_[0]]) # List of hopping rates for ribosomes that can jump
        self.rib_pos_ = np.array([],dtype = int) # List of positions of ribosomes
        self.a = self.w_[0]
        # Assert their data types
        assert self.pos_jump_.dtype == np.int 
        assert self.rib_pos_.dtype == np.int
        assert self.pos_w_.dtype == np.double
        
        
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
            t, enter, do_exit, move, aux =  self.__Sim_run(t,aux)
            
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

        # Update class lists
        self.pos_jump_ = aux[0]
        self.pos_w_ = aux[1]
        self.rib_pos_ = aux[2]
            
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
            
            self.TIME, enter, do_exit, move, aux = self.__Sim_run(self.TIME,aux,steady_reached = True)
            
            if enter:
                
                no_rib_enter += 1
            
        
        # Calculate current
        self.current = no_rib_enter/self.TIME
        # Calculate density profile
        self.rho_ = self.rho_/self.TIME
        self.rho_ = self.rho_[1:]

        # Update Class lists
        self.pos_jump_ = aux[0]
        self.pos_w_ = aux[1]
        self.rib_pos_ = aux[2]

    cdef __Sim_run(self, double t, aux, bool steady_reached = False):

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

        cdef np.ndarray[np.int_t, ndim=1] pos_jump_cy = aux[0]
        cdef np.ndarray[np.double_t, ndim=1] pos_w_cy = aux[1]
        cdef np.ndarray[np.int_t, ndim=1] rib_pos_cy = aux[2]
        cdef np.ndarray[np.double_t, ndim=1] w_cy = aux[3]
        
        cdef bool enter = False
        cdef bool do_exit = False
        cdef bool move = False
        
        
        # Get time for this step
        cdef double tau = Get_tau(self.a)

        # Update simulation time
        t += tau
        #Get step index for the step
        cdef Py_ssize_t index = Get_mu(pos_w_cy,self.a)
        cdef Py_ssize_t rib_pos_index

        # If steady state reached, calculate density
        if steady_reached:
            self.__calcDensity(tau,rib_pos_cy,self.rho_)
    
        # Inisitiation step
        if pos_jump_cy[index] == 0:

            #Update values
            rib_pos_cy = np.insert(rib_pos_cy,0,2)
            self.a -= w_cy[0]

            # Check if the new entered ribosime can move in the next step and update values
            if len(rib_pos_cy) > 1:
                if rib_pos_cy[1] > self.L+2: 
                    pos_jump_cy[index] = 2
                    pos_w_cy[index] = w_cy[2]
                    self.a += w_cy[2]
                else:
                    pos_jump_cy= np.delete(pos_jump_cy, index)
                    pos_w_cy= np.delete(pos_w_cy, index)
            else:
                pos_jump_cy[index] = 2
                pos_w_cy[index] = w_cy[2]
                self.a += w_cy[2]
                
            enter = True
            
        # Exit step
        elif pos_jump_cy[index] == self.LEN:
            #Update values
            rib_pos_cy = np.delete(rib_pos_cy, -1)
            self.a -= w_cy[-1]

            # Check if second last robisome can move in next step and update values
            if len(rib_pos_cy) >= 1 and rib_pos_cy[-1] == self.LEN-self.L:
                pos_jump_cy[index] = rib_pos_cy[-1]
                pos_w_cy[index] = w_cy[rib_pos_cy[-1]]
                self.a += w_cy[rib_pos_cy[-1]]
            
            else:
                pos_jump_cy= np.delete(pos_jump_cy, index)
                pos_w_cy= np.delete(pos_w_cy, index)
                
            do_exit = True

        # Move step    
        else:
            #Update values
            rib_pos_index = np.where(rib_pos_cy == pos_jump_cy[index])[0][0]
            self.a -= pos_w_cy[index]
            rib_pos_cy[rib_pos_index] += int(1)
            
            #Check if moved ribosome can move in next step and update values
            if len(rib_pos_cy) > rib_pos_index+1:
                if rib_pos_cy[rib_pos_index+1] - rib_pos_cy[rib_pos_index] > self.L:
                    pos_jump_cy[index] = rib_pos_cy[rib_pos_index]
                    pos_w_cy[index] = w_cy[rib_pos_cy[rib_pos_index]]
                    self.a += w_cy[rib_pos_cy[rib_pos_index]]
                else:
                    pos_jump_cy= np.delete(pos_jump_cy, index)
                    pos_w_cy= np.delete(pos_w_cy, index)
            else:
                
                pos_jump_cy[index] = rib_pos_cy[rib_pos_index]
                pos_w_cy[index] = w_cy[rib_pos_cy[rib_pos_index]]
                self.a += w_cy[rib_pos_cy[rib_pos_index]]

            # Check if ribosome behing the moved ribosome can move in next step and update values
            if rib_pos_index != 0 and rib_pos_cy[rib_pos_index] - rib_pos_cy[rib_pos_index-1] == self.L+1 :
                
                pos_jump_cy = np.insert(pos_jump_cy,index,rib_pos_cy[rib_pos_index-1]) 
                pos_w_cy = np.insert(pos_w_cy,index,w_cy[rib_pos_cy[rib_pos_index-1]])
                self.a += w_cy[rib_pos_cy[rib_pos_index-1]]
            
                
            elif rib_pos_index == 0 and rib_pos_cy[rib_pos_index] == self.L+2:
                pos_jump_cy = np.insert(pos_jump_cy,0,0) 
                pos_w_cy = np.insert(pos_w_cy,0,w_cy[0])
                self.a += w_cy[0]
        
            move = True
    
 
        return t, enter, do_exit, move, [pos_jump_cy, pos_w_cy, rib_pos_cy,w_cy]



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
        
  
    cdef void __calcDensity(self,double tau, np.ndarray[np.int_t, ndim=1] rib_pos_cy,np.ndarray[np.double_t, ndim=1] rho_ ):

        """
        Update density profile after a simulation step

        """
        cdef Py_ssize_t i
        for i in rib_pos_cy:
            
            rho_[i-1] += tau
        
        return
    
    def getTot_den(self):

        """
        Calculate total density given the density profile

        """
        
        return np.sum(self.rho_)/len(self.rho_)



