
# coding: utf-8
#Imports
"""
Script.py

"""

import numpy as np
import math
from helper_methods_cy import *
from simulation_up_cy import *
from find_alpha import *
import sys
import os


def main(args):
    """
    Calculate alpha for gene names and corresponding experimental density read from a file
    and saves the calculated alpha in a new file.

    Inputs:
          args : Arguments in terminal. 
                 filein <string>: Name of file containing gene names and density
                 rates_path <string> : Path of the folder containing hopping rates for genes.
                 sweeps/rel_tol <string> <optional> : Specify whether to use sweeps or relative tolerance to determine steady state
                 sweeps/rel_tol <int>/<float> <optional> : Specify the number of sweeps or value of relative tolerance

    """
    #Create empty lists
    gene_names_ = []
    exp_den_ = []
    bracs_ = []
    sweeps = 0
    rel_tol = 1
    
    in_file_name = args[0]
    rates_path = args[1]

    #Open files
    filein = open(in_file_name,'r')
    fileout = open(in_file_name[:-4]+'_alpha.dat','a')
    fileout_check = open(in_file_name[:-4]+'_alpha.dat','r') #Check if a file has already been made and if yes, how many genes have been already calculated

    #read lines
    lines_check = fileout_check.readlines()
    
    if len(args) > 2:
        
        if args[2] == 'Sweeps' or args[2] == 'sweeps':
            
            sweeps = int(args[3])
            
        elif args[2] == 'Tolrance' or args[2] == 'tolrance':
            
            rel_tol = float(args[3])
        
        
    lines= filein.readlines()
    #Close files
    filein.close()
    fileout_check.close()
    
    
    line_no = 0
    for line in lines:
        
        line_no += 1
        
        if line_no > len(lines_check):

            #Read data from gene file
            data = line.split(',')
            
            gene_names_.append(data[0])
            exp_den_.append(float(data[3]))

            # Calculate bracs to use for brent
            bracs_.append([1e-4,1])
        
            if '#' not in data[4]:
                exp_alpha = float(data[4])
            
            else: 
                exp_alpha = -1
            
            if exp_alpha != -1 and exp_alpha > 0.2:
                bracs_[-1][0] = exp_alpha-0.2
        
            if exp_alpha != -1:
                bracs_[-1][1] = exp_alpha+0.2
           
            
    n = 0
    #Calculate alpha for the genes
    for (gene,bracs,exp_den) in zip(gene_names_,bracs_,exp_den_):

        # Print progress
        sys.stdout.write("Genes Computed: "+str(n+len(lines_check))+'\n')
        sys.stdout.flush()  
        
        rates_file = rates_path+'/'+gene+'_rates.dat'
        
        sim = Simulation(10,sweeps = sweeps, rel_tol= rel_tol)
        sim.ExtractRates(rates_file,0.1,18.03) # Extract hopping rates
        
        
        try:
            a = calcAlpha((exp_den,sim),bracs)
        
        except(ValueError):
            # If root not found by brent's, report error
            print("\n Wrong Brackets for gene number:  "+str(n+len(lines_check))+'\n')
            a = ['N/A']

        # Write Calculated alpha in a file
        fileout.write(gene+','+str(a[0])+str('\n'))
        
        n += 1
        

    #Close file
    fileout.close()
        
    

# If script is run, then execute the main method

if __name__ == '__main__':
    main(sys.argv[1:])

