"""
gene_rho_file_maker.py

"""

def MakeGeneFile(infile,low_codon_lim,up_codon_lim):
    """
    Make file containing genes with codon numbers within a certain bracket

    Inputs:
           infile <string> : Name of file containing all genes
           low_codon_lim <int> : Lower limit on codon number
           up_codon_lim <int> : Upper limit on codon number
    """
    #Open files
    filein = open(infile,'r')
    outfile_name = infile[:-4]+'_'+str(up_codon_lim)+'.dat'
    fileout = open(outfile_name,'w')
    
    lines = filein.readlines()
    
    for line in lines:
        
        data = line.split(',')
        # Write data in new file if number of codons falls within the bracket
        if int(data[2]) > low_codon_lim and int(data[2]) <= up_codon_lim:
            if '#' not in data[3]:
                fileout.write(line)
                
    #Close files
    filein.close()
    fileout.close()




