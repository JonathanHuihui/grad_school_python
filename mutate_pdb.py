# Mutate residues in a protein using pymol


# USAGE:
# pymol -c mutate_pdb.py -- <PDBFILE> <NAME> <RESIDUE NUMBER> <TARGET MUTATION>


# Import the libraries
import sys
# For some reason the path won't save
sys.path.append('/opt//local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')

import __main__
__main__.pymol_argv = ['pymol','-c'] # Pymol: supress GUI

#import os 
#os.chdir('/Users/jonathanhuihui/Documents/resurrected_trx/pdbs')


# Import pymol and initialize
from pymol import cmd
import pymol
pymol.finish_launching()


# Take values from command line
pdb = sys.argv[1]
name = sys.argv[2]
selection = sys.argv[3]
mutant = sys.argv[4]


mut = str(mutant)[0:3] + str(selection)[0]

print "pdb: " + pdb
print "name: " + name
print "selection: " + selection
print "mutant: " + mutant

pymol.cmd.reinitialize()


# Run the wizard and save
cmd.load(pdb)

cmd.wizard("mutagenesis")
cmd.refresh_wizard()
cmd.get_wizard().do_select(str(selection))
cmd.get_wizard().set_mode(str(mutant))
cmd.get_wizard().apply()
cmd.set_wizard()
# Save as <Protein Name>_<Target Residue><Residue Index>.pdb
cmd.save( str(name) + '_' + str(mut) + '.pdb', name)








