# Mutate residues in a protein using pymol
# Import the libraries
import sys
# For some reason the path won't save
sys.path.append('/opt//local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')

import __main__
__main__.pymol_argv = ['pymol','-c'] # Pymol: supress GUI

#import os 
#os.chdir('/Users/jonathanhuihui/Documents/resurrected_trx/pdbs')


from pymol import cmd
import pymol
pymol.finish_launching()



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

cmd.load(pdb)

cmd.wizard("mutagenesis")
cmd.refresh_wizard()
cmd.get_wizard().do_select(str(selection))
cmd.get_wizard().set_mode(str(mutant))
cmd.get_wizard().apply()
cmd.set_wizard()
cmd.save( str(name) + '_' + str(mut) + '.pdb', name)








