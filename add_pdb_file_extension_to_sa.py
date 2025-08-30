# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 10:39:18 2023

@author: tfarg
"""

import os
files = []
pdbs = []
os.chdir("Y:\\CHEM_Zhang_Students\\tfarg\\Modelling\\TEMPLATES_FOR_AMBER\\FF14IDPSFF_TIP3P\\HYPER_8-7-24")

for i in os.listdir():
    if i.find(".sa")==-1:
        pass
    elif i.find(".sa.viols")>0:
        pass
    else:
        files.append(i)
        pdbs.append("%s.pdb"%i)
        os.rename("%s"%i, "%s.pdb"%i)
b = " ".join(pdbs)
b
#to open all files, go to vmd and paste "vmd %s"%b into command line
        
