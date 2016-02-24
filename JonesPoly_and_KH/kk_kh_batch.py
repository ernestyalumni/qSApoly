## kk_kh_batch.py
## knotkit Khovanov Homology batch file
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160222
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
## (at your option) any later version.                                        
##     
## This program is distributed in the hope that it will be useful,               
## but WITHOUT ANY WARRANTY; without even the implied warranty of              
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 
## GNU General Public License for more details.                             
##                                                         
## linkedin     : ernestyalumni                                                    
## wordpress    : ernestyalumni                                                    
############################################################################ 
"""
kk_kh_batch.py

I wanted to use knotkit (kk) to compute the Khovanov Homology (kh) of various
Torus knots, and to save the output.  This Python script automates this 
process, creating a "batch" of (text) files that contain the Khovanov Homology
of Torus knots as polynomials.

I assume that kk_kh_batch.py is placed in the file directory "above" your 
local copy of knotkit, i.e. ./kk_kh_batch.py and ./knotkit
"""

import subprocess # cf. http://stackoverflow.com/questions/2473655/how-to-make-a-call-to-an-executable-from-python-script

from itertools import combinations

import fractions
from fractions import gcd # cf. http://stackoverflow.com/questions/11175131/code-for-greatest-common-divisor-in-python

# list of torus knots; go ahead and edit to add/subtract what torus knots
# you'd like
TKNOTSLIST = ["unknot", "T(2,3)", "T(2,5)","T(2,7)","T(2,9)","T(2,11)","T(2,13)","T(3,4)","T(3,5)"] # list of knots with less than 14 crossings
# the fact that these knots have less than 14 crossings can be checked by first making a list of all possible coprime pair of numbers
[i for i in combinations(range(1,14),2) if gcd(i[0],i[1])==1]
# and then checking each pair of numbers in SnapPy to read out the crossings:
# in SnapPy: T0305 = Link("T(3,5)")
# T0305 # <Link: 1 comp; 10 cross>


# go ahead and edit args to suit your needs

def direct_output(knotlist=TKNOTSLIST,filename="khTknotsless14"):
    """
    direct_output = direct_output(knotlist=TKNOTSLIST)
    """
    knotpolys = []
    for knot in knotlist:
        args = ("./knotkit/kk","khplainout","-v",knot)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        popen.wait()
        output = popen.stdout.read()
        knotpolys.append(output)
    f = open(filename,'w')
    for item in knotpolys:
        f.write(item)
    return knotpolys

if __name__ == "__main__":
    defaultTknotpolys_14 = direct_output()
    print "khTknotsless14 file, containing all categorified Jones polynomials from Khovanov homology for torus knots of less than 14 crossings, has been created (as a default)."
