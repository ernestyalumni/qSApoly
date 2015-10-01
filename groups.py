## groups.py
## I implement groups, matrix groups, Lie groups, Lie Algebras using sympy for Python
## The main reference I use, besides my own qsuperApoly.pdf write up is 
##
## Yvette Kosmann-Schwarzbach, Groups and Symmetries: From Finite Groups to Lie Groups, 
## Springer, 2010. e-ISBN 978-0-387-78866-1 
## http://www.caam.rice.edu/~yad1/miscellaneous/References/Math/Groups/Groups\%20and\%20Symmetries.pdf
##
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
##                                                            
## 20150923
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
## You can have received a copy of the GNU General Public License             
## along with this program; if not, write to the Free Software Foundation, Inc.,  
## S1 Franklin Street, Fifth Floor, Boston, MA                      
## 02110-1301, USA                                                             
##                                                                  
## Governing the ethics of using this program, I default to the Caltech Honor Code: 
## ``No member of the Caltech community shall take unfair advantage of        
## any other member of the Caltech community.''                               
##                                                         
## Donate, and support my other scientific and engineering endeavors at 
## ernestyalumni.tilt.com                                                      
##                          
## Facebook     : ernestyalumni                                                   
## linkedin     : ernestyalumni                                                    
## Tilt/Open    : ernestyalumni                                                    
## twitter      : ernestyalumni                                                   
## youtube      : ernestyalumni                                                   
## wordpress    : ernestyalumni                                                    
##                                                                                  
############################################################################ 

import itertools
from itertools import product, permutations

import sympy
from sympy import I, LeviCivita
from sympy import Rational as Rat

from sympy.physics.matrices import msigma # <class 'sympy.matrices.dense.MutableDenseMatrix'>

#######################################################
##### Chapter 5 Representations of Finite Groups of 
####  Kosmann-Schwarzbach (2010) 
#######################################################

###########################################################################
## Bases of Lie Algebra su(2)\equiv su(2,\mathbb{C}) of Lie Group SU(2)
###########################################################################

def commute(A,B):
    """
    commute = commute(A,B)
    commute takes the commutator of A and B
    """
    return (A*B - B*A)

def xi(i):
    """
    xi = xi(i)
    xi is a function that returns the independent basis for 
    Lie algebra su(2)\equiv su(2,\mathbb{C}) of Lie group SU(2) of 
    traceless anti-Hermitian matrices, based on msigma of sympy
    cf. http://docs.sympy.org/dev/_modules/sympy/physics/matrices.html#msigma
    """
    if i not in [1,2,3]:
        raise IndexError("Invalid Pauli index")
    elif i==1:
        return I/Rat(2)*msigma(1)
    elif i==2:
        return -I/Rat(2)*msigma(2)
    elif i==3:
        return I/Rat(2)*msigma(3)

## check anti-Hermitian property and commutation relations with xi 
# xi is indeed anti-Hermitian
xi(1) == -xi(1).adjoint() # True
xi(2) == -xi(2).adjoint() # True
xi(3) == -xi(3).adjoint() # True

# xi obeys the commutation relations

for i,j in product([1,2,3],repeat=2): print i,j

for i,j in product([1,2,3],repeat=2): print i,j, "\t Commutator: ", commute(xi(i),xi(j))

## check traceless Hermitian property and commutation relations with Pauli matrices
# Pauli matrices i.e. msigam is indeed traceless Hermitian

msigma(1) == msigma(1).adjoint() # True
msigma(2) == msigma(2).adjoint() # True
msigma(3) == msigma(3).adjoint() # True

msigma(1).trace() == 0 # True
msigma(2).trace() == 0 # True
msigma(3).trace() == 0 # True

# Pauli matrices obey commutation relation
print "For Pauli matrices, the commutation relations are :\n"
for i,j in product([1,2,3],repeat=2): print i,j, "\t Commutator: ", commute(msigma(i),msigma(j))

for i,j,k in permutations([1,2,3],3): print "Commute: ", i,j,k, msigma(i), msigma(j),  ": and is 2*i of ", msigma(k), commute(msigma(i),msigma(j)) == 2*I*msigma(k)*LeviCivita(i,j,k)
