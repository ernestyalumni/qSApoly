## Lie.sage
## I implement matrix groups, Lie groups, Lie Algebras in Sage Math
## 
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
##                                                            
## 20150916
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

GL(2,RR)
GL(2,RR).matrix_space().gens()

GL(2,CC)

SL(2,RR)
SL(2,CC)

SO(3,RR)
SO(3,RR).invariant_bilinear_form()
SO(3,RR).invariant_quadratic_form()

SU(2,CC)
SU(2,CC).algebra(CC)
SU(2,CC).algebra(CC).zero()

#########################
##### Pauli matrices
#########################
Paulimat = [matrix(CC,[[1,0],[0,1]]),matrix(CC,[[0,1],[1,0]]),matrix(CC,[[0,-I],[I,0]]),matrix(CC,[[1,0],[0,-1]])]

SU2matrep = MatrixGroup(Paulimat)



