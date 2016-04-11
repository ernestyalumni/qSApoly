## Presentations_Finitely_Presented_Groups.py
## This is my implementation of Presentations in Abstract Algebra
##   utilizing
## Sage Math
##  
## The main references that I'll liberally copy from is from are 
##   Advanced Modern Algebra
##   by J. Rotman
## and Abstract Algebra and Sage (aata)
################################################################################
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>
## 
## 20160411
## 
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, 
## or (at your option) any later version. 
## 
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
## GNU General Public License for more details.
## 
## linkedin     : ernestyalumni 
## twitter      : ernestyalumni 
## wordpress    : ernestyalumni 
##  
################################################################################

# Advanced Modern Algebra - Rotman
# cf. Rotman, Chapter 5, Groups II, Proposition 5.80
# generalized quaternion group Q_n, for n >= 3

def make_Q_n(n):
    assert n >= 3
#    F.<a,b> = FreeGroup()
    F = FreeGroup(2,'a, b')
    a,b = F.gens()
    G = F/[a**(2**(n-1)),a*b*a*b**(-1),a**(2**(n-2))*b**(-2)]
    return G

def make_D_2n(n):
    assert n > 0
    F = FreeGroup(2,'a, b')
    a,b = F.gens()
    G = F/[a**n,b**2,b*a*b*a]
    return G

