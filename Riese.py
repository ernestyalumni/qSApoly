## Riese.py
## This is my implementation in sympy of 
##   qMultiSum - A Package for Proving q-Hypergeometric Multiple Summation Identities
##   by Axel Riese
################################################################################
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>
## 20160107
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
################################################################################
import sympy

from sympy.abc import a,i,j,k,n,q
from sympy import Product, Symbol
from sympy import Rational as Rat

n_    = Symbol('n_',positive=True)
k_pos = Symbol('k_pos',positive=True)
k_neg = Symbol('k_neg',positive=False)

def qPochhammer(a,q,k):    
    try:
        if k == 0:
            return 1
        elif k<0:
            aq_kinv = Product(Rat(1)-a*q**(-i),(i,1,-k))
            return Rat(1)/aq_kinv
        elif k>0:
            aq_k = Product(Rat(1)-a*q**i,(i,0,k))
            return aq_k
    except TypeError:
        aq_k = Product(Rat(1)-a*q**i,(i,0,k))
        return aq_k

def qbinom(n,k):
    try:
        if n >= k and k >= 0:
            n_k_q = qPochhammer(q,q,n)/(qPochhammer(q,q,n-k)*qPochhammer(q,q,k))
            return n_k_q
        elif n<k or k<0:
            return 0
    except TypeError:
        n_k_q = qPochhammer(q,q,n)/(qPochhammer(q,q,n-k)*qPochhammer(q,q,k))
        return n_k_q
        
# q-Vandermonde identity

F__nk = (-Rat(1))**k*q**((n-k)**2)*(qbinom(2*n,k))**2



