## groups.sage
## I implement groups, matrix groups, Lie groups, Lie Algebras in Sage Math
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

#######################################################
##### Chapter 2 Representations of Finite Groups of 
####  Kosmann-Schwarzbach (2010) 
#######################################################

##############################
## Section 1 Representations
##############################

#############################################
## Example 1.2. pp. 10 of 1.1 General Facts
#############################################

S_3 = SymmetricGroup(3)
S_3.list() # see the elements of $S_3$; there are 6 of them
t = S_3((2,3))
t.tuple() # clearly t does what we want, 123 \mapsto 132
c = S_3((1,2,3))
c.tuple() # clearly, c does what we want, 123 \mapsto 231
j = exp(2*I*pi/3)

tc_3 = S_3.subgroup((t,c))
rhot = Matrix(CC,[[0,1],[1,0]])
rhoc = Matrix(CC,[[j,0],[0,j**2]])
rhotc = MatrixGroup((rhot,rhoc))

# Matfin = sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_generic



S_4 = SymmetricGroup(4)
A_4 = AlternatingGroup(4)
D_2 = DihedralGroup(2)

################################################## 
## Exercise 2.15(c) of Kosmann-Schwarzbach (2010)
##################################################

P2CC.<x,y> = PolynomialRing(CC,2) # this declares a PolynomialRing of field of complex numbers, of order 2 (i.e. only 2 variables for a polynomial, such as x, y)

x**2 in P2CC # True, sanity check

A = var('A')
assume(A,''complex'')
B = var('B')
assume(B,''complex'')
C = var('C')
assume(C,''complex'')

f(x,y) = A*x**2 +2*B*x*y + C*y**2
P2CC.<x,y> = PolynomialRing(CC,2)
P2CC.gens() # (x, y) so x and y are truly the generators of this polynomial ring

x>y # True, lexicographic order
x**2 > x*y # True, lexicographic order
x*y > y**2 # True, lexicographic order

P2CCbasis = [ x**2,x*y,y**2]

def T(f,basis=P2CCbasis):
    coeffs = [f.coefficient(e_i) for e_i in basis]
    return coeffs

a = var('a')
assume(a,"complex")
b = var('b')
assume(b,"complex")
c = var('c')
assume(c,"complex")
d = var('d')
assume(d,"complex")
g = Matrix([[a,b],[c,d]] )
X = Matrix([[x],[y]])

f( (g.inverse()*X)[0,0], (g.inverse()*X)[1,0] ).expand().coefficient(x^2).full_simplify()
# (C*c^2 - 2*B*c*d + A*d^2)/(b^2*c^2 - 2*a*b*c*d + a^2*d^2)
f( (g.inverse()*X)[0,0], (g.inverse()*X)[1,0] ).expand().coefficient(x*y).full_simplify()
# -2*(C*a*c + A*b*d - (b*c + a*d)*B)/(b^2*c^2 - 2*a*b*c*d + a^2*d^2)
f( (g.inverse()*X)[0,0], (g.inverse()*X)[1,0] ).expand().coefficient(y^2).full_simplify()
# (C*a^2 - 2*B*a*b + A*b^2)/(b^2*c^2 - 2*a*b*c*d + a^2*d^2)

f((g*X)[0,0],(g*X)[1,0])
# (a*x + b*y)^2*A + 2*(a*x + b*y)*(c*x + d*y)*B + (c*x + d*y)^2*C
f((g*X)[0,0],(g*X)[1,0]).expand()
# A*a^2*x^2 + 2*B*a*c*x^2 + C*c^2*x^2 + 2*A*a*b*x*y + 2*B*b*c*x*y + 2*B*a*d*x*y + 2*C*c*d*x*y + A*b^2*y^2 + 2*B*b*d*y^2 + C*d^2*y^2
f((g*X)[0,0],(g*X)[1,0]).expand().coefficient(x^2)
# A*a^2 + 2*B*a*c + C*c^2
f((g*X)[0,0],(g*X)[1,0]).expand().coefficient(x*y)
# 2*A*a*b + 2*B*b*c + 2*B*a*d + 2*C*c*d
f((g*X)[0,0],(g*X)[1,0]).expand().coefficient(y^2)
# A*b^2 + 2*B*b*d + C*d^2

T( f((g*X)[0,0],(g*X)[1,0]).expand() )
#[A*a^2 + 2*B*a*c + C*c^2,
# 2*A*a*b + 2*B*b*c + 2*B*a*d + 2*C*c*d,
# A*b^2 + 2*B*b*d + C*d^2]

phi = var('phi')
assume(phi,"real")

########################################
### More on Exercise 2.15 G = SU(2)
########################################

U = exp(I*phi) * Matrix( [[a,b],[-b.conjugate(),a.conjugate()]])
U_try = Matrix( [[a,b],[-b.conjugate(),a.conjugate()]]) # trying without exp(I*phi)

f( ((exp(-I*phi)*U_try.adjoint())*X)[0,0],((exp(-I*phi)*U_try.adjoint())*X)[1,0] ).expand()
# C*a^2*y^2*e^(-2*I*phi) - 2*B*a*b*y^2*e^(-2*I*phi) + A*b^2*y^2*e^(-2*I*phi) + 2*B*a*x*y*conjugate(a)*e^(-2*I*phi) - 2*A*b*x*y*conjugate(a)*e^(-2*I*phi) + A*x^2*conjugate(a)^2*e^(-2*I*phi) + 2*C*a*x*y*conjugate(b)*e^(-2*I*phi) - 2*B*b*x*y*conjugate(b)*e^(-2*I*phi) + 2*B*x^2*conjugate(a)*conjugate(b)*e^(-2*I*phi) + C*x^2*conjugate(b)^2*e^(-2*I*phi)
f( ((exp(-I*phi)*U_try.adjoint())*X)[0,0],((exp(-I*phi)*U_try.adjoint())*X)[1,0] ).expand().coefficient(x^2)
# A*conjugate(a)^2*e^(-2*I*phi) + 2*B*conjugate(a)*conjugate(b)*e^(-2*I*phi) + C*conjugate(b)^2*e^(-2*I*phi)
f( ((exp(-I*phi)*U_try.adjoint())*X)[0,0],((exp(-I*phi)*U_try.adjoint())*X)[1,0] ).expand().coefficient(y^2)
# C*a^2*e^(-2*I*phi) - 2*B*a*b*e^(-2*I*phi) + A*b^2*e^(-2*I*phi)
f( ((exp(-I*phi)*U_try.adjoint())*X)[0,0],((exp(-I*phi)*U_try.adjoint())*X)[1,0] ).expand().coefficient(x*y)
# 2*B*a*conjugate(a)*e^(-2*I*phi) - 2*A*b*conjugate(a)*e^(-2*I*phi) + 2*C*a*conjugate(b)*e^(-2*I*phi) - 2*B*b*conjugate(b)*e^(-2*I*phi)


######################################################################
##### John Baez, Javier P Muniain.  Gauge Fields, Knots and Gravity
######################################################################

#####################################################################################
### Exercise 22 G = SU(2), fundamental representation and P^{(1)}(\mathbb{C}^2)
#####################################################################################

def f1(x):
    """
    f1 : \mathbb{C}^2 \to P^{(1)}(\mathbb{C}^2)
    f1 : (x,y) \mapsto Ax + By

    INPUT:
    (x,y) is a 2x1 Matrix of complex entries
    OUTPUT:
    Ax+By = A*x + B*y is a (homogeneous) polynomial of degree 1
    """
    return A*x[0,0] + B*x[1,0]
 
Paulimat = [Matrix([[1,0],[0,1]]),Matrix([[0,1],[1,0]]),Matrix([[0,-I],[I,0]]),Matrix([[1,0],[0,-1]])]

U_try1 = Matrix( [[a.conjugate(),-b],[b.conjugate(),a]] )

f1( U_try1*X).coefficient(x) # A*conjugate(a) + B*conjugate(b)
f1( U_try1*X).coefficient(y) # B*a - A*b

Paulimat[3] * Paulimat[1]*U_try*Paulimat[1] * Paulimat[3]
# [conjugate(a) conjugate(b)]
# [          -b            a]
