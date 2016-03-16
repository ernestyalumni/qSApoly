/* AeqB.gp
 * Implementation in Pari/GP of material in
 * Marko Petkovsek, Herbert S. Wilf, Doron Zeilberger. A = B. A K Peters, Ltd. 1997. 
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160315
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
 */
\\ cf. http://math.mit.edu/~brubaker/PARI/PARIrefcard.pdf
\\ To run this script, you can go into the gp calculator with the command
\\ gp
\\ at the command prompt, and then in the gp command prompt (?), type 
\\ \r AeqB.gp
\\ using the read file command \r and the filename qPochh.gp


/* Chapter 2. Tightening the Target */

\\ pp. 19 famous identity with binomial coefficients

coeffsbinom1(n) = sum(k=0,n,binomial(n,k))
coeffsbinom2(n) = sum(k=0,n,binomial(n,k)^2)

\\ 2.3 Human and computer proofs; an example pp. 23

\\ 2.4 A Mathematica session
\\ 2.5 A Maple session

factorial_sanity_check(n) = {
			  return factorial(n+1)/factorial(n) }

R(n,k) = -k^2*(3*n+3-2*k)/(2*(n+1-k)^2*(2*n+1) )
F(n,k) = factorial(n)^4/(factorial(k)^2*factorial(n-k)^2*factorial(2*n) )
G(n,k) = R(n,k)*F(n,k)
\\ F(n+1,k)-F(n,k) -G(n,k+1)+G(n,k)

\\ v2 is version 2 - using the binomial function in Pari/GP instead of factorial's
Fv2(n,k) = binomial(n,k)^2/binomial(2*n,n)
Gv2(n,k) = R(n,k)*Fv2(n,k)
\\ Fv2(n+1,k)-Fv2(n,k)-Gv2(n,k+1)+Gv2(n,k) 



stdWZprf(n,k,summand,rhs,R) =
		      {
			F(n,k) = summand(n,k)/rhs(n);
			return(F)
	 }

