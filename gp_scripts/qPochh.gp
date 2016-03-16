/* qPochh.gp
 * Implementation of the qPochhammer symbol in Pari/GP
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
\\ \r qPochh.gp
\\ using the read file command \r and the filename qPochh.gp

\\ sanity check
Pi

qPochh(a,q,k)
