# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array
from pyfem.fem.Assembly import assembleMassMatrix, assembleTangentStiffness

class ModalSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol     = 1.0e-3
    self.iterMax = 10 

    BaseModule.__init__( self , props )

    self.fext  = zeros( len(globdat.dofs) )  
    self.fext[23] = 1000.

  def run( self , props , globdat ):

    globdat.cycle = 1
    
    K,fint  = assembleTangentStiffness( props, globdat )
    M,Mlump = assembleMassMatrix( props , globdat )
          
    globdat.vecs , globdat.vals = globdat.dofs.eigensolve( K, M )

    globdat.active = False 
