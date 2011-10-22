#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import ref, man


#----------------------------------------------------------------------------
def _createFields( runTime, mesh ):
    ref.ext_Info() << "Reading field T\n" << ref.nl

    T = man.volScalarField( man.IOobject( ref.word( "T" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    ref.ext_Info() << "Reading transportProperties\n" << ref.nl

    transportProperties = man.IOdictionary( man.IOobject( ref.word( "transportProperties" ),
                                                          ref.fileName( runTime.constant() ),
                                                          mesh,
                                                          ref.IOobject.MUST_READ,
                                                          ref.IOobject.NO_WRITE ) )

    ref.ext_Info() << "Reading diffusivity DT\n" << ref.nl

    DT = ref.dimensionedScalar( transportProperties.lookup( ref.word( "DT" ) ) )
        
    return T, transportProperties, DT


#--------------------------------------------------------------------------------------
def write( runTime, mesh, T ):
    if runTime.outputTime():
        gradT = ref.fvc.grad(T)
       
        gradTx = ref.volScalarField( ref.IOobject( ref.word( "gradTx" ), 
                                                   ref.fileName( runTime.timeName() ),
                                                   mesh,
                                                   ref.IOobject.NO_READ,
                                                   ref.IOobject.AUTO_WRITE ),
                                     gradT.component( ref.vector.X ) )
      
        gradTy = ref.volScalarField( ref.IOobject( ref.word( "gradTy" ),
                                                   ref.fileName( runTime.timeName() ),
                                                   mesh,
                                                   ref.IOobject.NO_READ,
                                                   ref.IOobject.AUTO_WRITE ),
                                     gradT.component( ref.vector.Y ) )

        gradTz = ref.volScalarField( ref.IOobject( ref.word( "gradTz" ),
                                           ref.fileName( runTime.timeName() ),
                                           mesh,
                                           ref.IOobject.NO_READ,
                                           ref.IOobject.AUTO_WRITE ),
                                 gradT.component( ref.vector.Z ) )
        runTime.write()
        pass
    

#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    T, transportProperties, DT = _createFields( runTime, mesh )
    
    simple = man.simpleControl( mesh )

    ref.ext_Info() << "\nCalculating temperature distribution\n" << ref.nl
    
    while runTime.loop() :
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        for nonOrth in range( simple.nNonOrthCorr() + 1 ):
            ref.solve( ref.fvm.ddt( T ) - ref.fvm.laplacian( DT, T ) )
            pass
        
        write( runTime, mesh, T )
        
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
else:
   from Foam.OpenFOAM import ext_Info
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher \n "
   pass


#--------------------------------------------------------------------------------------

