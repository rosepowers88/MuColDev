#####################################
#
# simple script to create lcio files with single particle
# events - modify as needed
# @author F.Gaede, DESY
# @date 1/07/2014
#
# initialize environment:
#  export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
#
#####################################
import math
import random
from array import array
from g4units import deg

# --- LCIO dependencies ---
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

#---- number of events per momentum bin -----
nevt = 200

outfile = "LATEST.slcio"

#--------------------------------------------


wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )

wrt.open( outfile , EVENT.LCIO.WRITE_NEW ) 

random.seed()


#========== particle properties ===================

#momenta = [ 1. , 3., 5., 10., 15., 25., 50., 100. ]
#momenta = [ 5. ]

#generate flat in pT
pT = []
for i in range(500):
    pT.append(random.random()*300.+20)

genstat  = 1
pdg = 15
mass =  1.77686 #tau mass
charge = -1.

decayLen = 1.e22 
#=================================================



for pt in pT:
    
    for j in range( 0, nevt ):

        col = IMPL.LCCollectionVec( EVENT.LCIO.MCPARTICLE ) 
        evt = IMPL.LCEventImpl() 

        evt.setEventNumber( j )

        evt.addCollection( col , "MCParticle" )

        phi =  random.random() * math.pi * 2. #flat in phi
        theta = math.acos(2*random.random()-1) #flat in eta

        #energy   = math.sqrt( mass*mass  + p * p ) 
                
        px = pt * math.cos( phi )
        py = pt * math.sin( phi )
        pz = pt / math.tan( theta ) 

        momentum  = array('f',[ px, py, pz ] )  

        epx = decayLen * math.cos( phi ) * math.sin( theta ) 
        epy = decayLen * math.sin( phi ) * math.sin( theta )
        epz = decayLen * math.cos( theta ) 

        endpoint = array('d',[ epx, epy, epz ] )  
        

#--------------- create MCParticle -------------------
        
        mcp = IMPL.MCParticleImpl() 

        mcp.setGeneratorStatus( genstat ) 
        mcp.setMass( mass )
        mcp.setPDG( pdg ) 
        mcp.setMomentum( momentum )
        mcp.setCharge( charge ) 

        if( decayLen < 1.e9 ) :   # arbitrary ...
            mcp.setEndpoint( endpoint ) 


#-------------------------------------------------------

      

#-------------------------------------------------------


        col.addElement( mcp )
        
        wrt.writeEvent( evt )


wrt.close() 


#
#  longer format: - use ".hepevt"
#

#
#    int ISTHEP;   // status code
#    int IDHEP;    // PDG code
#    int JMOHEP1;  // first mother
#    int JMOHEP2;  // last mother
#    int JDAHEP1;  // first daughter
#    int JDAHEP2;  // last daughter
#    double PHEP1; // px in GeV/c
#    double PHEP2; // py in GeV/c
#    double PHEP3; // pz in GeV/c
#    double PHEP4; // energy in GeV
#    double PHEP5; // mass in GeV/c**2
#    double VHEP1; // x vertex position in mm
#    double VHEP2; // y vertex position in mm
#    double VHEP3; // z vertex position in mm
#    double VHEP4; // production time in mm/c
#
#    inputFile >> ISTHEP >> IDHEP 
#    >> JMOHEP1 >> JMOHEP2
#    >> JDAHEP1 >> JDAHEP2
#    >> PHEP1 >> PHEP2 >> PHEP3 
#    >> PHEP4 >> PHEP5
#    >> VHEP1 >> VHEP2 >> VHEP3
#    >> VHEP4;
