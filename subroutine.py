#contains some important subroutines to be run in background

import math

def sumup(RX,RY,RZ,box,n,rcut,rmin,sigma):
	sigmsq =sigma*sigma
	boxinv=1/box
	rcutsq=rcut*rcut
	rminsq=rmin*rmin
	v=0.0
	w=0.0
	
	for i in range (0,n-1):
		for j in range (i+1,n):
			rxij, ryij, rzij=RX[i]-RX[j],RY[i]-RY[j],RZ[i]-RZ[j]
	
			rxij=rxij-box*round(rxij*boxinv)
			ryij=ryij-box*round(ryij*boxinv)
			rzij=rzij-box*round(rzij*boxinv)

			rijsq=rxij**2+ryij**2+rzij**2

			if (rijsq<rminsq):
				print ("overlap",i,j,rijsq,rminsq)
				return
			elif(rijsq<rcutsq):
				sr2=sigmsq/rijsq
				sr6=sr2**3
				vij=sr6*(sr6-1.0)
				v=v+vij
				
	v=4.0*v
	return(v)


# THIS IS NAMED AS ENERGY WHICH CALCULATES THE ENERGY OF SINGLE PARTICLE
def energy(N,RXI,RYI,RZI,I,RCUT,SIGMA,box,RX,RY,RZ):
    RCUTSQ=RCUT**2
    SIGSQ=SIGMA*SIGMA
    boxinv=1/box
    V=0.0
    
    for J in range (0,N):
        if (J != I):
            rxij,ryij,rzij=RXI-RX[J],RYI-RY[J],RZI-RZ[J]
            
            rxij=rxij-box*round(rxij*boxinv)
            ryij=ryij-box*round(ryij*boxinv)
            rzij=rzij-box*round(rzij*boxinv)

            RIJSQ=rxij**2+ryij**2+rzij**2
            if (RIJSQ<RCUTSQ):
                SR2 = SIGSQ / RIJSQ
                SR6 = SR2 * SR2 * SR2
                VIJ = SR6 * ( SR6 - 1.0 )
                V   = V + VIJ

    V = 4.0 * V
    return(V)


def check(I,N,BOX,BOXINV,RXI,RYI,RZI,RMIN,RX,RY,RZ):
	NOVER=0.0
	for J in range (0,N):
		if (J!=I):
			RXIJ=RXI-RX[J]
			RYIJ=RYI-RY[J]
			RZIJ=RZI-RZ[J]

			RXIJ=RXIJ-BOX*round(RXIJ*BOXINV)
			RYIJ=RYIJ-BOX*round(RYIJ*BOXINV)
			RZIJ=RZIJ-BOX*round(RZIJ*BOXINV)

			DIS=math.sqrt(RXIJ**2+RYIJ**2+RZIJ**2)
			if (DIS>RMIN):
				NOVER=NOVER+1
	return(NOVER)


