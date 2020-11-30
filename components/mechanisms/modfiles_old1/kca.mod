
COMMENT

kca.mod

Calcium-dependent potassium channel
Based on
Pennefather (1990) -- sympathetic ganglion cells
taken from
Reuveni et al (1993) -- neocortical cells

Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
	
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE n, gk, gbar
	RANGE ninf, ntau
	GLOBAL Ra, Rb, caix
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 10   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
	cai  		(mM)
	caix = 1	
									
	Ra   = 0.01	(/ms)		: max act rate  
	Rb   = 0.02	(/ms)		: max deact rate 

	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		: original temp 	
	q10  = 2.3			: temperature sensitivity

	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau 		(ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	rates(cai)
	n = ninf
}

BREAKPOINT {
        SOLVE states
	gk = tadj*gbar*n
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp

PROCEDURE states() {   :Computes state variable n 
        rates(cai)      :             at the current v and dt.
        n = n + nexp*(ninf-n)

        VERBATIM
        //return 0;
        ENDVERBATIM
}

PROCEDURE rates(cai(mM)) {  

        LOCAL tinc

        a = Ra * cai^caix
        b = Rb
        ntau = 1/(a+b)
	ninf = a*ntau

        tadj = q10^((celsius - temp)/10)

        tinc = -dt * tadj
        nexp = 1 - exp(tinc/ntau)
}











