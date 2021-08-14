TITLE slow calcium-dependent potassium channel IAHP

: kinetics taken from hippocampal CA3 neuron model (Migliore et al. 1995)

NEURON {
	SUFFIX ahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE G, g
	RANGE minf, taum, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	g = 0.012 (siemens/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	cai
	ek (mV)
	i
	ik (mA/cm2)
	G (siemens/cm2)
	minf
	taum (ms)
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	G = g*m*m
	i = G*(v-ek)
	ik = i
}

INITIAL {
	rate(v)
	m = minf
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/taum
}


PROCEDURE rate(v (mV)) {

	LOCAL alpha, beta

	UNITSOFF
	
	alpha = 1.3e13*cai^4

	
	beta = 0.005

	minf = alpha/(alpha+beta)
	taum = 1/(alpha+beta)
	
	
	UNITSON
}
