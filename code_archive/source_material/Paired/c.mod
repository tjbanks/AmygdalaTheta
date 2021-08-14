TITLE fast calcium and voltage-dependent potassium channel IKC

: kinetics taken from LA neuron model (Li et al. 2009), which is in turn taken from a PFC model (Dursewitz et al. 2000)

NEURON {
	SUFFIX c
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

	LOCAL vm, alpha, beta, va

	UNITSOFF
	
	vm = v + logcafun(cai)  : Takes cai in mM and multiplies by 1e3 to convert to equivalent uM for this equation
	
	va = vm+18
	
	if(va<0) {va = -va}
	
	if(va<0.05) { vm = vm+0.1 }
	
	
	alpha = alphafun(vm)
	beta = betafun(vm)

	minf = alpha/(alpha+beta)
	if(1/(alpha+beta) > 1.1) {
	taum = 1/(alpha+beta)
	} else { taum = 1.1
	}
	
	
	UNITSON
}

FUNCTION logcafun(cai) {
	TABLE FROM 0.00005 TO 0.01 WITH 1000
	logcafun = 40*log10(1000*cai)
}

FUNCTION alphafun(vm(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	alphafun = (-0.00642*vm-0.1152)/(exp((vm+18)/-12)-1)
}

FUNCTION betafun(vm(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	betafun = 1.7*exp((vm+152)/-30)
}


