TITLE perirhinal a-type potassium current IA

: kinetics taken from study of A-current in PRC neurons (Biella et al. 2007)

NEURON {
	SUFFIX a
	USEION k READ ek WRITE ik
	RANGE g, G
	RANGE minf, taum, hinf, tauh, i
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	g = 0.012 (siemens/cm2) <0,1e9>
	sh1 = 0
	sh2 = 0
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	G (siemens/cm2)
	minf
	taum (ms)
	hinf
	tauh (ms)
	i
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	G = g*m*m*m*h
	i = G*(v-ek)
	ik=i
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/taum
	h' = (hinf-h)/tauh
}

PROCEDURE rate(v (mV)) {
	UNITSOFF
	
	minf = minffun(v)
	taum = taumfun(v)
	
	hinf = hinffun(v)
	tauh = tauhfun(v)
	
	UNITSON
}

FUNCTION minffun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	minffun = 1/(1+exp((v+46.9)/-14.1))
}

FUNCTION taumfun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	taumfun = sh1+0.3+0.2*exp(-v/27)
}

FUNCTION hinffun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	hinffun = 1/(1+exp((v+79.1)/9.3))
}

FUNCTION tauhfun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	tauhfun = 14+exp(v/30)
}