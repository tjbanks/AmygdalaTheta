TITLE hyperpolarization-activated cation current IH

: Kinetics based off of rat LA model (Li et al. 2009), with half-activation and k slope
: parameter taken from a study of IH in rat entorhinal cortex (Richter et al. 2000)

NEURON {
	SUFFIX h2
	USEION h READ eh WRITE ih
	VALENCE 1
	RANGE G, g
	RANGE minf, taum, minfs, taums, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	g = 0.012 (siemens/cm2) <0,1e9>
	delM = 0
	epsiM = 0
	zetaM = 1
}

ASSIGNED {
	v (mV)
	eh (mV)
	i
	ih (mA/cm2)
	G (siemens/cm2)
	minf
	taum (ms)
	minfs
	taums (ms)
}

STATE {
	m ms
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	G = g*(m+ms)
	i = G*(v-eh)
	ih = i
}

INITIAL {
	rate(v)
	m = minf
	ms = minfs
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/taum
	ms' = (minfs-ms)/taums
}


PROCEDURE rate(v (mV)) {

	UNITSOFF
		
	minf = minffun(v)
	taum = taumfun(v)
	
	minfs = minfsfun(v)
	taums = taumsfun(v)
	
	UNITSON
}

FUNCTION minffun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	minffun = 1/(1+exp((v+74.2)/9.78))^1.36
}

FUNCTION taumfun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	taumfun = 0.51/(exp((v-1.7)/10)+exp((v+340)/-52))
}

FUNCTION minfsfun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	minfsfun = 1/(1+exp((v+72.83)/15.9))^58.5	
}

FUNCTION taumsfun(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	taumsfun = 5.6/(exp((v-17)/14)+exp((v+260)/-43))
}
