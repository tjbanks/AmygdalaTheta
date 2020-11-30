: spike-generating sodium channel (interneuron)

NEURON {
	SUFFIX nainter_WB
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
	RANGE minf, hinf, mtau, htau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	dt (ms)
	gnabar = 0.035 (mho/cm2) <0,1e9>
	ena = 55.0 (mV)
}

STATE {
	m h
}

ASSIGNED {
	ina (mA/cm2)
	minf 
	hinf 
	mtau (ms)
	htau (ms)
	gna (mho/cm2)
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*m*m*m*h
	ina = gna*(v-ena)
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

UNITSOFF

FUNCTION malf(v(mV)) {
	malf = -1000.0*0.1*(v+35.0)/(exp(-0.1*(v+35.0))-1.0)
}

FUNCTION mbet(v(mV)) {
	mbet = 1000.0*4.0*exp(-(v+60.0)/18.0)
}	

FUNCTION half(v(mV)) {
	half = 5.0*0.07*exp(-(v+58.0)/20.0)
}

FUNCTION hbet(v(mV)) {
	hbet = 5.0*1.0/(exp(-0.1*(v+28.0))+1.0)
}

PROCEDURE rate(v(mV)) { LOCAL msum, hsum, ma, mb, ha, hb

	ma = malf(v) 
	mb = mbet(v) 
	ha = half(v) 
	hb = hbet(v)
	
	msum = ma+mb
	minf = ma/msum
	mtau = 1.0/(msum)
	
	hsum = ha+hb
	hinf = ha/hsum
	htau = 1.0/(hsum)
}

UNITSON

