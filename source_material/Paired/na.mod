: spike-generating sodium channel (Pyramid)

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE gnabar, gna
	RANGE minf, hinf, mtau, htau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar = 0.12 (siemens/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
	gna (siemens/cm2)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*m*m*m*h
	ina = gna*(v-ena)
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

FUNCTION malf(v (mV)) (/ms) {
	LOCAL va 
	UNITSOFF
	va = v+25
	:if (fabs(va)<1e-04){
	:	va = va+0.0001
	:}
	malf = -0.2816*va/((exp(-va/9.3))-1)
	UNITSON
}

FUNCTION mbet(v (mV)) (/ms) {
	LOCAL vb 
	UNITSOFF
	vb=v-2
	:if (fabs(vb)<1e-04){
	:   	vb = vb+0.0001
	:}
	mbet = 0.2464*vb/((exp(vb/6))-1)
	UNITSON
}

FUNCTION half(v (mV)) (/ms) {
	UNITSOFF
	half = 0.098*exp(-(v+40.1)/20)
	UNITSON
}

FUNCTION hbet(v (mV)) (/ms) {
	UNITSOFF
	hbet=1.4/(exp(-(v+10.1)/10)+1)
	UNITSON
}

PROCEDURE rate(v (mV)) {
	LOCAL msum, hsum, ma, mb, ha, hb
	UNITSOFF
	ma = malf(v) mb = mbet(v) ha = half(v) hb = hbet(v)
	
	msum = ma+mb
	minf = ma/msum
	mtau = 1/msum
	
	hsum = ha+hb
	hinf = ha/hsum
	htau = 1/hsum
	UNITSON
}
