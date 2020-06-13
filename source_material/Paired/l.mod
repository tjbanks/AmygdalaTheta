TITLE l-calcium channel  
: l-type calcium channel from Hemond et al. 2008, "Distinct classes of pyramidal cells exhibit mutually exclusive 
: firing patterns in hippocampal area CA3b"

: This file obtained from ModelDB: http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=101629

: This model was originally developed and explicitly described in Jaffe et al. 1994, "A Model For Dendritic Ca2+ Accumulation in Hippocampal
: Pyramidal Neurons Based on Fluorescence Imaging Measurements"

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER {
	v (mV)
	cai
	cao
	celsius 	(degC)
	g=.003 (mho/cm2)
	ki=.001 (mM)
	q10 = 5
	mmin=0.2
	tfa = 1
	a0m =0.1
	zetam = 2
	vhalfm = -25
	gmm=0.1	
	ggk
}


NEURON {
	SUFFIX cal
	USEION ca READ cai,cao
	USEION cal WRITE ical VALENCE 2
        RANGE G,cai, ica, g, ggk, i
        RANGE minf,tau
}

STATE {
	m
}

ASSIGNED {
	ical (mA/cm2)
        G (mho/cm2)
        minf
        tau   (ms)
		i
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	G = g*m*m*h2(cai)
	ggk=ghk(v,cai,cao)
	i = G*ggk
	ical = i

}

FUNCTION h2(cai(mM)) {
TABLE FROM 0 TO 0.01 WITH 1000
	h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alp(v(mV)) (1/ms) {
TABLE FROM -150 TO 150 WITH 500
	alp = 15.69*(-1.0*v+81.5)/(exp((-1.0*v+81.5)/10.0)-1.0)
}

FUNCTION bet(v(mV)) (1/ms) {
TABLE FROM -150 TO 150 WITH 500
	bet = 0.29*exp(-(v)/10.86)
}

FUNCTION alpmt(v(mV)) {
TABLE FROM -150 TO 150 WITH 500
  alpmt = exp(0.0378*zetam*(v-vhalfm)) 
}

FUNCTION betmt(v(mV)) {
TABLE FROM -150 TO 150 WITH 500
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm)) 
}

DERIVATIVE state {  
        rate(v)
        m' = (minf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a, b, qt
        qt=q10^((celsius-25)/10)
        a = alp(v)
        b = 1/((a + bet(v)))
        minf = a*b
	tau = betmt(v)/(qt*a0m*(1+alpmt(v)))
	if (tau<mmin/qt) {tau=mmin/qt}
}
