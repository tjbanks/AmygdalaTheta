TITLE CaT channel alpha-1G from McRory et al, 2001
: Reversal potential described by Nernst equation
: M.Migliore Jan 2003

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v (mV)
	celsius		(degC)
	gbar=0.008 (mho/cm2) :0.008
    vhalfn=-45   (mV):-51.73:-36
    vhalfl=-74   (mV):-85.4:-74
    kn=6.53   (1)
    kl=-5.4   (1)
	q10=2.3
	cai 	= .00005 (mM)	: initial [Ca]i = 50 nM
	cao 	= 2	(mM)	: [Ca]o = 2 mM
	eca
}


NEURON {
	SUFFIX cav31
	USEION ca READ eca WRITE ica
        RANGE gbar, carev
        GLOBAL ninf,linf,taul,taun, q10
}

STATE {
	n
        l
}

ASSIGNED {
	ica	(mA/cm2)		: current
	carev	(mV)			: rev potential
        ninf
        linf      
        taul
        taun
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gbar*n*l*(v-carev)
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-22)/10)
        ninf = 1/(1 + exp(-(v-vhalfn)/kn))
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
        taun = (0.5+0.124*exp(-v/15.8))/qt
        taul = (10.4+0.0118*exp(-v/7.85))/qt
}




















