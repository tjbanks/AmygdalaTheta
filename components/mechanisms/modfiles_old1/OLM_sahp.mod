:  iC   fast Ca2+/V-dependent K+ channel

NEURON {
	SUFFIX sAHPOLM
	USEION k READ ek WRITE ik
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
        RANGE ik, gk, gsAHPbar
}

UNITS {
        (mM) = (milli/liter)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gsAHPbar= 0.0001	(mho/cm2)
	cai (mM)
}

ASSIGNED {
	v (mV)
	ek (mV)
	lcai (mM)
	tcai (mM)
	ik (mA/cm2)
	cinf 
	ctau (ms)
	gk (mho/cm2)
}

STATE {
	c
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gsAHPbar*c       
	ik = gk*(v-ek)
}



INITIAL {
	cai = lcai + tcai
	rate(v,cai)
	c = cinf
}

DERIVATIVE states {
		cai = lcai + tcai
        rate(v,cai)
	c' = (cinf-c)/ctau
}

UNITSOFF


FUNCTION calf(v (mV), cai (mM)) (/ms) { LOCAL vs, va
	UNITSOFF
	vs=10*log10(1000*cai)
	calf = 0.0048/exp(-0.5*(vs-35))
	UNITSON
}

FUNCTION cbet(v (mV), cai (mM))(/ms) { LOCAL vs, vb 
	UNITSOFF
	  vs=10*log10(1000*cai)
	  cbet = 0.012/exp(0.2*(vs+100))
	UNITSON
}

UNITSON

PROCEDURE rate(v (mV), cai (mM)) {LOCAL  csum, ca, cb
	UNITSOFF
	ca=calf(v, cai) 
	cb=cbet(v, cai)		
	csum = ca+cb
	cinf = ca/csum
	ctau = 48
	UNITSON
}
