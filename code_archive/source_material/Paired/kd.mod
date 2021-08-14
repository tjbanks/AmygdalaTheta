TITLE Borg-Graham type generic K-DR channel IKD

: changes to alpha_n and beta_n made to match voltage-clamp data given in Biella et al. 2007
: This was done by substituting the alpha & beta functions from Li et al. 2009
: and multipling both functions by 0.17 to reproduce the approximate time of activation for non-A current shown in Biella et al. 2007


NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
        RANGE G,g, minf,hinf,taum,tauh,i
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
        ek		 (mV)
	celsius		(degC)
	g=.003 (mho/cm2)
        vhalfn=-32   (mV)
        vhalfl=-61   (mV)
        a0l=0.001      (/ms)
        a0n=0.03      (/ms)
        zetan=-5    (1)
        zetal=2    (1)
        gmn=0.4   (1)
        gml=0.88   (1)   : 1.0
		sh1 = 15
		sh2 = 15
}


STATE {
	m
    h
}

ASSIGNED {
	i
	ik (mA/cm2)
        minf
        hinf      
        G
        taum
        tauh
}

INITIAL {
        rates(v)
        m=minf
        h=hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	G = g*m^3*h
	i = G*(v-ek)
	ik = i

}

FUNCTION alpn(v(mV)) {
TABLE FROM -150 TO 150 WITH 500
:  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+32)))
  alpn = -0.17*0.036*(v-13+sh1)/(exp((v-13+sh1)/-20)-1)
}

FUNCTION betn(v(mV)) {
TABLE FROM -150 TO 150 WITH 500
:  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+32))) 
   betn = 0.17*0.0108*(v-23+sh2)/(exp((v-23+sh2)/12)-1) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+32))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+32))) 
}

DERIVATIVE states {  
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=3^((32-30)/10)
        a = alpn(v)
        minf = a/(a+betn(v))
        taum = 1/(a+betn(v))
        :a = alpl(v)
        hinf = 1    : 1/(1+a)
        tauh = 1000 : betl(v)/(q10*a0l*(1 + a))
}

