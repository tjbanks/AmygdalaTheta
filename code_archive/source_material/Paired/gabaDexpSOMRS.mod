NEURON {
	POINT_PROCESS gabaDexpSOMRS
	RANGE tau1a, tau2a, ea, i, Ga, Gb, igabaa, igabab, tau1b, tau2b, eb, ki, source, dest
	RANGE D, tauD, U, tau1g, tau2g, gaba, kg, initW
	NONSPECIFIC_CURRENT i
	RANGE F, f, tauF, D1, d1, tauD1, D2, d2, tauD2
	RANGE facfactor

	RANGE ga, gb
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1a=2 (ms) <1e-9,1e9>  : Cf Garden et al. 2002
	tau2a = 25 (ms) <1e-9,1e9>
	ea=-75	(mV)
	
	tau1b= 50 (ms) <1e-9,1e9>  :25
	tau2b = 500 (ms) <1e-9,1e9> :50
	eb=-92	(mV)
	tauD = 250 :250
	
	tau1g = 175
	tau2g = 200
	
	Ga = 0.1
	Gb = 0.002
	source
	dest
	U= 0.0 :0.4
	
	kg = 0.25
	
	initW = 1
	
	facfactor = 1
	: the (1) is needed for the range limits to be effective
        f = 1.2 (1) < 0, 1e9 >    : facilitation
        tauF = 20 (ms) < 1e-9, 1e9 >
        d1 = 0.5 (1) < 0, 1 >     : fast depression
        tauD1 = 25 (ms) < 1e-9, 1e9 >
        d2 = 0.5 (1) < 0, 1 >     : slow depression
        tauD2 = 20 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	v (mV)
	i (nA)
	ga (uS)
	factora
	gb (uS)
	factorb
	igabaa (nA)
	igabab (nA)
	D
	told
	gaba
	factorg
	tsyn
	F
	D1
	D2
}

STATE {
	Aa (uS)
	Ba (uS)
	Ab (uS)
	Bb (uS)
	Ag
	Bg
	w
	can
	cam
}

INITIAL {
	LOCAL tpa, tpb, tpg
	if (tau1a/tau2a > .9999) {
		tau1a = .9999*tau2a
	}
	if (tau1b/tau2b > .9999) {
		tau1b = .9999*tau2b
	}
	
		
	if (tau1g/tau2g > .9999) {
		tau1g = .9999*tau2g
	}
	
	Aa = 0
	Ba = 0
	Ab = 0
	Bb = 0
	Ag = 0
	Bg = 0
	
	tpa = (tau1a*tau2a)/(tau2a - tau1a) * log(tau2a/tau1a)
	factora = -exp(-tpa/tau1a) + exp(-tpa/tau2a)
	factora = 1/factora
	
	tpb = (tau1b*tau2b)/(tau2b - tau1b) * log(tau2b/tau1b)
	factorb = -exp(-tpb/tau1b) + exp(-tpb/tau2b)
	factorb = 1/factorb
	
	tpg = (tau1g*tau2g)/(tau2g - tau1g) * log(tau2g/tau1g)
	factorg = -exp(-tpg/tau1g) + exp(-tpg/tau2g)
	factorg = 1/factorg
	
	igabaa = 0
	igabab = 0
	w=initW
	can = 5e-5
	cam = 5e-5
	D=1
	told = -1e30
	F = 1
	D1 = 1
	D2 = 1
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gaba = (Bg - Ag)
	ga = (Ba - Aa)
	gb = (Bb - Ab)
	igabaa = Ga*ga*(v - ea)*facfactor
	igabab = Gb*gb*(v - eb)*facfactor
	
	i = igabaa + igabab
}
 
DERIVATIVE state { LOCAL caterm
	Aa' = -Aa/tau1a
	Ba' = -Ba/tau2a
	Ab' = -Ab/tau1b
	Bb' = -Bb/tau2b
	Ag' = -Ag/tau1g
	Bg' = -Bg/tau2g
	w'=0
	can'=0
	cam'=0
	
}

NET_RECEIVE(weight (uS)) {

	D = 1-(1-D*(1-U))*exp((t-told)/-tauD)
	
	Aa = Aa + w*factora*D*(kg/(gaba+kg))
	Ba = Ba + w*factora*D*(kg/(gaba+kg))
	
	Ab = Ab + w*factorb*D*(kg/(gaba+kg))
	Bb = Bb + w*factorb*D*(kg/(gaba+kg))
	
	Ag = Ag + w*factorg*D*(kg/(gaba+kg))
	Bg = Bg + w*factorg*D*(kg/(gaba+kg))

	told = t
	
	F  = 1 + (F-1)* exp(-(t - tsyn)/tauF)
	:D1 = 1 - (1-D1)*exp(-(t - tsyn)/tauD1)
	:D2 = 1 - (1-D2)*exp(-(t - tsyn)/tauD2)
	
	tsyn = t
	
	facfactor = F * D1 * D2
	:printf("%g\n",facfactor)
	F = F*f  :F+f  :F * f
	
	if (F > 3) { 
	F=3	}
	:if (facfactor < 0.5) { 
	:facfactor=0.5
	:}	
	:D1 = D1 * d1
	:D2 = D2 * d2
	
}
