NEURON {
	POINT_PROCESS gabaDexp
	RANGE tau1a, tau2a, ea, i, Ga, Gb, igabaa, igabab, tau1b, tau2b, eb, ki, source, dest
	RANGE D, tauD, U, tau1g, tau2g, gaba, kg, initW
	NONSPECIFIC_CURRENT i

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
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gaba = (Bg - Ag)
	ga = (Ba - Aa)
	gb = (Bb - Ab)
	igabaa = Ga*ga*(v - ea)
	igabab = Gb*gb*(v - eb)
	
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
	
}
