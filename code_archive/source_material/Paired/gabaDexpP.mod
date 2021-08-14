NEURON {
	POINT_PROCESS gabaDexpP
	USEION ca READ eca
	USEION cal READ ical

	RANGE tau1a, tau2a, ea, i, Ga, Gb, igabaa, igabab, tau1b, tau2b, eb, ki, source, dest
	RANGE D, tauD, U, tau1g, tau2g, gaba, kg, initW
	RANGE eca, tauCa, Icatotal
	RANGE Cainf, pooldiam, z
	RANGE ICag, P0g, fCag
	RANGE fmax, fmin, wmax, wmin, w
	RANGE L2, L3
	NONSPECIFIC_CURRENT i

	RANGE ga, gb
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	pi = 3.141592 (1)
}

PARAMETER {
	tau1a=2 (ms) <1e-9,1e9>  : Cf Garden et al. 2002
	tau2a = 25 (ms) <1e-9,1e9>
	ea=-85	(mV)
	
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
	fmax = 3
	fmin = .2
	
	Cainf = 50e-6 (mM)
	pooldiam =  1.8172 (micrometer)
	z = 2
    L2 = 1    L3 = 0.05
	
	dth = 0.6
	pth = 0.8
	
	k = 0.01	
	
	tauCa = 50 (ms)
	
	P0g = .01
	fCag = .024
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
	wmin
	wmax
	ICag (mA)
	Afactor	(mM/ms/nA)
	Icatotal (mA)
	eca (mV)
	ical (nA)
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
	capoolcon
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
	
	wmax = fmax*initW
	wmin = fmin*initW
	
	capoolcon = Cainf
	Afactor	= 1/(z*FARADAY*4/3*pi*(pooldiam/2)^3)*(1e6)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gaba = (Bg - Ag)
	ga = (Ba - Aa)
	gb = (Bb - Ab)
	igabaa = w*Ga*ga*(v - ea)
	igabab = w*Gb*gb*(v - eb)
	
	i = igabaa + igabab
}
 
DERIVATIVE state { LOCAL caterm
	Aa' = -Aa/tau1a
	Ba' = -Ba/tau2a
	Ab' = -Ab/tau1b
	Bb' = -Bb/tau2b
	Ag' = -Ag/tau1g
	Bg' = -Bg/tau2g
	can'=0
	cam'=0
	
	w' = eta(capoolcon)*(L2*omega(capoolcon) - L3*(w-initW))
	
		:Weight value limits
	if (w > wmax) { 
		w = wmax
	} else if (w < wmin) {
 		w = wmin
	}
	
	ICag = P0g*i
	Icatotal = ICag + k*ical*4*pi*((15/2)^2)*(0.01)    :  icag+k*ical*Area of soma*unit change
	capoolcon'= -fCag*Afactor*Icatotal + (Cainf-capoolcon)/tauCa
	
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
FUNCTION omega(ca) {
                if(ca<dth) {
                        omega = 0
                } else if(ca<pth) {
                        omega = -sqrt((pth-dth)/2*(pth-dth)/2-(ca-(dth+pth)/2)*(ca-(dth+pth)/2))
                } else {
			omega = 1/(1+50*exp(-50*(ca-pth)))
                }
	
}

FUNCTION eta(Cani (mM)) {
	LOCAL taulearn, P1, P2, P4, Cacon
	P1 = 0.1
	P2 = P1*1e-4
	P4 = 1
	Cacon = Cani*1e3
	taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
	eta = 1/taulearn*0.001
}