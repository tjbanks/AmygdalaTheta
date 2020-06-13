
NEURON {
	POINT_PROCESS glu
	USEION ca READ eca
	USEION cal READ ical
	
	RANGE iampa, Ga, tau1a, tau2a, ea
	RANGE inmda, Gn, tau1n, tau2n, en
	RANGE tau1m2, tau2m2
	RANGE tau1m5, tau2m5
	RANGE factora, factorn, factorm2, factorm5, initW, fmax, fmin, wmax, wmin, w
	RANGE Vmax, Kl
	RANGE L1, L2, L3
	RANGE Vm5, Vm2
	RANGE fm, fn
	RANGE source, dest
	RANGE tauD, U
	RANGE Qip3, Minf, Hinf, TauM, TauH, IP3
	
	NONSPECIFIC_CURRENT i
}


PARAMETER {
	Gn = 0.001  :0.005 :0.01 
	Ga = 0.135  :0.16 :0.1   :Gn = 0.01
	Vm5 = 0.25 
	Vm2 = 1    :Vm5 = 0.25
	
	tau1a = 0.5 tau2a = 7
	tau1n = 8 tau2n = 125  :60
	tau1m5 = 100 tau2m5 = 400   :250, 2500
	tau1m2 = 100 tau2m2 = 200
	
	ea = 0 en = 0
	
    a = 0.00412      fcytratio = 0.00280899   :a = 0.05
	fcytbeta = 0.0095   ferratio = 0.015873   :fcytbeta = 1  :fcytbeta original .0035
	:fcytbeta= (.0088 original) (.0095 works!)   :ferratio = 0.015873					
	Vn = 0.04   Kn = 1.0  : (uM/ms), (uM)      :Vn = 0.04
	Vp = 0.04   Kp = 0.1  : (uM/ms), (uM)
	Vs = 13    Ks = 0.2  : (uM/ms), (uM)      :Vs = 1.3
	Vmax = 0.1  Kl = 0.1  : (uM/ms), (uM)
	Vsoc = 0.03 Ksoc = 100  : (uM/ms), (uM)     :Vsoc = 0.01  Ksoc = 1
		
	fn = 0.22
	fm = 6e5                      :fm = 20  //fm original value 6e5
						: fn 0.15
	Vr = 15   Vl = 0.0021  :(pL/ms), (pL/ms)
	
	bip3 = 4.1e-3     aip3 = 0.42          : (1/ms), (uM-ms)^-1
	ainh = 0.005       dinh = 0.525          : (uM-ms)^-1, (uM)    :ainh = 2e-4  dinh = 1.05
	dip3 = 0.13       ddis = 0.47          : (uM), (uM)            :ddis = 0.94
	dact = 0.6                          : (uM)                   :dact = 8.2e-2
	
	dth = 1  :1.6    
	pth = 1.2  	:1.75   :(uM), (uM) 0.3 0.6  0.8 may be good value for duh
	L1 = 0.5       L2 = 0.097    L3 = 0.05
			:L2 0.025
	tauD = 50   U = 0.4
	
	initW = 1
	fmax = 3
	fmin = .2
}

ASSIGNED {
	Qip3 Hinf Minf TauM TauH

	Jnmda	Jpmcan	Jncxn	Jleak
	
	factora factorn factorm2 factorm5
	i ga gn
	ical iampa inmda eca
	nwt 
	v
	w
	IP3
	
	Jip3r 	Jvgcc 	Jpmcam	Jncxm 	Jserca	Jsoc
	Jin Jout

	source dest
	D
	tlast
	
	wmin
	wmax
}

STATE {
	Aa Am5
	Ba Bm5
	An Am2
	Bn Bm2
	
	
	wn
	wm	
	can cam cae
	m h
}

INITIAL {
	Aa = 0  Ba = 0  An = 0  Bn = 0
	Am5 = 0 Bm5 = 0 Am2 = 0 Bm2 = 0
	
	calcfactors()
	w = initW
	wn = initW/2
	wm = initW/2
	wmax = fmax*initW
	wmin = fmin*initW
	Jnmda = 0
	can = 0.05   cam = 0.05  cae = 120   : (uM), (uM), (uM)   :cae = 4
	m = minf(IP3,cam)  h = hinf(IP3,cam)
	D = 1
	tlast = -1e30
}

BREAKPOINT {

	SOLVE state METHOD cnexp
	
	ga = (Ba - Aa)
	iampa = Ga*ga*w*(v - ea)

	gn = (Bn - An)
	inmda = Gn*gn*s(v)*(v - en)
	
	i = iampa + inmda	
}

DERIVATIVE state { 
	An' = -An/tau1n      Bn' = -Bn/tau2n
	Aa' = -Aa/tau1a      Ba' = -Ba/tau2a
	Am5' = -Am5/tau1m5   Bm5' = -Bm5/tau2m5
	Am2' = -Am2/tau1m2   Bm2' = -Bm2/tau2m2
	
	
	
:NMDA Ca2+ pool	
	if(nwt!=0) { Jnmda = -10*fn/nwt*(Gn*gn*s(v)*(v-eca)) }
	
	Jpmcan = Vp*can^2/(Kp^2+can^2)
	Jncxn = Vn*can^4/(Kn^4+can^4)
	Jleak = Vmax*(0.05)^4/(0.05^4+Kl^4)
	
	can' = Jnmda + Jleak - Jpmcan - Jncxn

:mGluR5 Ca2+ pool
    IP3 = Vm5*(Bm5 - Am5)
	
	m' = (minf(IP3,cam)-m)/taum(cam)
	h' = (hinf(IP3,cam)-h)/tauh(IP3,cam)
	
	Jip3r = (Vr*m^3*h^3 + Vl) *(cae-cam)
	Jvgcc = -fm*ical
	Jpmcam = Vp*cam^2/(Kp^2+cam^2)
	Jncxm = Vn*cam^4/(Kn^4+cam^4)	
	Jserca = Vs*cam^2/(Ks^2+cam^2)
	Jsoc =  -Vsoc*Ksoc^4/(Ksoc^4+cae^4)*(v-eca)
	
	Jin = a*(Jvgcc + Jsoc)
	Jout = Jpmcam + Jncxm
	
	cam' = fcytbeta*(Jin - Jout) + fcytratio*(Jip3r - Jserca)
	cae' = ferratio*(Jserca - Jip3r)
	
	:wn' = eta(can)*(L1*omegaN(can) - L3*(w-initW))  
	wn' = eta(can)*(L1*omegaN(can) - L3*(w-initW))  
	wm' = eta(cam/4)*(L2*omegaM(cam) - L3*(w-initW))
	
	w = wn+wm
	
		:Weight value limits
	if (w > wmax) { 
		w = wmax
	} else if (w < wmin) {
 		w = wmin
	}
	
	:For plotting IP3 processes:
	Qip3 = Q(IP3)
	Minf = minf(IP3,cam)
	Hinf = hinf(IP3,cam)
	TauM = taum(cam)
	TauH = tauh(IP3,cam)
}

NET_RECEIVE(weight (uS)) {

	nwt = w

	D = 1-(1-D*(1-U))*exp((t-tlast)/-tauD)
	tlast = t

	:AMPA
	Aa = Aa + w*factora*D      Ba = Ba + w*factora*D
	
	:NMDA
	An = An + w*factorn*D      Bn = Bn + w*factorn*D
	
	:mGluR5
	Am5 = Am5 + factorm5*D   	Bm5 = Bm5 + factorm5*D
	
	:mGluR2
	Am2 = Am2 + w*factorm2   Bm2 = Bm2 + w*factorm2
	
}

PROCEDURE calcfactors() { LOCAL tpa, tpn, tpm2, tpm5
	if (tau1a/tau2a > .9999) {
		tau1a = .9999*tau2a
	}
	if (tau1n/tau2n > .9999) {
		tau1n = .9999*tau2n
	}
	
	if (tau1m2/tau2m2 > .9999) {
		tau1m2 = .9999*tau2m2
	}
	
	if (tau1m5/tau2m5 > .9999) {
		tau1m5 = .9999*tau2m5
	}

	tpa = (tau1a*tau2a)/(tau2a - tau1a) * log(tau2a/tau1a)
	factora = -exp(-tpa/tau1a) + exp(-tpa/tau2a)
	factora = 1/factora
	
	tpn = (tau1n*tau2n)/(tau2n - tau1n) * log(tau2n/tau1n)
	factorn = -exp(-tpn/tau1n) + exp(-tpn/tau2n)
	factorn = 1/factorn
	
	tpm2 = (tau1m2*tau2m2)/(tau2m2 - tau1m2) * log(tau2m2/tau1m2)
	factorm2 = -exp(-tpm2/tau1m2) + exp(-tpm2/tau2m2)
	factorm2 = 1/factorm2	
	
	tpm5 = (tau1m5*tau2m5)/(tau2m5 - tau1m5) * log(tau2m5/tau1m5)
	factorm5 = -exp(-tpm5/tau1m5) + exp(-tpm5/tau2m5)
	factorm5 = 1/factorm5

}

FUNCTION s(v) {
	s = 1/(1+0.33*exp(-0.06*v))
}

FUNCTION eta(ca) {
	eta = .001/((.1/(.1/10000+ca*ca*ca))+1)
}

FUNCTION omegaN(ca) {
                if(ca<dth) {
                        omegaN = 0
                } else if(ca<pth) {
                        omegaN = -sqrt((pth-dth)/2*(pth-dth)/2-(ca-(dth+pth)/2)*(ca-(dth+pth)/2))
                } else {
			omegaN = 1/(1+50*exp(-50*(ca-pth)))
                }
	
}

FUNCTION omegaM(ca) {
                if(ca<2) {
                        omegaM = 0
                } else {
			omegaM = -1/(1+50*exp(-50*(ca-2))) :original value 2 then 1.5
                }
	
}

FUNCTION minf(IP3,ca) {
	minf = (IP3/(IP3+dip3))*(ca/(ca+dact))
}

FUNCTION hinf(IP3,ca) {
	hinf = Q(IP3)/(Q(IP3)+ca)
}

FUNCTION taum(ca) {
	taum = 1/(bip3 + aip3*ca)
}

FUNCTION tauh(IP3,ca) {
	tauh = 1/(ainh*(Q(IP3)+ca))
}

FUNCTION Q(IP3) {
	Q = dinh*(IP3+dip3)/(IP3+ddis)
}
