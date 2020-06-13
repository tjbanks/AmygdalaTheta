
NEURON {
	POINT_PROCESS Newglu
	USEION ca READ eca 
	USEION cal READ ical VALENCE 2
	
	RANGE Cdur_nmda, AlphaTmax_nmda, Beta_nmda, Erev_nmda, gbar_nmda, W_nmda, on_nmda, g_nmda
	RANGE Cdur_ampa, AlphaTmax_ampa, Beta_ampa, Erev_ampa, gbar_ampa, wn, on_ampa, g_ampa
	RANGE eca, P0, fCa, tauCa
	RANGE Cainf, pooldiam, z
	RANGE tau1m2, tau2m2
	RANGE tau1m5, tau2m5
	RANGE factora, factorn, factorm2, factorm5, initW, fmax, fmin, wmax, wmin, w
	RANGE Vmax, Kl
	RANGE L1, L2, L3, L4, L5
	RANGE Vm5, Vm2
	RANGE fm, fn
	RANGE source, dest
	RANGE tauD, U
	RANGE Qip3, Minf, Hinf, TauM, TauH, IP3
	RANGE Jnmda, Jpmcan,	Jncxn,	Jleak, dth, pth
	RANGE on_mglu
	RANGE F, f, tauF, D1, d1, tauD1, D2, d2, tauD2
	RANGE facfactor
	
	NONSPECIFIC_CURRENT i,iampa,inmda
}

UNITS {
	(mV) = (millivolt)
        (nA) = (nanoamp)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	pi = 3.141592 (1)
}

PARAMETER {

	Cdur_nmda = 16.7650 (ms)
	AlphaTmax_nmda = .2659 (/ms)
	Beta_nmda = 0.008 (/ms)
	Erev_nmda = 0 (mV)
	gbar_nmda = .5e-3 (uS)

	Cdur_ampa = 1.4210 (ms)
	AlphaTmax_ampa = 3.8142 (/ms)
	Beta_ampa = 0.1429 (/ms)
	Erev_ampa = 0 (mV)
	gbar_ampa = 1.18e-3 (uS)   :1e-3 (uS)
	


	Cainf = 50e-6 (mM)
	pooldiam =  1.8172 (micrometer)
	z = 2

	tauCa = 50 (ms)
	P0 = 0.13  :.015
	fCa = .024
	
	Vm5 = 0.25 
	Vm2 = 1    :Vm5 = 0.25

	tau1m5 = 100 tau2m5 = 400   :250, 2500
	tau1m2 = 100 tau2m2 = 200
	
	
    a = 0.00412      fcytratio = 0.00280899   :a = 0.05
	fcytbeta = 0.0095   ferratio = 0.015873   :fcytbeta = 1  :fcytbeta original .0035
	:fcytbeta= (.0088 original) (.0095 works!)   :ferratio = 0.015873					
	Vn = 0.04 (mM/ms)   Kn = 1.0  (mM)      :Vn = 0.04
	Vp = 0.04 (mM/ms)   Kp = 0.1  (mM)
	Vs = 13 (mM/ms)   Ks = 0.2   (mM)      :Vs = 1.3
	Vmax = 0.1 (mM/ms)  Kl = 0.1  (mM)
	Vsoc = 0.03 (mM/ms) Ksoc = 100 (mM)     :Vsoc = 0.01  Ksoc = 1
		
	fn = 0.22
	fm = 6e5                      :fm = 20  //fm original value 6e5
						: fn 0.15
	Vr = 15   Vl = 0.0021  :(pL/ms), (pL/ms)
	
	bip3 = 0.1e-3  :4.1e-3     
	aip3 = 0.42          : (1/ms), (uM-ms)^-1
	ainh = 0.015  :0.005       
	dinh = 0.525          : (uM-ms)^-1, (uM)    :ainh = 2e-4  dinh = 1.05
	dip3 = 0.13       ddis = 0.47          : (uM), (uM)            :ddis = 0.94
	dact = 0.6                          : (uM)                   :dact = 8.2e-2
	
	dth = 0.5  (uM):1.6    
	pth = 0.7 (uM) 	:1.75   :(uM), (uM) 0.3 0.6  0.8 may be good value for duh
	L1 = 0.5       L2 = 0.097    L3 = 0.05  L4 = 0.004  L5 = 0.006
			:L2 0.025
	tauD = 50   U = 0.4
	
	initW = 1
	fmax = 2  :1.75 :3
	fmin = .2
	
	facfactor = 1

	: the (1) is needed for the range limits to be effective
        f = 0 (1) < 0, 1e9 >    : facilitation
        tauF = 1 (ms) < 1e-9, 1e9 >
        d1 = 0.5 (1) < 0, 1 >     : fast depression
        tauD1 = 25 (ms) < 1e-9, 1e9 >
        d2 = 0.5 (1) < 0, 1 >     : slow depression
        tauD2 = 25 (ms) < 1e-9, 1e9 >		

}

ASSIGNED {
	Qip3 Hinf Minf TauM TauH

	Jnmda	Jpmcan	Jncxn	Jleak
	
	factorm2 factorm5
	i
	ical eca
	nwt 
	v
	w
	IP3
	
	Jip3r 	Jvgcc 	Jpmcam	Jncxm 	Jserca	Jsoc CanRa
	Jin Jout

	source dest
	D
	tlast
	
	wmin
	wmax
	
	inmda (nA)
	g_nmda (uS)
	on_nmda
	W_nmda
	iampa (nA)
	g_ampa (uS)
	on_ampa
	Afactor	(mM/ms/nA)

	t0 (ms)
	on_mglu
	tsyn
	F
	D1
	D2
}

STATE {
	r_nmda r_ampa
	can cam cae wm wn
	m h
	Am5 Bm5
	Am2 Bm2
}

INITIAL {
	on_nmda = 0
	r_nmda = 0
	W_nmda = initW

	on_ampa = 0
	r_ampa = 0
	wn = initW/2
	wmax = fmax*initW
	wmin = fmin*initW
	calcfactors()

	on_mglu = 0
	
	can = Cainf
	Afactor	= 1/(z*FARADAY*4/3*pi*(pooldiam/2)^3)*(1e6)

	Jnmda = 0
	cam = 0.05  cae = 120   : (uM), (uM), (uM)   :cae = 4
	m = minf(IP3,cam)  h = hinf(IP3,cam)
	D = 1
	tlast = -1e30
	t0 = -1
	wm = initW/2
	w = wn+wm
	
	F = 1
	D1 = 1
	D2 = 1

	}

BREAKPOINT {

	SOLVE state METHOD cnexp
	

}

DERIVATIVE state { 

	Am5' = -Am5/tau1m5   Bm5' = -Bm5/tau2m5
	Am2' = -Am2/tau1m2   Bm2' = -Bm2/tau2m2
	
	
	
:NMDA Ca2+ pool

if (t0>0) {
			
			if (t-t0 < Cdur_nmda) {
				on_nmda = 1
			} else {
				on_nmda = 0
			}
			if (t-t0 < Cdur_ampa) {
				on_ampa = 1
			} else {
				on_ampa = 0
			}
			
			if (t-t0 < 100) {
				on_mglu = 1
			} else {
				on_mglu = 0
			}
			
}
	r_nmda' = AlphaTmax_nmda*on_nmda*(1-r_nmda)-Beta_nmda*r_nmda
	r_ampa' = AlphaTmax_ampa*on_ampa*(1-r_ampa)-Beta_ampa*r_ampa
	

	:gbar_nmda = (0.85/(1+0.7*exp(-0.4*cam)))*1e-3
	gbar_nmda = (1.2/(1+1.4*exp(-0.4*cam)))*1e-3
	:gbar_nmda = (1.05/(1+1.4*exp(-0.4*cam)))*1e-3
	
	
	g_nmda = gbar_nmda*r_nmda*facfactor
	inmda = W_nmda*g_nmda*(v - Erev_nmda)*s(v)

	g_ampa = gbar_ampa*r_ampa*facfactor
	iampa = w*g_ampa*(v - Erev_ampa)
	
	i = iampa+inmda
	
	
	Jnmda = -fCa*Afactor*P0*g_nmda*(v - eca)*s(v)*1e3 :+ L4*cam
	
	Jpmcan = Vp*can^2/(Kp^2+can^2)
	Jncxn = Vn*can^4/(Kn^4+can^4)
	Jleak = Vmax*(0.05)^4/(0.05^4+Kl^4)
	
	can' = Jnmda + Jleak - Jpmcan - Jncxn 

:mGluR5 Ca2+ pool



    IP3 = Vm5*(Bm5 - Am5)
	
	m' = (minf(IP3,cam)-m)/taum(cam)
	h' = (hinf(IP3,cam)-h)/tauh(IP3,cam)
	
	Jip3r = (Vr*m^3*h^3 + Vl) *(cae-cam) 
	Jvgcc = -fm*ical*on_mglu 
	Jpmcam = Vp*cam^2/(Kp^2+cam^2) 
	Jncxm = Vn*cam^4/(Kn^4+cam^4)	
	Jserca = Vs*cam^2/(Ks^2+cam^2)
	Jsoc =  -Vsoc*Ksoc^4/(Ksoc^4+cae^4)*(v-eca)
	
	Jin = a*(Jvgcc + Jsoc)
	Jout = Jpmcam + Jncxm
	CanRa =  L5*can
	
	cam' = fcytbeta*(Jin - Jout) + fcytratio*(Jip3r - Jserca) + CanRa
	cae' = ferratio*(Jserca - Jip3r)
	  
	wn' = eta(can)*(L1*omegaN(can) - L3*(w-initW))  :- L4*(wm-initW/2)) 
	wm' = eta(cam/4)*(L2*omegaM(cam) - L3*(w-initW))  :- L5*(wn-initW/2)) 
	
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
	t0 = t
	nwt = w	

	D = 1-(1-D*(1-U))*exp((t-tlast)/-tauD)
	tlast = t

	
	:mGluR5
	Am5 = Am5 + factorm5*D   	Bm5 = Bm5 + factorm5*D
	
	:mGluR2
	Am2 = Am2 + w*factorm2   Bm2 = Bm2 + w*factorm2
	
	:F  = 1 + (F-1)* exp(-(t - tsyn)/tauF)
	D1 = 1 - (1-D1)*exp(-(t - tsyn)/tauD1)
	D2 = 1 - (1-D2)*exp(-(t - tsyn)/tauD2)

	tsyn = t
	
	facfactor = F * D1 * D2
	:F = F+f  :F * f
	
	if (F > 3) { 
	F=3	}
	:if (facfactor < 0.5) { 
	:facfactor=0.5
	:}	
	D1 = D1 * d1
	D2 = D2 * d2

	
}

PROCEDURE calcfactors() { LOCAL tpa, tpn, tpm2, tpm5
	
	if (tau1m2/tau2m2 > .9999) {
		tau1m2 = .9999*tau2m2
	}
	
	if (tau1m5/tau2m5 > .9999) {
		tau1m5 = .9999*tau2m5
	}

	
	tpm2 = (tau1m2*tau2m2)/(tau2m2 - tau1m2) * log(tau2m2/tau1m2)
	factorm2 = -exp(-tpm2/tau1m2) + exp(-tpm2/tau2m2)
	factorm2 = 1/factorm2	
	
	tpm5 = (tau1m5*tau2m5)/(tau2m5 - tau1m5) * log(tau2m5/tau1m5)
	factorm5 = -exp(-tpm5/tau1m5) + exp(-tpm5/tau2m5)
	factorm5 = 1/factorm5

}

FUNCTION s(v) {
	s = 1/(1+0.33*exp(-0.13*v))
	:s = 1/(1+0.33*exp(-0.2*v))
	:s = 1/(1+0.33*exp(-0.06*v))
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
                if(ca<1.5) {
                        omegaM = 0
                } else {
			omegaM = -1/(1+50*exp(-50*(ca-1.5))) :original value 2 then 1.5
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
