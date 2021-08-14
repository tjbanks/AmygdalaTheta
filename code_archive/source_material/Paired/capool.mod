: Calcium Pool dynamics including ca-current dependent influx, buffering


NEURON {
  SUFFIX pool
  USEION ca READ ica WRITE cai
  USEION cal READ ical VALENCE 2
  RANGE f, tau, ca0
}

INITIAL {
	cai = 0.05e-3  : (mM)
}

PARAMETER {
	f = 0   : mM/pA.ms
	tau = 100
	ca0 = 0.05e-3

}

STATE { cai

} 

ASSIGNED {

 ica
 RTbyZF
 ical
 
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  
}


DERIVATIVE states { 

	cai'= -f*(ica+ical) + (ca0 - cai)/tau

}
