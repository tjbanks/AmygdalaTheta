
: Current clamp

UNITS {
	(nA) = (nanoamp)
	
}

NEURON {
	POINT_PROCESS thetastim2
	RANGE  amp, i, ton, start, tdone
	ELECTRODE_CURRENT i
}

PARAMETER {
	amp = 0
	ton = 1e30
	start = 1e30

	freq = 8
	tdone = 1e30
}

ASSIGNED { 
	i  (nA)

}

INITIAL { 
i = 0 
}

BREAKPOINT {
	at_time(ton)
	at_time(ton+0.3)
	at_time(start+1000)
	at_time(tdone)	

	if(t>(ton+0.3)) {
		i = 0
		ton = ton + 1000/freq
	} else if(t>ton) {
		i = amp
	}



	if(t>(start+1000)){
		i=0
		start = start+1500
		ton = start
	}



	if(t>=tdone) {
	  i = 0
	  start = 1e30
	  ton = 1e30
	}

}
