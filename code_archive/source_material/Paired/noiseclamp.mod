: Current clamp

UNITS {
	(nA) = (nanoamp)
	
}

NEURON {
	POINT_PROCESS noiseclamp
	RANGE i, noise, invl
	ELECTRODE_CURRENT i
}

PARAMETER {
	noise = 0	(nA)
	invl = 1
}

ASSIGNED { 

	i  (nA)

}

INITIAL { 

i = 0 
noise = 0
  net_send(0, 1)
}

BREAKPOINT {

 i = i + 0

}

NET_RECEIVE (w) {
  : respond only to self-events with flag==1
  if (flag==1) {

	i = noise
  
    net_send(invl, 1)
  }
}
