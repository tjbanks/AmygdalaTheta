NEURON {
	POINT_PROCESS progress
	RANGE tstop
}

PARAMETER {
	tstop = 1000
}

INITIAL { 
net_send(tstop/100,1)
}
	
	
NET_RECEIVE (w) {
if (flag==1) {
printf("%g percent complete\n",t/tstop*100)
net_send(tstop/100,1)
}
}
