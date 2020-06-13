: passive leak current

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT il
	RANGE il, e, g
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	g = 3.333333e-5 (siemens/cm2) < 0, 1e9 >
	e = -75 (mV)
}

ASSIGNED {
	v (mV)
	i
	il (mA/cm2)
}

BREAKPOINT { 
	i = g*(v - e)
	il = i
}
