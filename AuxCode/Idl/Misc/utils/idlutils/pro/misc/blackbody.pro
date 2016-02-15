function blackbody, lambda, T, retfnu=retfnu

	; Lambda is in angstroms

	c = 2.99792e18 ; Ang/s	
	hc = 1.988e-6  ; erg * A

	wave = double(lambda)

	;
	;	Fnu is prop to ergs/s cm^-2 Hz^-1 
	;
	fnu = hc/(wave^3) * 1.0/(exp(1.43868e8/(wave*T)) - 1.0)

	if (keyword_set(retfnu)) then return, fnu

	flam = fnu*c/wave^2 

	;
	;	flam is prop to ergs/s cm^-2 Ang^-1 
	;

	return, flam
end 


