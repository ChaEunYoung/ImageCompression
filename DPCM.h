#pragma once

class dpcm {
public:
	// member data
	float size, amplitude;
	char code;
	// member function
	double getsize(double a);
	
	//constructors
	dpcm() {amplitude = 0; size = getsize(amplitude);};
	dpcm(double x, double y) { size = x; amplitude = y; };
	dpcm(double y) { amplitude = y; size = getsize(amplitude); };
	

	dpcm& operator=(const dpcm&);

};