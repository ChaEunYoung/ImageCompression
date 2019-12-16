#pragma once


class rlc {
public:
	// member data
	float size, amplitude, rl,symbol;
	char code;
	

	// member function
	float getsize(float a);
	float getsymbol();
	
	//constructors
	rlc() { amplitude = 0; size = getsize(amplitude); } //symbol = getsymbol(rl, size);
	
rlc(float x, float y) { rl = x; amplitude = y; size = getsize(y); } //symbol = getsymbol(y,size); }
rlc(float y) { amplitude = y; size = getsize(amplitude); } //symbol = getsymbol(rl, size);
	


	rlc& operator=(const rlc&);

};