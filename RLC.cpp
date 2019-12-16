#include "rlc.h"
#include <iostream>
using namespace std;

float rlc::getsize(float amplitude) {
	if (amplitude == 0) { size = 0; code = NULL; }
	else if (amplitude == -1 || amplitude == 1) size = 1;
	else if ((amplitude >= -3 && amplitude <= -2) || (amplitude >= 2 && amplitude <= 3)) size = 2;
	else if ((amplitude >= -7 && amplitude <= -4) || (amplitude >= 4 && amplitude <= 7)) size = 3;
	else if ((amplitude >= -15 && amplitude <= -8) || (amplitude >= 8 && amplitude <= 15)) size = 4;
	else if ((amplitude >= -31 && amplitude <= -16) || (amplitude >= 16 && amplitude <= 31)) size = 5;
	else if ((amplitude >= -63 && amplitude <= -32) || (amplitude >= 32 && amplitude <= 63)) size = 6;
	else if ((amplitude >= -127 && amplitude <= -64) || (amplitude >= 64 && amplitude <= 127)) size = 7;
	else if ((amplitude >= -255 && amplitude <= -128) || (amplitude >= 128 && amplitude <= 255)) size = 8;
	else if ((amplitude >= -511 && amplitude <= -256) || (amplitude >= 256 && amplitude <= 511)) size = 9;
	else if ((amplitude >= -1023 && amplitude <= -512) || (amplitude >= 512 && amplitude <= 1023)) size = 10;
	else { cout << "dpcm error!" << endl; }
	return size;
}

float rlc::getsymbol() {
	if (size > 10) {
		return rl * 100 + size;
	}
	else{
		return rl * 10 + size;
	}
	
}

rlc& rlc::operator=(const rlc& x)
{
	rl = x.rl;
	size = x.size;
	amplitude = x.amplitude;
	return *this;
}
