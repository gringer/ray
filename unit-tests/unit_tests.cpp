/*
 * test_cscodec.cpp
 *
 *  Created on: 22/06/2011
 *      Author: David Eccles (gringer) <david.eccles@mpi-muenster.mpg.de>
 *
 *  Checks to make sure colour space encoding / decoding works correctly
 */

#include<format/ColorSpaceCodec.h>
#include<structures/Read.h>

int main(){
	int result = 0;
	cout << "Checking ColorSpaceCodec:\n";
	result += ColorSpaceCodec::check()?0:1;
	result += Read::check()?0:1;
	return result;
}
