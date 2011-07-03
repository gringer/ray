/*
 	Ray
    Copyright (C)  2010  Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (COPYING).  
	see <http://www.gnu.org/licenses/>

*/

#ifndef _ColorSpaceCodec
#define _ColorSpaceCodec

#include<string>
#include<stdint.h>

using namespace std;

class ColorSpaceCodec{
  static const char csColors[5];
  static const char bsBases[5];
 public:
	ColorSpaceCodec();
	static uint8_t csChrToInt(char tChr);
	static uint8_t bsChrToInt(char tChr);
	static char csChrToDE(char tChr);
	static char bsChrToBS(char tChr);
	static char bsIntToBS(uint8_t tInt);
	static char csIntToCS(uint8_t tInt, bool doubleEncoding);
	static int complement(int codeB);
	static bool isColorSpace(string sequence);
	static string transformCStoDE(string csInput);
	static int mapBStoCS(int mapX, int mapY);
	static int mapCStoBS(int mapB, int mapC);
	static int revMapCStoBS(int mapC, int mapB);
	static string decodeCStoBS(string csInput);
	static string decodeCStoBS(string csInput, bool reverseComplement);
	static string encodeBStoCS(string bsInput);
	static bool check();
};

#endif
