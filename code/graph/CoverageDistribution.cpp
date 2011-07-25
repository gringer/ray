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

#include <graph/CoverageDistribution.h>
#include <iostream>
#include <fstream>
#include <map>
#ifdef ASSERT
#include <assert.h>
#endif
using namespace std;

CoverageDistribution::CoverageDistribution(map<int,uint64_t>*distributionOfCoverage,string*file){
	if(file!=NULL){
		ofstream f;
		f.open(file->c_str());
		for(map<int,uint64_t>::iterator i=distributionOfCoverage->begin();i!=distributionOfCoverage->end();i++){
			f<<""<<i->first<<" "<<i->second<<endl;
		}
		f.close();
	}
	
	int windowSize=10;
	int minimumX=1;
	uint64_t minimumY=2*4096;
	uint64_t minimumY2=55000;
	int maximumX=65535-1;
	int safeThreshold=256;

	vector<int> x;
	vector<uint64_t> y;
	uint64_t lastY = 0;
	bool slidePassed = false;

	for(map<int,uint64_t>::iterator i=distributionOfCoverage->begin();i!=distributionOfCoverage->end();i++){
		int tx = i->first;
		uint64_t ty = i->second;
		if (!slidePassed && (lastY != 0) && (lastY < ty)){
			slidePassed = true;
		}
		lastY = ty;
		if((tx >= minimumX) && (tx <= maximumX) && slidePassed){
			if((tx < safeThreshold) || (ty >= minimumY2)){
				x.push_back(tx);
				y.push_back(ty);
			}
		}
	}

	/** get the votes to find the peak */
	uint8_t votes[x.size()];
	for(int i = 0; i < (int)x.size(); i++){
		votes[i] = 0; // initialise to zero, if needed
	}
	for(int i=0;i<(int)x.size();i++){
		int largestPosition=i;
		for(int position=i;(position < (int)x.size()) && ((position - i) < windowSize); position++){
			if(y[position] > y[largestPosition])
				largestPosition=position;
		}
		if(y[largestPosition]>minimumY){
			cout << "y = " << y[largestPosition] << ", minimumY = " << minimumY << endl;
			votes[largestPosition]++;
		}
	}
	/** check votes */
	/** Select the last peak with any votes, or fall back to the highest value (i.e. max(y))
	 * if there are no votes */
	int largestPosition=0;
	for(int i=0; i < (int)x.size(); i++){
		//bool changed = false;
		if((votes[i] > votes[largestPosition]) || (y[i] > y[largestPosition])){
			largestPosition = i;
			//changed = true;
		}
		//cout << "x = " << x[i] << ", y = " << y[i] <<
		//		", vote = " << (int)votes[i] << (changed?"*":"") << endl;
	}
	
	m_minimumCoverage= distributionOfCoverage->begin()->first;
	m_peakCoverage=x[largestPosition];

	m_repeatCoverage=2*m_peakCoverage;
	int diff=m_peakCoverage-m_minimumCoverage;
	int candidate=m_peakCoverage+diff;

	if(candidate<m_repeatCoverage)
		m_repeatCoverage=candidate;
}

int CoverageDistribution::getPeakCoverage(){
	return m_peakCoverage;
}

int CoverageDistribution::getMinimumCoverage(){
	return m_minimumCoverage;
}

int CoverageDistribution::getRepeatCoverage(){
	return m_repeatCoverage;
}
