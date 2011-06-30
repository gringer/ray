/*
 * UnitTestHarness.h
 *
 *  Created on: 29/06/2011
 *      Author: bioinf
 */

#ifndef UNITTESTHARNESS_H_
#define UNITTESTHARNESS_H_

#include <string>

using namespace std;

class UnitTestHarness{
	int m_testsDone;
	string m_classDescription;
	bool m_preProcessed;
	bool m_successSoFar;
public:
	UnitTestHarness(string tDescription);
	bool preProcessTest(string tDescription);
	bool runTest(string tDescription, bool test, string input, string expectedOutput);
	bool compareTest(bool result, bool expectedResult);
	bool compareOutput(string input, string expectedOutput);
	bool compareOutput(char* input, string expectedOutput);
	bool getSuccess();
};


#endif /* UNITTESTHARNESS_H_ */
