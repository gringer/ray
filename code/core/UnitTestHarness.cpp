/*
 * UnitTestHarness.cpp
 *
 *  Created on: 29/06/2011
 *      Author: David Eccles (gringer) <david.eccles@mpi-muenster.mpg.de>
 */

#include <string>
#include <iostream>
#include <core/UnitTestHarness.h>

UnitTestHarness::UnitTestHarness(string tDesc){
	m_testsDone = 0;
	m_classDescription = tDesc;
	m_successSoFar = true;
	cout << "Checking " << m_classDescription << ":" << endl;
	m_preProcessed = false;
}

bool UnitTestHarness::preProcessTest(string tDescription){
	if(m_preProcessed){
		cerr << "error: a test has already been pre-processed" << endl;
		return false;
	}
	cout << ++m_testsDone << ": checking " << tDescription << "... ";
	m_preProcessed = true;
	return true;
}

bool UnitTestHarness::compareTest(bool result, bool expectedResult){
	if(!m_preProcessed){
		cerr << "error: test has not been pre-processed" << endl;
		return false;
	}
	bool test = (result == expectedResult);
	m_successSoFar = m_successSoFar && test;
	if(!test) {
		cout << "failed\n"
				<< "  result:   " << result << endl
				<< "  expected: " << expectedResult << endl;
	} else {
		cout << "success!\n";
	}
	m_preProcessed = false;
	return test;
}

bool UnitTestHarness::compareOutput(string input, string expectedOutput){
	if(!m_preProcessed){
		cerr << "error: test has not been pre-processed" << endl;
		return false;
	}
	bool test = (input.compare(expectedOutput) == 0);
	m_successSoFar = m_successSoFar && test;
	if(!test) {
		cout << "failed\n"
				<< "  input:    " << input << endl
				<< "  expected: " << expectedOutput << endl;
	} else {
		cout << "success!\n";
	}
	m_preProcessed = false;
	return test;
}

bool UnitTestHarness::compareOutput(char* input, string expectedOutput){
	return(compareOutput(string(input), expectedOutput));
}

bool UnitTestHarness::getSuccess(){
	return m_successSoFar;
}
