#include <assembler/Loader.h>
#include <structures/Read.h>
#include <stdint.h>
#include <iostream>
using namespace std;

int main(int argc,char**argv){
	if(argc!=2){
		cout<<"Provide an SFF file."<<endl;
		return 0;
	}
	string file=argv[1];
	Loader loader;
	loader.constructor("",false);
	loader.load(file,false);
	for(uint64_t i=0;i<loader.size();i++){
		string read = loader.at(i)->getSeq(false);
		cout<<">"<<i<<endl<<read<<endl;
	}
	return 0;
}
