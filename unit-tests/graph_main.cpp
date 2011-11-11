#include <stdlib.h>
#include <routing/ConnectionGraph.h>

int main(int argc,char**argv){

	int n=atoi(argv[1]);

	ConnectionGraph graph;

	graph.buildGraph(n,argv[2],true);

	if(atoi(argv[3]))
		graph.writeFiles("./");

	return 0;
}
