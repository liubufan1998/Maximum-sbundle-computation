#include "Graph.h"
#include "Utility.h"
#include "Timer.h"

using namespace std;

int main(int argc, char *argv[])
{

	char filename[150]; 
	assert(argc == 4);	
	strcpy(filename, argv[1]);
	int s = std::stoi(argv[2]);						  
	int number_parallel_threads = std::stoi(argv[3]); 
	Graph *graph = new Graph(filename, s);			
	graph->readBinFile();
	graph->SymBD_Parallel(number_parallel_threads);

	printf("everything is OK!!!");
	delete graph;
	return 0;
}
