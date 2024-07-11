#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage() {
	printf("Example usage: ./SymBK_H -g path_to_graph -s 3\n");
}

int main(int argc, char *argv[])
{
#ifndef NDEBUG
	printf("**** SymBK_H (Debug) build at %s %s ***\n", __TIME__, __DATE__);
	printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
	printf("**** SymBK_H (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	bool output = false;
	bool mode = false;
	string alg;

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto s_option = op.add<Value<int>>("s", "s", "\'the value of s for s-bundle\'");

	op.parse(argc, argv);

	if(help_option->is_set()||argc <= 1) {
		cout << op << endl;
		if(argc <= 1) {
			print_usage();
			return 0;
		}
	}
	if(!graph_option->is_set()) {
		printf("!!! The argument -g is required! Exit !!!\n");
		print_usage();
		return 0;
	}
	if(!s_option->is_set()) {
		printf("!!! The argument -s is required! Exit !!!\n");
		print_usage();
		return 0;
	}

	Graph *graph = new Graph(graph_option->value().c_str(), s_option->value());
    graph->read_binary_file();

#ifndef NDEBUG
	printf("\t*** Finished reading graph\n");
#endif

	graph->MSBP_SymBD_H(); // For obtaining the exact maximum sbundle

	delete graph;

	printf("\n");
	return 0;	
}
