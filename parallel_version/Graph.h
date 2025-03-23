#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph
{
private:
	string dir; // input graph directory
	ui n;		// number of nodes of the graph
	ept m;		// number of edges of the graph
	ui S;		// the value of s in s-bundle

	ept *pstart; // offset of neighbors of nodes
	ept *pend;	 // used in search
	ept *pend_buf;
	ui *edges; // adjacent ids of edges
	ui *edgelist_pointer;

	vector<ui> sbundle;

public:
	Graph(const char *_dir, const int _K);
	~Graph();

	string get_file_name_without_suffix(const string &file_path);
	void readBinFile();
	void output_one_sbundle();
	void verify_sbundle();
	void SymBD_Parallel(const int number_parallel_threads);

private:
	void reorganize_adjacency_lists(ui n, vector<ui> &peel_sequence, vector<ui> &rid, ui *pstart, ui *pend, ui *edges);
	void extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const vector<ui> &p_rid, vector<ui> &degree, vector<ui> &ids, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges);
	void extract_subgraph_wo_prune(ui u, const vector<ui> &p_rid, vector<ui> &ids, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges);

	void write_subgraph(ui n, const vector<pair<int, int>> &edge_list);
	void extract_subgraph(ui u, ui *ids, ui &ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, ui *Q, vector<ui> &degree, vector<char> &exists, ept *pend, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_full(const ui *ids, ui ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);

	ui degen(ui n, vector<ui> &peel_sequence, vector<ui> &core, ept *pstart, ui *edges, vector<ui> &degree, char *vis, ListLinearHeap *heap, bool output);
	void ego_degen(ui n, ui m, vector<ui> &peel_sequence, ept *pstart, ui *edges, vector<ui> &degree, vector<ui> &rid, char *vis, ListLinearHeap *heap, bool output);
	void core_shrink_graph(ui &n, ept &m, vector<ui> &peel_sequence, vector<ui> &core, vector<ui> &out_mapping, ui *in_mapping, vector<ui> &rid, ept *&pstart, ui *&edges, bool output);
	void orient_graph(ui n, ui m, vector<ui> &peel_sequence, ept *pstart, ept *pend, ui *edges, vector<ui> &rid);
	void oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj);
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf);
	ept peeling(ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, vector<ui> &degree, ept *pstart, ept *pend, ui *edges, vector<char> &exists);
	char find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges);
};
#endif
