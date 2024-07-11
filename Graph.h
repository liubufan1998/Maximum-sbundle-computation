#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph
{
private:
	std::string dir; // input graph directory
	ui n;			 // number of nodes of the graph
	ept m;			 // number of edges of the graph
	ui S;			 // the value of s in s-bundle

	ept *pstart; // offset of neighbors of nodes
	ept *pend;	 // used in search
	ept *pend_buf;
	ui *edges; // adjacent ids of edges
	ui *edgelist_pointer;

	std::vector<ui> sbundle;

public:
	Graph(const char *_dir, const int _K);
	~Graph();

	void read_binary_file();
	std::string get_file_name_without_suffix(const std::string &file_path);

	void output_one_sbundle();
	void verify_sbundle();

	void MSBP_SymBD_H();

private:
	void reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges);
	void extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges);
	void extract_subgraph_wo_prune(ui u, const ui *p_rid, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges);

	void write_subgraph(ui n, const std::vector<std::pair<int, int>> &edge_list);
	void extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);

	ui degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output);
	void degen2(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output);
	void core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output);
	void orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid);
	void oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj);
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf);
};
#endif
