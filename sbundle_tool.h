#ifndef _SBUNDLE_TOOL_H
#define _SBUNDLE_TOOL_H

#include <stdio.h>
#include <utility>
#include <fstream>
#include <bits/stdc++.h>


extern unsigned int *head;
extern unsigned int *que, *nV, n;
const unsigned int INF = 0x3f3f3f3f;
using ui = unsigned int; // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges

class MaximumFlow
{
private:
    struct Node
    {
        ui from, to, next;
        ui cap;
    };
    Node *edge;
    ui *cap_backup;
    ui tol;
    ui *head;
    ui *dep;
    ui *gap, *que;
    ui *cur;
    ui *S;
    ui n;
    void BFS(ui start, ui end);
    ui *nV, *oID;
    bool flag;

public:
    MaximumFlow();
    ~MaximumFlow();
    void reserve(ui n, ui m);
    void init(ui _n);
    void addedge(ui u, ui v, ui w);
    ui SAP(ui start, ui end, ui nn);
    void BFSMinCut(ui start, std::vector<bool> &visited);
    void get_local_vertexConnectivity(ui R_end, ui *SR, char *matrix, long long row, ui *local_vc, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected_new, int s); 
    void prepare(ui n, ui m);                                                                                                                                    
    bool verify_SBundle_by_MaxFlowAlg0(ui* ids, ui R_end, char *matrix, long long row, int s, bool is_global);
    bool verify_SBundle_by_MaxFlowAlg(std::vector<ui> ids, ui R_end, char *matrix, long long row, int s, bool is_global);
    bool verify_Sbundle_by_maxFlowAlg2(ui* ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s);
    bool verify_Sbundle_by_maxFlowAlg2(ui* ids, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, int s, ui u, ui lb);
    bool verify_Sbundle_by_maxFlowAlg3(std::vector<ui> ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s);
    bool is_connected_graph(ui* ids, ui R_end, char *matrix, long long row);
    bool is_connected_graph(std::vector<ui> ids, ui R_end, char *matrix, long long row);
    void delete_vertices_based_on_lower_bound(ui *SR, ui S_end, ui R_end, char *matrix, long long row, int s, ui lower_bound, std::vector<ui> &results);
    void delete_vertices_based_on_lower_bound(ui *SR, ui S_end, ui R_end, ept *pstart, ept *pend, ui *edges, ui s, ui lower_bound, std::vector<ui> &results);
    void get_branching_case(ui *SR, ui &u1, ui &u2, std::vector<char> &s_connected, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, ui &min_vc, int s, ui n);
    void fill_all_s_connected_in_S(ui *SR, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected_new, int s, ui n);
};
#endif