#ifndef _SBUNDLE_TOOL_H
#define _SBUNDLE_TOOL_H

#include <stdio.h>
#include <utility>
#include <fstream>
#include <bits/stdc++.h>

const unsigned int INF = 0x3f3f3f3f;
using ui = unsigned int; 
using ept = unsigned int; 
using namespace std;

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
    ui m; 
    ui *nV, *oID;
    bool flag; 
    vector<bool> visited;
    vector<ui> real_minCut; 

public:
    MaximumFlow(); 
    ~MaximumFlow();
    void BFS(ui start, ui end);
    void reserve(ui n, ui m);
    void init(ui _n);
    void addedge(ui u, ui v, ui w);
    ui SAP(ui start, ui end, ui nn, vector<ui>& real_minCut);
    void BFSMinCut(ui start, vector<bool> &visited);
    void get_local_vertexConnectivity(ui R_end, vector<ui> &SR, char *matrix, long long row, vector<ui> &local_vc, map<pair<ui, ui>, vector<char>> &all_s_connected_new, int s); 
    void prepare(ui n, ui m);                                                                                                                                    
    bool verify_SBundle_by_MaxFlowAlg0(vector<ui> &ids, ui R_end, char *matrix, long long row, int s, bool is_global); 
    bool verify_Sbundle_by_maxFlowAlg2(vector<ui> &ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s); 
    bool verify_Sbundle_by_maxFlowAlg2(vector<ui> &ids, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, int s, ui u, ui lb); 
    bool verify_Sbundle_by_maxFlowAlg3(vector<ui> ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s); 
    bool is_connected_graph(vector<ui> &ids, ui R_end, char *matrix, long long row); 
    void delete_vertices_based_on_lower_bound(vector<ui> &SR, ui S_end, ui R_end, char *matrix, long long row, int s, ui lower_bound, vector<ui> &results); 
    void delete_vertices_based_on_lower_bound(vector<ui> &SR, ui S_end, ui R_end, ept *pstart, ept *pend, ui *edges, ui s, ui lower_bound, vector<ui> &results); 
    void get_branching_case(vector<ui> &SR, ui &u1, ui &u2, vector<char> &s_connected, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, ui &min_vc, int s, ui n); 
    void fill_all_s_connected_in_S(vector<ui> &SR, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, map<pair<ui, ui>, vector<char>> &all_s_connected_new, int s, ui n); 
    
};
#endif