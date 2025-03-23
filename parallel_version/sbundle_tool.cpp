#include <stdio.h>
#include "sbundle_tool.h"

MaximumFlow::MaximumFlow()
{
    edge = nullptr;
    cap_backup = nullptr;
    head = nullptr;
    dep = nullptr;
    gap = nullptr;
    que = nullptr;
    cur = nullptr;
    S = nullptr;
    nV = nullptr;
    oID = nullptr;
}

MaximumFlow::~MaximumFlow()
{
    if (edge != nullptr)
    {
        delete[] edge;
        edge = nullptr;
    }
    if (cap_backup != nullptr)
    {
        delete[] cap_backup;
        cap_backup = nullptr;
    }
    if (head != nullptr)
    {
        delete[] head;
        head = nullptr;
    }
    if (dep != nullptr)
    {
        delete[] dep;
        dep = nullptr;
    }
    if (gap != nullptr)
    {
        delete[] gap;
        gap = nullptr;
    }
    if (que != nullptr)
    {
        delete[] que;
        que = nullptr;
    }
    if (cur != nullptr)
    {
        delete[] cur;
        cur = nullptr;
    }
    if (S != nullptr)
    {
        delete[] S;
        S = nullptr;
    }
    if (nV != nullptr)
    {
        delete[] nV;
        nV = nullptr;
    }
    if (oID != nullptr)
    {
        delete[] oID;
        oID = nullptr;
    }
}

void MaximumFlow::BFS(ui start, ui end)
{
    memset(dep, -1, n * sizeof(ui));
    memset(gap, 0, n * sizeof(ui));
    gap[0] = 1;
    ui front, rear;
    front = rear = 0;
    dep[end] = 0;
    que[rear++] = end;
    while (front != rear)
    {
        ui u = que[front++];
        for (ui i = head[u]; i != -1; i = edge[i].next)
        {
            ui v = edge[i].to;
            if (dep[v] != -1)
                continue;
            que[rear++] = v;
            dep[v] = dep[u] + 1;
            ++gap[dep[v]];
        }
    }
}

void MaximumFlow::reserve(ui n, ui m)
{
    n += 3;
    edge = new Node[m * 4 + n * 2];
    cap_backup = new ui[m * 4 + n * 2];
    head = new ui[n + n];
    dep = new ui[n + n];
    gap = new ui[n + n];
    que = new ui[n + n];
    cur = new ui[n + n];
    S = new ui[n + n];
}

void MaximumFlow::init(ui _n)
{
    tol = 0;
    n = _n + 2;
    memset(head, -1, n * sizeof(ui));
}

void MaximumFlow::addedge(ui u, ui v, ui w)
{
    edge[tol].from = u;
    edge[tol].to = v;
    edge[tol].cap = w;
    cap_backup[tol] = w;
    edge[tol].next = head[u];
    head[u] = tol++;

    edge[tol].from = v;
    edge[tol].to = u;
    edge[tol].cap = 0;
    cap_backup[tol] = 0;
    edge[tol].next = head[v];
    head[v] = tol++;
}

ui MaximumFlow::SAP(ui start, ui end, ui new_size, vector<ui> &real_minCut)
{
    real_minCut.clear();

    for (ui i = 0; i < tol; ++i)
        edge[i].cap = cap_backup[i];
    ui res = 0;
    BFS(start, end);

    ui top = 0;
    memcpy(cur, head, n * sizeof(ui));
    ui u = start;
    ui i;
    while (dep[start] < n)
    {
        if (u == end)
        {

            ui temp = INF + 1;
            ui inser = 0;
            for (i = 0; i < top; ++i)
                if (temp > edge[S[i]].cap)
                {
                    temp = edge[S[i]].cap;
                    inser = i;
                }

            for (i = 0; i < top; ++i)
            {
                edge[S[i]].cap -= temp;
                edge[S[i] ^ 1].cap += temp;
            }
            res += temp;
            top = inser;
            u = edge[S[top]].from;
        }
        if (u != end && gap[dep[u] - 1] == 0)
            break;
        for (i = cur[u]; i != -1; i = edge[i].next)
            if (edge[i].cap != 0 && dep[u] == dep[edge[i].to] + 1)
                break;

        if (i != -1)
        {
            cur[u] = i;
            S[top++] = i;
            u = edge[i].to;
        }
        else
        {
            ui min = n;
            for (i = head[u]; i != -1; i = edge[i].next)
            {
                if (edge[i].cap == 0)
                    continue;
                if (min > dep[edge[i].to])
                {
                    min = dep[edge[i].to];
                    cur[u] = i;
                }
            }
            --gap[dep[u]];
            dep[u] = min + 1;
            ++gap[dep[u]];

            if (u != start)
                u = edge[S[--top]].from;
        }
    }

    vector<bool> visited(n, false);
    BFSMinCut(start, visited);
    for (ui i = 0; i < tol; i += 2)
    {
        ui u = edge[i].from;
        ui v = edge[i].to;

        if (visited[u] && !visited[v])
        {
            assert(u < v);
            real_minCut.push_back(u);
        }
    }
    assert(real_minCut.size() == res);
    return res;
}

void MaximumFlow::BFSMinCut(ui start, vector<bool> &visited)
{
    queue<ui> q;
    q.push(start);
    visited[start] = true;
    while (!q.empty())
    {
        ui u = q.front();
        q.pop();
        for (ui i = head[u]; i != -1; i = edge[i].next)
        {
            ui v = edge[i].to;

            if (!visited[v] && edge[i].cap > 0)
            {
                visited[v] = true;
                q.push(v);
            }
        }
    }
}

void MaximumFlow::prepare(ui _n, ui _m)
{
    n = _n;
    m = _m;
}

void MaximumFlow::get_local_vertexConnectivity(ui R_end, vector<ui> &SR, char *matrix, long long row, vector<ui> &local_vc, map<pair<ui, ui>, vector<char>> &all_s_connected_new, int s)
{

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    ui new_size = 0;
    for (ui i = 0; i < R_end; i++)
        nV[SR[i]] = new_size, oID[new_size++] = SR[i];

    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
        mf.addedge(nV[SR[i]], nV[SR[i]] + new_size, 1);
    for (ui i = 0; i < R_end; i++)
    {
        char *t_matrix = matrix + SR[i] * row;
        for (ui j = 0; j < R_end; j++)
        {
            if (!t_matrix[SR[j]])
                continue;
            mf.addedge(nV[SR[i]] + new_size, nV[SR[j]], INF);
        }
    }

    for (ui i = 0; i < new_size; ++i)
    {
        char *t_matrix = matrix + oID[i] * row;
        for (ui j = i + 1; j < new_size; ++j)
        {

            if (t_matrix[oID[j]])
                local_vc[oID[i] * row + oID[j]] = local_vc[oID[j] * row + oID[i]] = row;
            else
            {
                ui t_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                assert(t_vc == real_minCut.size());
                local_vc[oID[i] * row + oID[j]] = local_vc[oID[j] * row + oID[i]] = t_vc;

                ui temp_u = oID[i], temp_v = oID[j];
                if (temp_u > temp_v)
                {
                    ui ttemp = temp_u;
                    temp_u = temp_v;
                    temp_v = ttemp;
                }

                assert(temp_u < temp_v);
                pair<ui, ui> key = make_pair(temp_u, temp_v);
                vector<char> value(row, 0);
                for (ui k = 0; k < real_minCut.size(); k++)
                    value[oID[real_minCut[k]]] = 1;
                all_s_connected_new[key] = value;
            }
        }
    }
}

bool MaximumFlow::verify_SBundle_by_MaxFlowAlg0(vector<ui> &SR, ui R_end, char *matrix, long long row, int s, bool is_global)
{
    if (R_end <= s)
        return true;

    assert(R_end > s);

    if (!is_connected_graph(SR, R_end, matrix, row))
        return false;

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }

    mf.init(new_size + new_size);

    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }

    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        char *t_matrix = matrix + u * row;
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];
            if (!t_matrix[v])
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    if (is_global)
    {
        for (ui i = 0; i < new_size; ++i)
        {
            char *t_matrix = matrix + oID[i] * row;
            for (ui j = i + 1; j < new_size; ++j)
                if (!t_matrix[oID[j]] && mf.SAP(new_size + i, j, new_size, real_minCut) < R_end - s)
                    return false;
        }
    }
    else
    {
        ui u = oID[new_size - 1];
        assert(u == SR[R_end - 1]);
        char *t_matrix = matrix + u * row;
        for (ui i = 0; i < new_size - 1; ++i)
        {
            assert(oID[i] != u);
            if (!t_matrix[oID[i]] && mf.SAP(new_size + new_size - 1, i, new_size, real_minCut) < R_end - s)
            {
                return false;
            }
        }
    }
    return true;
}

bool MaximumFlow::verify_Sbundle_by_maxFlowAlg2(vector<ui> &SR, ui R_end, ui *pstart, ui *pend, ui *edges, int s)
{
    if (R_end <= s)
    {
        return true;
    }

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }
    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    for (ui i = 0; i < new_size; ++i)
        for (ui j = i + 1; j < new_size; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj && mf.SAP(new_size + i, j, new_size, real_minCut) < R_end - s)
                return false;
        }
    return true;
}

bool MaximumFlow::verify_Sbundle_by_maxFlowAlg2(vector<ui> &SR, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, int s, ui to_add_vertex, ui lb)
{
    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    ui new_size = 0;
    ui record_to_add_vertex = -1;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        if (u == to_add_vertex)
            record_to_add_vertex = new_size;
        nV[u] = new_size, oID[new_size++] = u;
    }
    assert(record_to_add_vertex != -1);
    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    for (ui i = 0; i < new_size; ++i)
    {
        ui u = oID[i];
        for (ui j = i + 1; j < new_size; ++j)
        {
            ui v = oID[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }

            if (!is_adj && mf.SAP(new_size + i, j, new_size, real_minCut) + s < R_end)
                return false;
        }
    }

    return true;
}

bool MaximumFlow::verify_Sbundle_by_maxFlowAlg3(vector<ui> SR, ui R_end, ui *pstart, ui *pend, ui *edges, int s)
{
    if (R_end <= s)
    {
        return true;
    }

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }

    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    for (ui i = 0; i < new_size; ++i)
        for (ui j = i + 1; j < new_size; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj && mf.SAP(new_size + i, j, new_size, real_minCut) < R_end - s)
                return false;
        }
    return true;
}

bool MaximumFlow::is_connected_graph(vector<ui> &SR, ui R_end, char *matrix, long long row)
{
    if (R_end == 0)
    {
        return true;
    }

    vector<bool> visited(R_end, false);
    queue<ui> queue_for_BFS;

    queue_for_BFS.push(SR[0]);
    visited[0] = true;

    while (!queue_for_BFS.empty())
    {
        ui current_vertex = queue_for_BFS.front();
        queue_for_BFS.pop();
        for (ui i = 0; i < R_end; ++i)
        {
            if (SR[i] != current_vertex && !visited[i] && matrix[current_vertex * row + SR[i]])
            {
                visited[i] = true;
                queue_for_BFS.push(SR[i]);
            }
        }
    }

    for (ui i = 0; i < R_end; ++i)
    {
        if (!visited[i])
        {
            return false;
        }
    }
    return true;
}

void MaximumFlow::get_branching_case(vector<ui> &SR, ui &u1, ui &u2, vector<char> &s_connected, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, ui &min_vc, int s, ui n)
{

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }
    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }
    assert(new_size == R_end);

    for (ui i = 0; i < S_end; ++i)
    {
        for (ui j = i + 1; j < S_end; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < R_end - s)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    min_vc = temp_vc;
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    s_connected = value;
                    assert(!s_connected.empty());
                    value.clear();
                    return;
                }
            }
        }
    }

    for (ui i = 0; i < S_end; ++i)
    {
        for (ui j = S_end; j < R_end; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < R_end - s)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    min_vc = temp_vc;
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    s_connected = value;
                    assert(!s_connected.empty());
                    value.clear();
                    return;
                }
            }
        }
    }

    for (ui i = S_end; i < R_end; ++i)
    {
        for (ui j = i + 1; j < R_end; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < R_end - s)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    min_vc = temp_vc;
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    s_connected = value;
                    assert(!s_connected.empty());
                    value.clear();
                    return;
                }
            }
        }
    }
}

void MaximumFlow::delete_vertices_based_on_lower_bound(vector<ui> &SR, ui S_end, ui R_end, char *matrix, long long row, int s, ui lower_bound, vector<ui> &results)
{

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }

#ifndef NDEBUG

    for (ui i = 0; i < R_end; i++)
    {
        assert(oID[i] == SR[i]);
    }
#endif

    mf.init(new_size + new_size);

    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }

    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        char *t_matrix = matrix + u * row;
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];
            if (!t_matrix[v])
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    assert(new_size == R_end);
    for (ui i = S_end; i < R_end; i++)
    {
        ui u = oID[i];
        char *t_matrix = matrix + u * row;
        assert(u == SR[i]);
        for (ui j = 0; j < S_end; j++)
        {
            ui v = oID[j];
            if (!t_matrix[v] && mf.SAP(new_size + i, j, new_size, real_minCut) + s < lower_bound + 1)
            {
                results.push_back(u);
                break;
            }
        }
    }
}

void MaximumFlow::delete_vertices_based_on_lower_bound(vector<ui> &SR, ui S_end, ui R_end, ept *pstart, ept *pend, ui *edges, ui s, ui lower_bound, vector<ui> &results)
{

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }

#ifndef NDEBUG

    for (ui i = 0; i < R_end; i++)
    {
        assert(oID[i] == SR[i]);
    }
#endif

    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }

    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }

    assert(new_size == R_end);
    for (ui i = S_end; i < R_end; i++)
    {
        ui u = oID[i];
        assert(u == SR[i]);
        for (ui j = 0; j < S_end; j++)
        {
            ui v = oID[j];
            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }

            if (!is_adj && mf.SAP(new_size + i, j, new_size, real_minCut) + s < lower_bound + 1)
            {
                results.push_back(u);
                break;
            }
        }
    }
}

void MaximumFlow::fill_all_s_connected_in_S(vector<ui> &SR, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, map<pair<ui, ui>, vector<char>> &all_s_connected_new, int s, ui n)
{

    MaximumFlow mf;
    mf.reserve(n, m);
    nV = new ui[n + 3];
    oID = new ui[n + 3];

    int new_size = 0;
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        nV[u] = new_size, oID[new_size++] = u;
    }
    mf.init(new_size + new_size);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        mf.addedge(nV[u], nV[u] + new_size, 1);
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui u = SR[i];
        for (ui j = i + 1; j < R_end; j++)
        {
            ui v = SR[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[u] + new_size, nV[v], INF);
            mf.addedge(nV[v] + new_size, nV[u], INF);
        }
    }
    assert(new_size == R_end);

    for (ui i = 0; i < S_end; ++i)
    {
        ui u = oID[i];
        for (ui j = i + 1; j < S_end; ++j)
        {
            ui v = oID[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                assert(!is_adj);
                ui temp_u = oID[i], temp_v = oID[j];
                if (temp_u > temp_v)
                {
                    ui ttemp = temp_u;
                    temp_u = temp_v;
                    temp_v = ttemp;
                }

                assert(temp_u < temp_v);
                pair<ui, ui> key = make_pair(temp_u, temp_v);
                vector<char> value(n, 0);
                for (ui k = 0; k < real_minCut.size(); k++)
                {
                    assert(oID[real_minCut[k]] < n);
                    value[oID[real_minCut[k]]] = 1;
                }
                all_s_connected_new[key] = value;
                value.clear();
            }
        }
    }
}