#include "Graph.h"
#include "SBundle_BB_matrix.h"
#include "SBundle_BB.h"
#include "CTPrune.h"
#include <fstream>

using namespace std;
ui pre_n;
ui pre_m;

Graph::Graph(const char *_dir, const int _S)
{
	dir = string(_dir);
	S = _S;
	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	edgelist_pointer = nullptr;

	sbundle.clear();
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != nullptr)
	{
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (edgelist_pointer != nullptr)
	{
		delete[] edgelist_pointer;
		edgelist_pointer = nullptr;
	}
}

void Graph::readBinFile()
{
	ifstream in(dir);
	if (!in.is_open())
	{
		printf("Failed to open %s \n", dir.c_str());
		fflush(stdout);
		exit(1);
	}
	string suffix = "bin";
	if (suffix == "bin")
	{
		FILE *in = fopen(dir.c_str(), "rb");
		if (in == nullptr)
		{
			printf("Failed to open %s \n", dir.c_str());
			exit(1);
		}
		ui size_int;
		fread(&size_int, sizeof(ui), 1, in);
		if (size_int != sizeof(ui))
		{
			printf("sizeof int is different: graph_file(%u), machine(%u)\n", size_int, (int)sizeof(ui));
			exit(1);
		}
		fread(&n, sizeof(ui), 1, in);
		fread(&m, sizeof(ui), 1, in);
		cout << "File: " << get_file_name_without_suffix(dir) << " n= " << n << " m= " << m / 2 << " s = " << S << endl;
		vector<ui> degree(n);
		if (pstart == nullptr)
			pstart = new ept[n + 1];
		if (edges == nullptr)
			edges = new ui[m];
		fread(degree.data(), sizeof(ui), n, in);
		fread(edges, sizeof(ui), m, in);
		pstart[0] = 0;
		for (ui i = 1; i <= n; i++)
			pstart[i] = pstart[i - 1] + degree[i - 1];
	}

	printf("Graph init ok\n");
	fflush(stdout);
	in.close();
}

string Graph::get_file_name_without_suffix(const string &file_path)
{
	size_t last_slash = file_path.find_last_of("/\\");
	size_t last_dot = file_path.find_last_of(".");
	return file_path.substr(last_slash + 1, last_dot - last_slash - 1);
}

void Graph::output_one_sbundle()
{
	FILE *fout = Utility::open_file("sbundles.txt", "w");
	fprintf(fout, "%lu\n", sbundle.size());
	sort(sbundle.begin(), sbundle.end());
	for (ui i = 0; i < sbundle.size(); i++)
		fprintf(fout, " %u", sbundle[i]);
	fprintf(fout, "\n");
	fclose(fout);
}

void Graph::verify_sbundle()
{
	vector<ui> to_test = {
		253,
		322,
		342,
		381,
		607,
		1097,
		5852,
		6098,
		6250,
		9969,
		12612,
		14416,
		16987,
		19010,
		25814,
		27537,
		33547,
		33622,
		33902,
		33984,
		34124,
		34323,
		34737,
		57743,
	};
	ui n = to_test.size();

	char *matrix = new char[n * n];
	memset(matrix, 0, sizeof(char) * n * n);
	for (ui i = 0; i < n; i++)
	{
		ui t_u = to_test[i];
		for (ui j = pstart[t_u]; j < pstart[t_u + 1]; j++)
		{
			for (ui k = i + 1; k < n; k++)
			{
				if (edges[j] == to_test[k])
				{
					matrix[k * n + i] = matrix[i * n + k] = 1;
				}
			}
		}
	}

	vector<ui> forFlowArr;
	for (ui i = 0; i < n; i++)
	{
		forFlowArr.push_back(i);
	}
	MaximumFlow *vc = new MaximumFlow();
	vc->prepare(n, n * n);

	printf("now the S is %u\n", S);
	printf("Now the size of forFlowArr is %u\n", forFlowArr.size());
	if (vc->verify_SBundle_by_MaxFlowAlg0(forFlowArr, forFlowArr.size(), matrix, n, S, true))
		printf("Congratulation! The instance is a qualified %u-bundle!\n", S);
	else
		printf("Sorry, but the instance is not a qualified %u-bundle.\n", S);

	delete[] matrix;
	delete vc;
}

void Graph::SymBD_Parallel(const int number_parallel_threads)
{
	if (S <= 1 or S >= n)
	{
		printf("\tIncorrect parameter for s! s should be in [2, n - 1] \n");
		return;
	}

	Timer t;
	long long branch_node = 0;
	sbundle.resize(S);
	pre_n = n;
	pre_m = m;
	vector<ui> peel_sequence(n, 0);
	vector<ui> core(n, 0);
	vector<ui> degree(n, 0);
	vector<char> vis(n, 0);
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis.data(), heap, false);
	if (sbundle.size() < UB)
	{
		ui old_size = sbundle.size();
		vector<ui> out_mapping(n, 0);
		vector<ui> rid(n);

		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
		if (sbundle.size() + 1 > 2 * S)
			CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, false);

#ifdef _HEURISTIC_SECOND_STAGE_
		pre_n = n;
		pre_m = m;
		ego_degen(n, m, peel_sequence, pstart, edges, degree, rid, vis.data(), heap, true);
		if (sbundle.size() > old_size)
		{
			old_size = sbundle.size();
			for (ui i = 0; i < sbundle.size(); i++)
			{
				assert(sbundle[i] < n);
				sbundle[i] = out_mapping[sbundle[i]];
			}
			if (sbundle.size() + 1 > 2 * S)
				CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, false);
			else
				core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, false);
		}
#endif

#ifdef _HEURISTIC_THIRD_STAGE_
		vector<ui> peel_sequence_rid(n);
		for (ui i = 0; i < n; i++)
			peel_sequence_rid[peel_sequence[i]] = i;

		if (pend == nullptr)
			pend = new ept[n + 1];
		vis.clear();

		reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);

		omp_set_num_threads(number_parallel_threads);
		vector<double> thread_times(number_parallel_threads, 0.0);
		vector<vector<char>> vis_p(number_parallel_threads, vector<char>(n, 0));
		vector<vector<ui>> rid_p(number_parallel_threads, vector<ui>(n, 0));
		vector<vector<ui>> degree_p(number_parallel_threads, vector<ui>(n, 0));
		vector<vector<ui>> ids_p(number_parallel_threads);
		vector<vector<pair<ui, ui>>> vp_p(number_parallel_threads);
		SBUNDLE_BB_matrix **sbundle_solver_m_p = new SBUNDLE_BB_matrix *[number_parallel_threads];
		for (ui i = 0; i < number_parallel_threads; ++i)
		{
			sbundle_solver_m_p[i] = new SBUNDLE_BB_matrix();
			sbundle_solver_m_p[i]->allocateMemory(n);
		}

		Timer t3;

#pragma omp parallel for schedule(static, 1)
		for (ui i = n; i > 0; i--)
		{
			ui u = peel_sequence[i - 1];
			if (pend[u] - pstart[u] + S <= sbundle.size() || n - i < sbundle.size() || sbundle.size() == UB)
				continue;

			fflush(stdout);
			ui thread_id = omp_get_thread_num();
			double start_time = omp_get_wtime();

			assert(thread_id < number_parallel_threads);
			ui t_sbundle_size = sbundle.size();
			if (t_sbundle_size >= 2 * S - 1)
				extract_subgraph_with_prune(u, t_sbundle_size + 1 - S, t_sbundle_size + 1 - 2 * S, t_sbundle_size + 3 - 2 * S, peel_sequence_rid, degree_p[thread_id], ids_p[thread_id], rid_p[thread_id], vp_p[thread_id], vis_p[thread_id], pstart, pend, edges);
			else
				extract_subgraph_wo_prune(u, peel_sequence_rid, ids_p[thread_id], rid_p[thread_id], vp_p[thread_id], vis_p[thread_id], pstart, pend, edges);

			if (ids_p[thread_id].empty() || ids_p[thread_id].size() <= sbundle.size())
			{
				double end_time = omp_get_wtime();
#pragma omp atomic
				thread_times[thread_id] += (end_time - start_time);
				continue;
			}

			sbundle_solver_m_p[thread_id]->load_graph_heuristic(ids_p[thread_id].size(), vp_p[thread_id], degree_p[thread_id]);
			sbundle_solver_m_p[thread_id]->sBundle_heuristic(S, sbundle, true, vp_p[thread_id].size(), branch_node, degree_p[thread_id]);
			double end_time = omp_get_wtime();

#pragma omp atomic
			thread_times[thread_id] += (end_time - start_time);
		}

		assert(branch_node == 0);
		printf("The time spent on the parallel aspect I is: %s (microseconds)\n", Utility::integer_to_string(t3.elapsed()).c_str());

		if (sbundle.size() > old_size)
		{
			old_size = sbundle.size();
			if (sbundle.size() + 1 > 2 * S)
				CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, false);
			else
				core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, false);
		}

		for (ui i = 0; i < n; i++)
			peel_sequence_rid[peel_sequence[i]] = i;
		reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
		vis.clear();

		Timer t4;
#pragma omp parallel for schedule(static, 1)
		for (ui i = n; i > 0; i--)
		{
			ui u = peel_sequence[i - 1];
			if (pend[u] - pstart[u] + S <= sbundle.size() || n - i < sbundle.size() || sbundle.size() == UB)
				continue;

			fflush(stdout);
			ui thread_id = omp_get_thread_num();
			double start_time = omp_get_wtime();

			assert(thread_id < number_parallel_threads);
			ui t_sbundle_size = sbundle.size();
			if (t_sbundle_size >= 2 * S - 1)
				extract_subgraph_with_prune(u, t_sbundle_size + 1 - S, t_sbundle_size + 1 - 2 * S, t_sbundle_size + 3 - 2 * S, peel_sequence_rid, degree_p[thread_id], ids_p[thread_id], rid_p[thread_id], vp_p[thread_id], vis_p[thread_id], pstart, pend, edges);
			else
				extract_subgraph_wo_prune(u, peel_sequence_rid, ids_p[thread_id], rid_p[thread_id], vp_p[thread_id], vis_p[thread_id], pstart, pend, edges);

			bool can_continue = true;
			if (ids_p[thread_id].empty() || ids_p[thread_id].size() <= sbundle.size())
				can_continue = false;

			if (can_continue)
			{
				sbundle_solver_m_p[thread_id]->load_graph(ids_p[thread_id].size(), vp_p[thread_id], degree_p[thread_id]);
				sbundle_solver_m_p[thread_id]->sBundle(S, UB, sbundle, true, vp_p[thread_id].size(), branch_node, degree_p[thread_id]);
			}

			double end_time = omp_get_wtime();
#pragma omp atomic
			thread_times[thread_id] += (end_time - start_time);
		}

		printf("The time spent on the parallel aspect II is: %s (microseconds)\n", Utility::integer_to_string(t4.elapsed()).c_str());
		for (ui i = 0; i < number_parallel_threads; ++i)
			delete sbundle_solver_m_p[i];
		delete[] sbundle_solver_m_p;

		SBundle_BB *sbundle_solver = new SBundle_BB();
		sbundle_solver->allocateMemory(n);
		if (n > sbundle.size() && UB > sbundle.size() && sbundle.size() < 2 * S - 2)
		{
			if (sbundle.size() > old_size)
			{
				old_size = sbundle.size();
				core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
			}
			sbundle_solver->load_graph(n, pstart, pstart + 1, edges);
			if (2 * S - 2 < UB)
				UB = 2 * S - 2;
			sbundle_solver->sBundle(S, UB, sbundle, false, branch_node);
			if (sbundle.size() > 2 * S - 2)
				printf("!!! WA in sBundle_exact!\n");
		}
		delete sbundle_solver;
#endif
	}

	delete heap;

	printf("\tMaximum %d-bundle Size: %lu, Total Time: %s (microseconds)\n", S, sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, vector<ui> &peel_sequence, vector<ui> &rid, ui *pstart, ui *pend, ui *edges)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;

	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		for (ui j = pstart[u]; j < pend[u] && rid[edges[j]] > rid[u]; j++)
		{
			ui v = edges[j];
			edges[pend[v]++] = u;
			assert(pend[v] <= pstart[v + 1]);
		}
	}
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		while (end < pstart[i + 1] && rid[edges[end]] > rid[i])
			++end;
	}
}

void Graph::extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const vector<ui> &p_rid, vector<ui> &degree, vector<ui> &ids, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(!exists[i]);
#endif

	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);

	ui *Q = rid.data();
	ui Q_n = 0;
	for (ui i = 1; i < ids.size(); i++)
	{
		ui v = ids[i];
		degree[v] = 0;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			if (exists[edges[j]])
				++degree[v];
		}
		if (degree[v] < triangle_threshold)
			Q[Q_n++] = v;
	}

	for (ui i = 0; i < Q_n; i++)
	{
		ui v = Q[i];
		exists[v] = 3;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if (degree[edges[j]] == triangle_threshold)
					Q[Q_n++] = edges[j];
				--degree[edges[j]];
			}
	}

	assert(Q_n < ids.size());

	if (ids.size() - Q_n - 1 < degree_threshold)
	{
		for (ui i = 0; i < ids.size(); i++)
			exists[ids[i]] = 0;
		ids.clear();
		return;
	}

	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
		if (exists[ids[i]] == 2)
		{
			ui v = ids[i];
			for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids.push_back(edges[j]);
					exists[edges[j]] = 1;
					degree[edges[j]] = 1;
				}
				else
					++degree[edges[j]];
			}
		}

	ui new_size = 1;
	for (ui i = 1; i < old_size; i++)
	{
		if (exists[ids[i]] == 3)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);

	for (ui i = old_size; i < ids.size(); i++)
	{
		if (degree[ids[i]] < cn_threshold)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for (ui i = 0; i < ids.size(); i++)
		rid[ids[i]] = i;

	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
			{
				assert(rid[v] < ids.size() && rid[edges[j]] < ids.size());
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
			}
	}

	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_wo_prune(ui u, const vector<ui> &p_rid, vector<ui> &ids, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(!exists[i]);
#endif
	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	rid[u] = 0;

	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 1;
		rid[edges[i]] = ids.size() - 1;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);
	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			ui w = edges[j];
			if (exists[w])
				continue;
			ids.push_back(w);
			exists[w] = 1;
			rid[w] = ids.size() - 1;
		}
	}
	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}
	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::write_subgraph(ui n, const vector<pair<int, int>> &edge_list)
{
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for (ui i = 0; i < edge_list.size(); i++)
		fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	ids_n = 0;
	vp.clear();
	ids[ids_n++] = u;
	exists[u] = 1;
	rid[u] = 0;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			rid[v] = ids_n;
			ids[ids_n++] = v;
			exists[v] = 1;
		}
	pend[u] = u_n;
	ui old_size = ids_n;
	for (ui i = 1; i < old_size; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				ui v = edges[j];
				if (exists[v])
					continue;
				rid[v] = ids_n;
				ids[ids_n++] = v;
				exists[v] = 1;
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < old_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (edges[j] > u)
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = old_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_full(const ui *ids, ui ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, vector<char> &exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	for (ui i = 0; i < ids_n; i++)
	{
		rid[ids[i]] = i;
		exists[ids[i]] = 1;
	}
	for (ui i = 0; i < ids_n; i++)
	{
		ui u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]] && edges[j] > u && exists[edges[j]])
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, vector<ui> &rid, vector<pair<ui, ui>> &vp, ui *Q, vector<ui> &degree, vector<char> &exists, ept *pend, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	ids_n = 0;
	ids[ids_n++] = u;
	exists[u] = 1;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			ids[ids_n++] = v;
			exists[v] = 2;
		}
	pend[u] = u_n;

	ui Q_n = 0;
	for (ui i = 1; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (exists[edges[j]] == 2)
					++degree[u];
			}
		pend[u] = u_n;
		if (degree[u] + 2 * S <= sbundle.size())
			Q[Q_n++] = u;
	}
	for (ui i = 0; i < Q_n; i++)
	{
		u = Q[i];
		exists[u] = 10;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if ((degree[edges[j]]--) + 2 * S == sbundle.size() + 1)
				{
					assert(Q_n < m / 2);
					Q[Q_n++] = edges[j];
				}
			}
	}
	assert(Q_n <= ids_n);
	if (ids_n - 1 - Q_n + S <= sbundle.size())
	{
		for (ui i = 0; i < ids_n; i++)
			exists[ids[i]] = 0;
		ids_n = 0;
		return;
	}

	ui nr_size = ids_n;
	for (ui i = 1; i < nr_size; i++)
		if (exists[ids[i]] == 2)
		{
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids[ids_n++] = edges[j];
					exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (exists[edges[j]] == 3)
					++degree[edges[j]];
			}
		}

#ifndef NDEBUG

#endif

	ui new_size = 1;
	for (ui i = 1; i < nr_size; i++)
	{
		if (exists[ids[i]] == 10)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if (new_size + Q_n != nr_size)
	{
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for (ui i = 0; i < new_size; i++)
			printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for (ui i = 0; i < Q_n; i++)
			printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	ui old_nr_size = nr_size;
	nr_size = new_size;
	for (ui i = old_nr_size; i < ids_n; i++)
	{
		if (degree[ids[i]] + 2 * S <= sbundle.size() + 2)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for (ui i = 1; i < nr_size; i++)
		assert(exists[ids[i]] == 2);
	for (ui i = nr_size; i < ids_n; i++)
		assert(exists[ids[i]] == 3);
#endif

	for (ui i = 0; i < ids_n; i++)
	{
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for (ui i = 0; i < nr_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] && edges[j] > u)
			{
				assert(!deleted[edgelist_pointer[j]]);
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = nr_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(exists[i] == 0);
#endif
}

ui Graph::degen(ui n, vector<ui> &peel_sequence, vector<ui> &core, ept *pstart, ui *edges, vector<ui> &degree, char *vis, ListLinearHeap *heap, bool output)
{
	Timer t;
	ui threshold = (sbundle.size() + 1 > S ? sbundle.size() + 1 - S : 0);

	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			peel_sequence[queue_n++] = i;

	for (ui i = 0; i < queue_n; i++)
	{
		ui u = peel_sequence[i];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
	}

	ui UB = n;
	if (queue_n == n)
		UB = sbundle.size();

	memset(vis, 0, sizeof(char) * n);
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
			peel_sequence[queue_n + (new_size++)] = i;
		else
		{
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	if (new_size != 0)
	{
		heap->init(new_size, new_size - 1, peel_sequence.data() + queue_n, degree.data());
		ui max_core = 0;
		ui idx = n;
		ui idx_splex = n;
		UB = 0;
		for (ui i = 0; i < new_size; i++)
		{
			ui u, key;
			heap->pop_min(u, key);
			if (key > max_core)
				max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + S;
			if (new_size - i < t_UB)
				t_UB = new_size - i;
			if (t_UB > UB)
				UB = t_UB;

			if (idx == n && key + S >= new_size - i)
			{
				if (idx_splex == n)
					idx_splex = i;

				ui u2 = n, second_min = n;
				heap->get_min(u2, second_min);
				assert(u2 != n && second_min != n);
				if (key + second_min >= 2 * (new_size - i) - S - 2)
					idx = i;
			}

			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 0)
					heap->decrement(edges[j], 1);
		}

		assert(idx_splex <= idx);
#ifdef _HEURISTIC_FIRST_STAGE_
		if (new_size - idx_splex > sbundle.size())
		{
			if (idx_splex == idx)
			{
				sbundle.clear();
				for (ui i = idx; i < new_size; i++)
					sbundle.pb(peel_sequence[queue_n + i]);
				if (!output)
					printf("Find an s-bundle of size: %u\n", new_size - idx);
			}
			else
			{
				assert(idx_splex < idx);
				MaximumFlow *vc2 = new MaximumFlow();
				vc2->prepare(pre_n, pre_m);
				vector<ui> test_sbundle;
				while (new_size - idx_splex > sbundle.size())
				{
					test_sbundle.clear();

					for (ui t = idx_splex; t < new_size; t++)
						test_sbundle.push_back(peel_sequence[queue_n + t]);

					bool is_sbundle = vc2->verify_Sbundle_by_maxFlowAlg3(test_sbundle, test_sbundle.size(), pstart, pstart + 1, edges, S);
					if (is_sbundle)
					{
						assert(test_sbundle.size() > sbundle.size());
						sbundle.clear();
						sbundle = test_sbundle;

						printf("Now we find a larger solution, and its size is %u\n", sbundle.size());
						if (!output)
							printf("Find a s-bundle of size: %u\n", new_size - idx_splex);
						break;
					}
					else
						idx_splex++;
				}

				delete vc2;
			}
		}
		if (output)
			printf("*** Degeneracy s-bundle size: %u, max_core: %u, UB: %u, Time: %s (microseconds)\n", new_size - idx, max_core, UB, Utility::integer_to_string(t.elapsed()).c_str());
#endif
	}
	return UB;
}

void Graph::ego_degen(ui n, ui m, vector<ui> &peel_sequence, ept *pstart, ui *edges, vector<ui> &degree, vector<ui> &rid, char *vis, ListLinearHeap *heap, bool output)
{

	Timer t;
	if (pend == nullptr)
		pend = new ept[n + 1];

	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);

	if (pend_buf == nullptr)
		pend_buf = new ept[n + 1];
	if (edgelist_pointer == nullptr)
		edgelist_pointer = new ui[m];
	ui *pstart_s = pend_buf;
	ui *pend_s = rid.data();
	ui *edges_s = edgelist_pointer;

	vector<ui> Q;
	vector<ui> vs;
	memset(vis, 0, sizeof(char) * n);
	for (ui i = n; i > 0; i--)
	{

		ui u = peel_sequence[i - 1];
		if (pend[u] - pstart[u] < sbundle.size())
			continue;

		vs.clear();
		for (ui j = pstart[u]; j < pend[u]; j++)
		{
			vs.push_back(edges[j]);

			vis[edges[j]] = 1;
			degree[edges[j]] = 0;
		}

		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					++degree[vs[j]];
					++degree[edges[k]];
				}

		pend_s[vs[0]] = pstart_s[vs[0]] = 0;
		for (ui j = 1; j < vs.size(); j++)
			pend_s[vs[j]] = pstart_s[vs[j]] = pstart_s[vs[j - 1]] + degree[vs[j - 1]];

		for (ui j = 0; j < vs.size(); j++)
			for (ui k = pstart[vs[j]]; k < pend[vs[j]]; k++)
				if (vis[edges[k]])
				{
					edges_s[pend_s[vs[j]]++] = edges[k];
					edges_s[pend_s[edges[k]]++] = vs[j];
				}

		ui threshold = (sbundle.size() > S ? sbundle.size() - S : 0);
		Q.clear();

		for (ui j = 0; j < vs.size(); j++)
			if (degree[vs[j]] < threshold)
			{
				Q.push_back(vs[j]);
				vis[vs[j]] = 0;
			}

		for (ui j = 0; j < Q.size(); j++)
			for (ui k = pstart_s[Q[j]]; k < pend_s[Q[j]]; k++)
				if (vis[edges_s[k]])
				{
					if ((degree[edges_s[k]]--) == threshold)
					{
						Q.push_back(edges_s[k]);
						vis[edges_s[k]] = 0;
					}
				}

		ui cnt = 0;
		for (ui j = 0; j < vs.size(); j++)
			if (vis[vs[j]])
				vs[cnt++] = vs[j];
		assert(cnt + Q.size() == vs.size());
		vs.resize(cnt);
		if (cnt == 0)
			continue;

		heap->init(vs.size(), vs.size() - 1, vs.data(), degree.data());
		bool found = false;

		ui idx_sbundle = n;
		vector<ui> t_degeneracy_ordering;

		for (ui ii = 0; ii < vs.size(); ii++)
		{
			ui v, key;
			heap->pop_min(v, key);
			t_degeneracy_ordering.push_back(v);
			if (found)
			{
				sbundle.push_back(v);
				continue;
			}

			if (vs.size() - ii + 1 <= sbundle.size())
				break;

			if (key + S >= vs.size() - ii)
			{
				ui u2 = n, second_min = n;
				heap->get_min(u2, second_min);
				assert(u2 != n && second_min != n);
				if (key + second_min >= 2 * (vs.size() - ii) - S - 2)
				{
					sbundle.clear();
					sbundle.push_back(u);
					sbundle.push_back(v);
					found = true;
					continue;
				}
			}
			vis[v] = 0;
			for (ept j = pstart_s[v]; j < pend_s[v]; j++)
				if (vis[edges_s[j]])
					heap->decrement(edges_s[j], 1);
		}

		if (found)
		{
			assert(t_degeneracy_ordering.size() == vs.size());
			vector<ui> sbundle_copy;
			for (ui i = 1; i < sbundle.size(); i++)
				sbundle_copy.push_back(sbundle[i]);
			MaximumFlow *vc2 = new MaximumFlow();
			vc2->prepare(pre_n, pre_m);
			ui t_R_end = vs.size() - idx_sbundle;
			for (ui kk = t_R_end - 1; kk >= 0 && kk < vs.size(); kk--)
			{
				sbundle_copy.push_back(t_degeneracy_ordering[kk]);
				bool can_add_v = vc2->verify_Sbundle_by_maxFlowAlg3(sbundle_copy, sbundle_copy.size(), pstart_s, pend_s, edges_s, S);
				if (!can_add_v)
					sbundle_copy.pop_back();

				if (sbundle_copy.size() + 1 + kk < sbundle.size())
					break;
			}

			if (sbundle_copy.size() + 1 > sbundle.size())
			{
				sbundle.clear();
				sbundle = sbundle_copy;
				sbundle.push_back(u);
				printf("Now we find a bigger initial sbundle! It's size is %u\n", sbundle_copy.size() + 1);
			}
			delete vc2;
		}
		for (ui j = 0; j < vs.size(); j++)
			vis[vs[j]] = 0;
	}

	for (ui i = 0; i < n; i++)
		pend_buf[i] = pend[i];

	for (ui i = 0; i < n; i++)
		for (ept j = pstart[i]; j < pend[i]; j++)
			edges[pend_buf[edges[j]]++] = i;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend_buf[i] == pstart[i + 1]);
#endif

	if (output)
		printf("*** EGo-Degen s-bundle size: %lu, Time: %s (microseconds)\n", sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::core_shrink_graph(ui &n, ept &m, vector<ui> &peel_sequence, vector<ui> &core, vector<ui> &out_mapping, ui *in_mapping, vector<ui> &rid, ept *&pstart, ui *&edges, bool output)
{

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] + S > sbundle.size())
		{
			rid[i] = cnt;
			if (in_mapping == nullptr)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] + S > sbundle.size())
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] + S > sbundle.size())
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + S <= sbundle.size());
		assert(cnt == 0 || core[peel_sequence[n - cnt]] + S > sbundle.size());

		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			core[i] = core[out_mapping[i]];
		}

		if (pos > 0 && pos < m / 2)
		{

			ept *pstart_new = new ept[cnt + 1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
			memcpy(edges_new, edges, sizeof(ui) * pos);

			delete[] pstart;
			pstart = pstart_new;
			delete[] edges;
			edges = edges_new;
		}
		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

void Graph::orient_graph(ui n, ui m, vector<ui> &peel_sequence, ept *pstart, ept *pend, ui *edges, vector<ui> &rid)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for (int i = 0; i < n; i++)
		sum += pend[i] - pstart[i];
	assert(sum * 2 == m);
#endif
}

void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
{
	memset(adj, 0, sizeof(ui) * n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui) * m);
	for (ui u = 0; u < n; u++)
	{
		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = j + 1;

		for (ept j = pstart[u]; j < pend[u]; j++)
		{
			ui v = edges[j];
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (adj[edges[k]])
				{
					++tri_cnt[j];
					++tri_cnt[k];
					++tri_cnt[adj[edges[k]] - 1];
					++cnt;
				}
		}

		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = 0;
	}

#ifndef NDEBUG

#endif
}

void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf)
{
	for (ui i = 0; i < n; i++)
		pend2[i] = pend[i];
	ept pos = 0;
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pstart[i]; j < pend[i]; j++)
		{
			tri_cnt[pos >> 1] = edgelist_pointer[j];
			edge_list[pos++] = i;
			edge_list[pos++] = edges[j];

			ept &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos >> 1) - 1;
			edges[k++] = i;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend2[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pend2[i]; j < pstart[i + 1]; j++)
		{
			ept &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k++] = i;
		}
	}

	ept *ids = pend2;
	for (ui i = 0; i < n; i++)
	{
		if (pend[i] == pstart[i] || pend[i] == pstart[i + 1])
			continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while (j < pend[i] && k < pstart[i + 1])
		{
			if (edges[j] < edges[k])
			{
				ids[pos] = edges[j];
				buf[pos++] = edgelist_pointer[j++];
			}
			else
			{
				ids[pos] = edges[k];
				buf[pos++] = edgelist_pointer[k++];
			}
		}
		while (j < pend[i])
		{
			ids[pos] = edges[j];
			buf[pos++] = edgelist_pointer[j++];
		}
		while (k < pstart[i + 1])
		{
			ids[pos] = edges[k];
			buf[pos++] = edgelist_pointer[k++];
		}
		for (ept j = 0; j < pos; j++)
		{
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}

char Graph::find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges)
{
	if (b >= e)
		return 0;

	while (b + 1 < e)
	{
		idx = b + (e - b) / 2;
		if (edges[idx] > w)
			e = idx;
		else
			b = idx;
	}

	if (edges[b] == w)
	{
		idx = edgelist_pointer[b];
		if (!deleted[idx])
			return 1;
	}

	return 0;
}

ept Graph::peeling(ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, vector<ui> &degree, ept *pstart, ept *pend, ui *edges, vector<char> &exists)
{
	ept Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	if (initialize_Qe)
	{
		ept active_edgelist_newn = 0;
		for (ept j = 0; j < active_edgelist_n; j++)
			if (!deleted[active_edgelist[j]])
			{
				if (tri_cnt[active_edgelist[j]] < t_threshold)
					Qe[Qe_n++] = active_edgelist[j];
				else
					active_edgelist[active_edgelist_newn++] = active_edgelist[j];
			}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	ept deleted_edges_n = 0;
	ui Qv_idx = 0;
	while (Qv_idx < Qv_n || Qe_n)
	{
		if (Qe_n == 0)
		{

			ui u = Qv[Qv_idx++];
			ept u_n = pstart[u];
			for (ept k = pstart[u]; k < pend[u]; k++)
				if (!deleted[edgelist_pointer[k]])
				{
					edges[u_n] = edges[k];
					edgelist_pointer[u_n++] = edgelist_pointer[k];
					exists[edges[k]] = 1;
				}
			pend[u] = u_n;

			for (ept k = pstart[u]; k < pend[u]; k++)
				deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			linear_heap->decrement(u, degree[u]);
			degree[u] = 0;

			for (ept k = pstart[u]; k < pend[u]; k++)
			{
				ui v = edges[k];
#ifndef NO_TRUSS_PRUNE
				ept v_n = pstart[v];
				for (ept x = pstart[v]; x < pend[v]; x++)
					if (!deleted[edgelist_pointer[x]])
					{
						edges[v_n] = edges[x];
						edgelist_pointer[v_n++] = edgelist_pointer[x];
						if (edges[x] > v && exists[edges[x]])
						{
							if ((tri_cnt[edgelist_pointer[x]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[x];
						}
					}
				pend[v] = v_n;
#endif

				linear_heap->decrement(v, 1);
				if ((degree[v]--) == d_threshold)
					Qv[Qv_n++] = v;
			}

			for (ept k = pstart[u]; k < pend[u]; k++)
				exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for (ept j = 0; j < Qe_n; j++)
		{
			ept idx = Qe[j];
			ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
			ui tri_n = tri_cnt[idx];

			deleted[idx] = 1;
			linear_heap->decrement(u, 1);
			linear_heap->decrement(v, 1);
			if ((degree[u]--) == d_threshold)
				Qv[Qv_n++] = u;
			if ((degree[v]--) == d_threshold)
				Qv[Qv_n++] = v;
			deleted_edges_n++;

			if (degree[u] < degree[v])
				swap(u, v);

			if (degree[u] > degree[v] * 2)
			{

				ept v_n = pstart[v], start = pstart[u];
				for (ept k = pstart[v]; k < pend[v]; k++)
					if (!deleted[edgelist_pointer[k]])
					{
						edges[v_n] = edges[k];
						edgelist_pointer[v_n++] = edgelist_pointer[k];

						if (tri_n && find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges))
						{
							--tri_n;
							if ((tri_cnt[idx]--) == t_threshold)
								Qe[Qe_n++] = idx;
							if ((tri_cnt[edgelist_pointer[k]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[k];
						}
					}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else
			{
				ept ii = pstart[u], jj = pstart[v];
				ept u_n = pstart[u], v_n = pstart[v];

				while (true)
				{
					while (ii < pend[u] && deleted[edgelist_pointer[ii]])
						++ii;
					while (jj < pend[v] && deleted[edgelist_pointer[jj]])
						++jj;
					if (ii >= pend[u] || jj >= pend[v])
						break;

					if (edges[ii] == edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if ((tri_cnt[edgelist_pointer[ii]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[ii];
						if ((tri_cnt[edgelist_pointer[jj]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[jj];

						++ii;
						++jj;
					}
					else if (edges[ii] < edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ii;
					}
					else
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++jj;
					}
				}
				while (ii < pend[u])
				{
					if (!deleted[edgelist_pointer[ii]])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ii;
				}
				while (jj < pend[v])
				{
					if (!deleted[edgelist_pointer[jj]])
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++jj;
				}
				pend[u] = u_n;
				pend[v] = v_n;
			}
		}
		Qe_n = 0;
	}
#ifndef NDEBUG
	printf("*** Truss removed %s undirected edges\n", Utility::integer_to_string(deleted_edges_n).c_str());
#endif
	return deleted_edges_n;
}
