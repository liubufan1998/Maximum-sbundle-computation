#include "Graph.h"
#include "SBundle_BB_matrix.h"
#include "SBundle_BB.h"
#include "CTPrune.h"
#include <fstream>

using namespace std;
ui pre_n; // record the number of vertices
ui pre_m; // record the number of edges

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

void Graph::read_binary_file()
{
	ifstream in(dir);
	if (!in.is_open())
	{
		printf("Failed to open %s \n", dir.c_str());
		fflush(stdout);
		exit(1);
	}
	string suffix = "bin"; // Note that here the function only accepts the data file in binary form; if not, please use our provided transform function!
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
		ui *degree = new ui[n];
		if (pstart == nullptr)
			pstart = new ept[n + 1];
		if (edges == nullptr)
			edges = new ui[m];
		fread(degree, sizeof(ui), n, in);
		fread(edges, sizeof(ui), m, in);
		pstart[0] = 0;
		for (ui i = 1; i <= n; i++)
			pstart[i] = pstart[i - 1] + degree[i - 1];
		delete[] degree;
	}

	printf("Graph init ok\n");
	fflush(stdout);
	in.close();
}

std::string Graph::get_file_name_without_suffix(const std::string &file_path)
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
	std::vector<ui> to_test = {}; // Please enter the original vertex number of the s-bundle to be detected into the to_test
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

	std::vector<ui> forFlowArr;
	for (ui i = 0; i < n; i++)
	{
		forFlowArr.push_back(i); 
	}
	MaximumFlow *vc = new MaximumFlow(); 
	vc->prepare(n, n * n);				 
	printf("now the S is %u\n", S);
	printf("Now the size of forFlowArr is %u\n", forFlowArr.size());
	if (vc->verify_SBundle_by_MaxFlowAlg(forFlowArr, forFlowArr.size(), matrix, n, S, true))
	{
		printf("Congratulations!!! The instance is a qualified %u-bundle!\n", S);
	}
	else
	{
		printf("Sorry!!! But the instance is not a qualified %u-bundle.\n", S);
	}

	if (matrix != nullptr)
	{
		delete[] matrix;
		matrix = nullptr;
	}
	if (vc != nullptr)
	{
		delete vc;
		vc = nullptr;
	}
}

void Graph::MSBP_SymBD_H()
{
	Timer t;
	long long branch_node = 0; // Record the generated branches during BnB

	assert(S > 0);
	if (S <= 1)
	{
		printf("\tFor k <= 1, please invoke clique computation algorithms\n");
		return;
	}

	if (S >= n)
	{
		printf("\tMaximum sBundle Size: %u, Total Time: %s (microseconds)\n", n, Utility::integer_to_string(t.elapsed()).c_str());
		return;
	}

	// Randomly generate an initial s-bundle of size S for subsequent operations
	sbundle.clear(); 
	for (int i = 0; i < S; i++)
	{
		sbundle.push_back(i);
	}

	pre_n = n; 
	pre_m = m;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

	/*****************The first stage of the heuristic solution: based on degeneracy ordering*********************/
	Timer record_preprocess; // A variable used to record the pre-processing time
	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);
	assert(sbundle.size() >= S);

	if (sbundle.size() < UB)
	{
		ui old_size = sbundle.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n]; 

		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);

		if (sbundle.size() + 1 > 2 * S)
		{
			CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
		}

		/***************The second stage of the heuristic solution: degeneracy ordering based on 1-hop subgraphs****************/
#ifdef _HEURISTIC_SECOND_STAGE_
		pre_n = n; 
		pre_m = m;
		// printf("Before degen2, now the size of the max %d-bundle is %d\n", S, sbundle.size());
		Timer t1;																	
		degen2(n, m, peel_sequence, pstart, edges, degree, rid, vis, heap, true);
		// printf("After degen2, now the size of the max %d-bundle is %d\n", S, sbundle.size());

		if (sbundle.size() > old_size) 
		{
			old_size = sbundle.size();
			for (ui i = 0; i < sbundle.size(); i++)
			{
				assert(sbundle[i] < n);
				sbundle[i] = out_mapping[sbundle[i]]; 
			}
			if (sbundle.size() + 1 > 2 * S)
				CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
			else
				core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
		}
		long long time_ego_degen = t1.elapsed(); 
#endif

		/***********The third stage of the heuristic solution: degeneracy ordering based on 2-hop subgraphs*************/
#ifdef _HEURISTIC_THIRD_STAGE_
		vector<ui> ids;
		vector<pair<ui, ui>> vp;

		ui *peel_sequence_rid = new ui[n]; 
		for (ui i = 0; i < n; i++)
			peel_sequence_rid[peel_sequence[i]] = i; 
		memset(vis, 0, sizeof(char) * n);

		SBUNDLE_BB_matrix *sbundle_solver_m = new SBUNDLE_BB_matrix(); 
		sbundle_solver_m->allocateMemory(n);

		if (pend == nullptr)
			pend = new ept[n + 1];

		reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
		Timer t2;
		for (ui i = n; i > 0 && sbundle.size() < UB; i--)
		{
			long long time_IE_2hop = t2.elapsed();
			// if(time_IE_2hop > 10*time_ego_degen){
			if (time_IE_2hop > 1000000) // If the heuristic takes longer than the time specified by the user, it is terminated
			{ 
				printf("2-hop-Degen heuristic timed out！！\n");
				break;
			}
			ui u = peel_sequence[i - 1];
			if (pend[u] - pstart[u] + S <= sbundle.size() || n - i < sbundle.size())
				continue;

			fflush(stdout);
			if (sbundle.size() >= 2 * S - 1)
				extract_subgraph_with_prune(u, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, sbundle.size() + 3 - 2 * S, peel_sequence_rid, degree, ids, rid, vp, vis, pstart, pend, edges);
			else
				extract_subgraph_wo_prune(u, peel_sequence_rid, ids, rid, vp, vis, pstart, pend, edges);

			if (ids.empty() || ids.size() <= sbundle.size())
				continue;

			ui t_old_size = sbundle.size();
			sbundle_solver_m->load_graph_heuristic(ids.size(), vp);
			sbundle_solver_m->sBundle_heuristic(S, sbundle, true, vp.size(), branch_node); 
			if (sbundle.size() > t_old_size)											 
			{
				for (ui j = 0; j < sbundle.size(); j++)
					sbundle[j] = ids[sbundle[j]]; 

				for (ui i = 0; i < sbundle.size(); i++)
				{
					sbundle[i] = out_mapping[sbundle[i]]; 
				}
			}
		}
		delete sbundle_solver_m;
		sbundle_solver_m = nullptr;
		delete[] peel_sequence_rid;
		peel_sequence_rid = nullptr;
		assert(branch_node == 0); 
		printf("*** 2-hop-Degen sbundle size: %lu, Time: %s (microseconds)\n", sbundle.size(), Utility::integer_to_string(t2.elapsed()).c_str());

		if (sbundle.size() > old_size) 
		{
			old_size = sbundle.size();
			if (sbundle.size() + 1 > 2 * S)
				CTPrune::core_truss_copruning(n, m, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
			else
				core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
		}
#endif

		/***************Now to start solving for the exact maximum s-bundle*****************/
		Timer tt;
		SBundle_BB *sbundle_solver = new SBundle_BB();
		// printf("Now the size of n is %u , the size of m is %u\n", n, m / 2);
		// printf("Now the size of lb is %u\n", sbundle.size());
		printf("*** The time for preprocessing is %s (microseconds)\n", Utility::integer_to_string(record_preprocess.elapsed()).c_str());
		sbundle_solver->allocateMemory(n);
		{
#ifdef _IE_FRAMEWORK_
			vector<ui> ids;
			vector<pair<ui, ui>> vp;

			ui *peel_sequence_rid = new ui[n]; 
			for (ui i = 0; i < n; i++)
				peel_sequence_rid[peel_sequence[i]] = i; 

			memset(vis, 0, sizeof(char) * n);

			SBUNDLE_BB_matrix *sbundle_solver_m = new SBUNDLE_BB_matrix(); 
			sbundle_solver_m->allocateMemory(n);

			ui search_cnt = 0;
			double min_density = 1, total_density = 0;

			if (pend == nullptr)
				pend = new ept[n + 1];

			reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges); 
			for (ui i = n; i > 0 && sbundle.size() < UB; i--) 
			{
				ui u = peel_sequence[i - 1]; 

				if (pend[u] - pstart[u] + S <= sbundle.size() || n - i < sbundle.size()) 
					continue;

				fflush(stdout);
				
				if (sbundle.size() >= 2 * S - 1) 
					extract_subgraph_with_prune(u, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, sbundle.size() + 3 - 2 * S, peel_sequence_rid, degree, ids, rid, vp, vis, pstart, pend, edges);
				else
					extract_subgraph_wo_prune(u, peel_sequence_rid, ids, rid, vp, vis, pstart, pend, edges);

				if (ids.empty() || ids.size() <= sbundle.size())
					continue;

				double density = vp.size() * 2 / (double)ids.size() / (ids.size() - 1);
				++search_cnt;
				total_density += density;
				if (density < min_density)
					min_density = density;

				ui t_old_size = sbundle.size();
				sbundle_solver_m->load_graph(ids.size(), vp);
				sbundle_solver_m->sBundle(S, sbundle.size() + 1, sbundle, true, vp.size(), branch_node);
				if (sbundle.size() > t_old_size)														
				{
					for (ui j = 0; j < sbundle.size(); j++)
						sbundle[j] = ids[sbundle[j]];
				}
			}
			delete sbundle_solver_m;
			delete[] peel_sequence_rid;
			peel_sequence_rid = nullptr;

			if (search_cnt == 0)
				printf("search_cnt: 0, ave_density: 1, min_density: 1\n");
			else
				printf("search_cnt: %u, ave_density: %.5lf, min_density: %.5lf\n", search_cnt, total_density / search_cnt, min_density);
#endif

			/************If the size of the current s-bundle is not more than 2s-2****************/
#ifdef _IE_FRAMEWORK_
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
					printf("!!! WA in kPlex_exact!\n");
			}
#else
			if (n > sbundle.size() && UB > sbundle.size()) 
			{
				sbundle_solver->load_graph(n, pstart, pstart + 1, edges);
				sbundle_solver->sBundle(S, UB, sbundle, false, branch_node);
			}
#endif

			delete sbundle_solver;
			if (sbundle.size() > old_size) // Update vertex index
			{
				for (ui i = 0; i < sbundle.size(); i++)
				{
					sbundle[i] = out_mapping[sbundle[i]];
				}
			}
			printf("Branch_nodes: %lld\n", branch_node);
			printf("*** Search time: %s\n", Utility::integer_to_string(tt.elapsed()).c_str());

			// choose to print the obtained maximum sbundle
			// std::sort(sbundle.begin(), sbundle.end()); 
			// for (ui num : sbundle)
			// {
			// 	printf("%u, ", num);
			// }
			// printf("\n");
		}

		delete[] out_mapping;
		delete[] rid;
	}

	delete heap;
	delete[] core;
	delete[] peel_sequence;
	delete[] vis;
	delete[] degree;

	printf("\tMaximum %d-bundle Size: %lu, Total Time: %s (microseconds)\n", S, sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges)
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

// each of u's neighbors must have at least triangle_threshold common neighbors with u
// each of u's non-neighbor must have at least cn_threshold common neighbors with u
// after pruning, u must have at least degree_threshold neighbors
void Graph::extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
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

	ui *Q = rid;
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

void Graph::extract_subgraph_wo_prune(ui u, const ui *p_rid, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
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

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
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

void Graph::extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
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

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer)
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
		// printf("Entire list: ");
		// for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
		// printf("\n");
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

	// for(ui i = 0;i < ids_n;i ++) printf(" %u", ids[i]);
	// printf("\n");

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

// degeneracy-based s-bundle and return an upper bound of the largest s-bundle
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output)
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
			{
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
			}
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
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
		ui max_core = 0; 
		ui idx = n;	// The index of the heuristic sbundle	 
		ui idx_splex = n; // The index of the heuristic splex
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
				{
					idx_splex = i;
				}
				ui u2 = n, second_min = n;	  
				heap->get_min(u2, second_min); 
				assert(u2 != n && second_min != n);
				if (key + second_min >= 2 * (new_size - i) - S - 2)
				{
					idx = i; 
				}
			}

			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 0) 
				{
					heap->decrement(edges[j], 1);
				}
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
				std::vector<ui> test_sbundle;
				while (new_size - idx_splex > sbundle.size())
				{
					test_sbundle.clear();
					for (ui t = idx_splex; t < new_size; t++)
					{
						test_sbundle.push_back(peel_sequence[queue_n + t]);
					}
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
					{
						idx_splex++;
					}
				}

				if (vc2 != nullptr)
				{
					delete vc2;
					vc2 = nullptr;
				}
			}
		}
		if (output)
			printf("*** Degeneracy s-bundle size: %u, max_core: %u, UB: %u, Time: %s (microseconds)\n", new_size - idx, max_core, UB, Utility::integer_to_string(t.elapsed()).c_str());
#endif
	}
	return UB;
}

void Graph::degen2(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output)
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
	ui *pend_s = rid;
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

		ui threshold = (sbundle.size() > S ? sbundle.size() - S : 0); // all vertices with degree < threshold can be pruned
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

		heap->init(vs.size(), vs.size() - 1, vs.data(), degree); 
		bool found = false;

		ui idx_sbundle = n;
		std::vector<ui> t_degeneracy_ordering;
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
			std::vector<ui> sbundle_copy;
			for (ui i = 1; i < sbundle.size(); i++) 
			{
				sbundle_copy.push_back(sbundle[i]);
			}
			MaximumFlow *vc2 = new MaximumFlow();
			vc2->prepare(pre_n, pre_m);
			ui t_R_end = vs.size() - idx_sbundle;					
			for (ui kk = t_R_end - 1; kk >= 0 && kk < vs.size(); kk--) 
			{
				sbundle_copy.push_back(t_degeneracy_ordering[kk]); 
				bool can_add_v = vc2->verify_Sbundle_by_maxFlowAlg3(sbundle_copy, sbundle_copy.size(), pstart_s, pend_s, edges_s, S);
				if (!can_add_v)
				{
					sbundle_copy.pop_back(); 
				}
				
				if (sbundle_copy.size() + 1 + kk < sbundle.size())
				{
					break; 
				}
			}

			if (sbundle_copy.size() + 1 > sbundle.size()) 
			{
				sbundle.clear();
				sbundle = sbundle_copy;
				sbundle.push_back(u); 
				printf("Now we find a bigger initial sbundle! It's size is %u\n", sbundle_copy.size() + 1);
			}

			if (vc2 != nullptr)
			{
				delete vc2;
				vc2 = nullptr;
			}
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
		printf("*** 1-hop subgraph degen s-bundle size: %lu, Time: %s (microseconds)\n", sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// in_mapping and out_mapping can be the same array
void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output)
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

// orient graph
void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid)
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

// oriented triangle counting
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
	// printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

// reorganize the adjacency lists and sort each adjacency list to be in increasing order
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
