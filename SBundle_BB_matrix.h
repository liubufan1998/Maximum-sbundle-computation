#ifndef _SBUNDLE_BB_MATRIX_
#define _SBUNDLE_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "sbundle_tool.h"

class SBUNDLE_BB_matrix
{
private:
	long long n;
	char *matrix;
	long long matrix_size;
	ui *local_vc;													// The array used for dynamic maintaining of local vertex connectivity
	std::vector<char> s_connected;									// The array used for storing s-connected vertices
	std::vector<ui> temp_deleted;									// Used to save temporarily deleted vertices
	MaximumFlow *vc;												// The maximum flow object for calculating vertec connectivity and assess sbundle
	long long branch_count;											// Counting the overall generated branches during BnB
	std::map<std::pair<ui, ui>, int> all_not_adj_inS;				// For storing non-adjacent vertex pairs in S for caculating our proposed upper bound
	std::map<std::pair<ui, ui>, std::vector<char>> all_s_connected; // Used to hold the set of s-connected vertices for each non-adjacent vertex pair

	ui *degree;
	ui *degree_in_S;

	ui s;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;

	ui *neighbors;
	ui *nonneighbors;
	ui *S2;

	ui *SR;		// union of s and R, where s is at the front
	ui *SR_rid; // reverse ID for SR
	std::queue<ui> Qv;
	ui *level_id;
	ui max_level;

	std::vector<std::pair<ui, ui>> vp;
	std::vector<ui> non_adj;

public:
	SBUNDLE_BB_matrix()
	{
		n = 0;
		matrix = nullptr;
		matrix_size = 0;

		degree = degree_in_S = nullptr;

		best_solution = nullptr;
		s = best_solution_size = _UB_ = 0;

		neighbors = nonneighbors = nullptr;
		S2 = nullptr;

		SR = SR_rid = nullptr;
		level_id = nullptr;
		max_level = 0;
	}

	~SBUNDLE_BB_matrix()
	{
		if (matrix != NULL)
		{
			delete[] matrix;
			matrix = NULL;
		}
		if (degree != NULL)
		{
			delete[] degree;
			degree = NULL;
		}
		if (degree_in_S != NULL)
		{
			delete[] degree_in_S;
			degree_in_S = NULL;
		}
		if (best_solution != NULL)
		{
			delete[] best_solution;
			best_solution = NULL;
		}
		if (SR != NULL)
		{
			delete[] SR;
			SR = NULL;
		}
		if (SR_rid != NULL)
		{
			delete[] SR_rid;
			SR_rid = NULL;
		}
		if (neighbors != NULL)
		{
			delete[] neighbors;
			neighbors = NULL;
		}
		if (nonneighbors != NULL)
		{
			delete[] nonneighbors;
			nonneighbors = NULL;
		}
		if (level_id != NULL)
		{
			delete[] level_id;
			level_id = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;

		matrix_size = 1;
		matrix = new char[matrix_size];
		degree = new ui[n];
		degree_in_S = new ui[n];
		best_solution = new ui[n];
		SR = new ui[n];
		SR_rid = new ui[n];
		neighbors = new ui[n];
		nonneighbors = new ui[n];
		S2 = new ui[n];
		level_id = new ui[n];
	}

	// 
	void load_graph_heuristic(ui _n, const std::vector<std::pair<ui, ui>> &vp)
	{
		n = _n;
		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
		}
		memset(matrix, 0, sizeof(char) * ((long long)n) * n);
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			degree[a]++;
			degree[b]++;
			if (matrix[a * n + b])
				printf("Duplicate edge in SBundle_BB_matrix.load_graph()\n");
			matrix[a * n + b] = matrix[b * n + a] = 1;
		}
	}

	void load_graph(ui _n, const std::vector<std::pair<ui, ui>> &vp)
	{
		n = _n;
		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
		}
		local_vc = new ui[matrix_size];
		memset(matrix, 0, sizeof(char) * ((long long)n) * n);
		memset(local_vc, 0, sizeof(ui) * ((long long)n) * n);
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			degree[a]++;
			degree[b]++;
			if (matrix[a * n + b])
			{
				printf("Duplicate edge in SBundle_BB_matrix.load_graph()\n");
			}

			matrix[a * n + b] = matrix[b * n + a] = 1;
		}
	}

	void sBundle(ui S_, ui UB_, std::vector<ui> &sbundle, bool must_include_0, const int m, long long &branch_node)
	{
		s = S_;
		branch_count = 0;
		_UB_ = UB_;
		if (s <= 1)
		{
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return;
		}
		best_solution_size = sbundle.size();
		ui R_end;
		initialization(R_end, must_include_0, true);
		vc = new MaximumFlow();
		vc->prepare(n, m);
		vc->get_local_vertexConnectivity(R_end, SR, matrix, n, local_vc, all_s_connected, s);

		if (R_end && best_solution_size < _UB_)
			SymBD_H(0, R_end, 1, must_include_0);

		// if (R_end && best_solution_size < _UB_)
		// 	SymBD_L(0, R_end, 1, must_include_0); // Branches based entirely on local vertex connectivity

		branch_node += branch_count;

		if (vc != nullptr)
		{
			delete vc;
			vc = nullptr;
		}

		if (local_vc != nullptr)
		{
			delete[] local_vc;
			local_vc = nullptr;
		}

		all_s_connected.clear();

		if (best_solution_size > sbundle.size()) // Update vertex index
		{
			sbundle.clear();
			for (int i = 0; i < best_solution_size; i++)
				sbundle.push_back(best_solution[i]);
		}
	}

	void sBundle_heuristic(ui S_, std::vector<ui> &sbundle, bool must_include_0, const int m, long long &branch_node)
	{
		s = S_;
		branch_count = 0;
		if (s <= 1)
		{
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return;
		}
		best_solution_size = sbundle.size();
		ui R_end;
		initialization(R_end, must_include_0, false);

		if (best_solution_size > sbundle.size())
		{
			sbundle.clear();
			for (int i = 0; i < best_solution_size; i++)
				sbundle.push_back(best_solution[i]);
		}
	}

	void initialization(ui &R_end, bool must_include_0, bool is_precise)
	{
		// the following computes a degeneracy ordering and a heuristic solution
		ui *peel_sequence = neighbors;
		ui *core = nonneighbors;
		ui *vis = SR;
		memset(vis, 0, sizeof(ui) * n);
		ui max_core = 0, idx = n;
		ui idx_splex = n;
		for (ui i = 0; i < n; i++)
		{
			ui u, min_degree = n;
			ui min2_degree = n; // the second minimum degree in the graph
			for (ui j = 0; j < n; j++)
			{
				if (!vis[j] && degree[j] < min_degree)
				{
					u = j;
					min2_degree = min_degree;
					min_degree = degree[j];
				}
				else if (!vis[j] && degree[j] < min2_degree)
				{
					min2_degree = degree[j];
				}
			}
			if (min_degree > max_core)
				max_core = min_degree;
			core[u] = max_core;
			peel_sequence[i] = u;
			vis[u] = 1;

			if (idx == n && min_degree + s >= n - i)
			{
				if (idx_splex == n)
				{
					idx_splex = i;
				}
				if (min_degree + min2_degree >= 2 * (n - i) - s - 2)
				{
					idx = i;
				}
			}

			for (ui j = 0; j < n; j++)
				if (!vis[j] && matrix[u * n + j])
					--degree[j];
		}

		assert(idx_splex <= idx);

#ifdef _HEURISTIC_THIRD_STAGE_
		// Try to extend the heuristic sbundle into a maximal solution
		if (n - idx_splex > best_solution_size)
		{
			std::vector<ui> sbundle_copy;
			for (ui i = idx; i < n; i++)
			{
				sbundle_copy.push_back(peel_sequence[i]);
			}
			MaximumFlow *vc2 = new MaximumFlow();
			vc2->prepare(n, n * n);
			for (ui kk = idx - 1; kk > 0 && kk < n; kk--)
			{
				assert(kk >= 0 && kk <= n);
				sbundle_copy.push_back(peel_sequence[kk]);
				bool can_add_v = vc2->verify_SBundle_by_MaxFlowAlg(sbundle_copy, sbundle_copy.size(), matrix, n, s, false);
				if (!can_add_v)
				{
					sbundle_copy.pop_back();
				}
				else if (sbundle_copy.size() > best_solution_size)
				{
					break;
				}

				if (sbundle_copy.size() + kk < best_solution_size)
				{
					break;
				}
			}

			if (sbundle_copy.size() > best_solution_size)
			{
				best_solution_size = sbundle_copy.size();
				for (ui i = 0; i < best_solution_size; i++)
					best_solution[i] = sbundle_copy[i];
				printf("2hop-Degen find a larger solution of size %u\n", best_solution_size);
			}

			if (vc2 != nullptr)
			{
				delete vc2;
				vc2 = nullptr;
			}
		}
#endif

		if (is_precise)
		{
			R_end = 0;
			for (ui i = 0; i < n; i++)
				SR_rid[i] = n;
			for (ui i = 0; i < n; i++)
				if (core[i] + s > best_solution_size)
				{
					SR[R_end] = i;
					SR_rid[i] = R_end;
					++R_end;
				}

			if ((must_include_0 && SR_rid[0] == n) || best_solution_size >= _UB_)
			{
				R_end = 0;
				return;
			}

			for (ui i = 0; i < R_end; i++)
			{
				ui u = SR[i];
				degree[u] = degree_in_S[u] = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[u * n + SR[j]])
						++degree[u];
			}

			memset(level_id, 0, sizeof(ui) * n);
			for (ui i = 0; i < R_end; i++)
				level_id[SR[i]] = n;

			if (!Qv.empty())
				printf("!!! Something wrong. Qv must be empty in initialization\n");
		}
	}

	void store_solution(ui size)
	{
		if (size <= best_solution_size)
		{
			printf("!!! the solution to store is no larger than the current best solution!");
			return;
		}
		best_solution_size = size;
		for (ui i = 0; i < best_solution_size; i++)
			best_solution[i] = SR[i];
	}

	bool is_splex(ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
			if (degree[SR[i]] + s < R_end)
				return false;
		return true;
	}

	// First locate the s-plex that meets the requirements and then call our branch method based on local vertex connectivity
	void SymBD_H(ui S_end, ui R_end, ui level, bool choose_zero)
	{

		if (S_end > best_solution_size && vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end, matrix, n, s, true))
			store_solution(S_end);
		if (R_end > best_solution_size && is_splex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		branch_count++;
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
			// Verify that the code for dynamically updating local vertex connectivity is correct
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			char *t_matrix = matrix + u * n;
			for (ui j = i + 1; j < R_end; j++)
			{
				ui v = SR[j];
				if (!t_matrix[v])
				{
					ui cnt = 0;
					std::vector<char> t_s_connected;
					if (u > v)
					{
						t_s_connected = all_s_connected[std::make_pair(v, u)];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(u, v)];
					}
					assert(!t_s_connected.empty());
					for (ui j = 0; j < R_end; j++)
					{
						int w = SR[j];
						if (t_s_connected[w])
							cnt++;
					}
					assert(local_vc[n * v + u] == cnt && local_vc[n * u + v] == cnt);
				}
			}
		}
#endif
		ui old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());

		if (choose_zero && SR_rid[0] < R_end && !move_u_to_S_with_prune(0, S_end, R_end, level))
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		choose_zero = false;
		ui S2_n = 0;
		for (ui i = 0; i < S_end; i++)
			if (R_end - degree[SR[i]] > s)
				S2[S2_n++] = SR[i];

#ifdef _SUPPORT_BOUND_
		if (S2_n >= 2)
		{
			collect_removable_vertices_based_on_total_edges(S2_n, S_end, R_end, level);
			if (!remove_vertices_and_prune(S_end, R_end, level))
			{

				restore_SR(S_end, R_end, old_S_end, old_R_end, level);
				return;
			}
		}
#endif

		if (R_end > best_solution_size && is_splex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif

#ifdef _VP_BOUND_
		if (S2_n >= 2)
		{
			bool ok = collect_removable_vertices_based_on_vp_bound(S_end, R_end, level);
			if (!ok)
			{

				restore_SR(S_end, R_end, old_S_end, old_R_end, level);
				return;
			}
		}
#endif

		ui u = n; // u is the branching vertex
		ui min_deg = n;
		assert(!choose_zero);
		for (ui i = 0; i < R_end; i++)
		{
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}
		}
		assert(u != n && SR_rid[u] < R_end && min_deg != n);

#ifdef _SBUNDLE_BRANCHING_
		if (min_deg >= R_end - s)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			ui pre_S_end = S_end, pre_R_end = R_end;
			SymBD_L(S_end, R_end, level + 1, choose_zero);
			assert(S_end == pre_S_end && R_end == pre_R_end);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}
#else
		if (min_deg >= R_end - s)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			ui pre_S_end = S_end, pre_R_end = R_end;
			sbundle_binary_BB(S_end, R_end, level + 1);
			assert(S_end == pre_S_end && R_end == pre_R_end);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}
#endif

		if (SR_rid[u] < S_end)
		{
			u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end, u);
		}
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
		assert(degree[u] + s > best_solution_size && degree[u] + s > S_end);

		swap_pos(S_end, SR_rid[u]);
		bool can_add_u = vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);

		if (can_add_u)
		{
			// the first branch includes u into s
			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			assert(Qv.empty());

#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < S_end; i++)
				assert(degree_in_S[SR[i]] + s >= S_end);
			for (ui i = S_end; i < R_end; i++)
				assert(degree_in_S[SR[i]] + s > S_end);
			for (ui i = 0; i < S_end; i++)
				if (degree_in_S[SR[i]] + s == S_end)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = S_end; j < R_end; j++)
						assert(t_matrix[SR[j]]);
				}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
			if (move_u_to_S_with_prune(u, S_end, R_end, level))
			{
#ifndef NDEBUG
				for (ui i = 0; i < R_end; i++)
				{
					ui d1 = 0, d2 = 0;
					for (ui j = 0; j < S_end; j++)
						if (matrix[SR[i] * n + SR[j]])
							++d1;
					d2 = d1;
					for (ui j = S_end; j < R_end; j++)
						if (matrix[SR[i] * n + SR[j]])
							++d2;
					assert(d1 == degree_in_S[SR[i]]);
					assert(d2 == degree[SR[i]]);
				}
				for (ui i = 0; i < S_end; i++)
					assert(degree_in_S[SR[i]] + s >= S_end);
				for (ui i = S_end; i < R_end; i++)
					assert(degree_in_S[SR[i]] + s > S_end);
				for (ui i = 0; i < S_end; i++)
					if (degree_in_S[SR[i]] + s == S_end)
					{
						char *t_matrix = matrix + SR[i] * n;
						for (ui j = S_end; j < R_end; j++)
							assert(t_matrix[SR[j]]);
					}
				for (ui i = 0; i < R_end; i++)
					assert(level_id[SR[i]] > level);
#endif
				SymBD_H(S_end, R_end, level + 1, false);
			}

			if (best_solution_size >= _UB_)
				return;
			assert(S_end == t_old_S_end + 1 && SR[S_end - 1] == u);
			restore_SR(S_end, R_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif

			// the second branch exclude u from s
			assert(Qv.empty());
			bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, R_end, level);

			if (succeed)
			{
				SymBD_H(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < S_end; i++)
				assert(degree_in_S[SR[i]] + s >= S_end);
			for (ui i = S_end; i < R_end; i++)
				assert(degree_in_S[SR[i]] + s > S_end);
			for (ui i = 0; i < S_end; i++)
				if (degree_in_S[SR[i]] + s == S_end)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = S_end; j < R_end; j++)
						assert(t_matrix[SR[j]]);
				}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
		}
		else
		{
			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			assert(Qv.empty());
			bool succeed = remove_u_from_C_with_prune(u, S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, R_end, level);

			if (succeed)
			{
				SymBD_H(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
	}

	void sbundle_binary_BB(ui S_end, ui R_end, ui level)
	{
		if (S_end > best_solution_size && vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end, matrix, n, s, true))
			store_solution(S_end);
		if (R_end > best_solution_size && is_splex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1)
			return;

		branch_count++;

		ui u = n; // u is the branching vertex
		ui min_deg = n;
		for (ui i = S_end; i < R_end; i++)
		{
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}
		}

		assert(u != n && min_deg != n);

		// The first branch
		assert(SR_rid[u] < R_end && SR_rid[u] >= S_end);
		delfrC(u, R_end, level);
		sbundle_binary_BB(S_end, R_end, level + 1);
		addtoC(u, R_end, level);

		// The second branch
		assert(level_id[u] > level);
		if (SR_rid[u] != S_end)
		{
			swap_pos(S_end, SR_rid[u]);
		}
		bool can_add = vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);
		if (can_add)
		{
			addtoS(u, S_end);
			sbundle_binary_BB(S_end, R_end, level + 1);
			delfrS(u, S_end);
		}
	}

	void SymBD_L(ui S_end, ui R_end, ui level, bool choose_zero)
	{
		if (S_end > best_solution_size)
			store_solution(S_end);
		if (R_end > best_solution_size && is_splex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;
		branch_count++;

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			char *t_matrix = matrix + u * n;
			for (ui j = i + 1; j < R_end; j++)
			{
				ui v = SR[j];
				if (!t_matrix[v])
				{
					ui cnt = 0;
					std::vector<char> t_s_connected;
					if (u > v)
					{
						t_s_connected = all_s_connected[std::make_pair(v, u)];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(u, v)];
					}
					assert(!t_s_connected.empty());
					for (ui j = 0; j < R_end; j++)
					{
						int w = SR[j];
						if (t_s_connected[w])
							cnt++;
					}
					assert(local_vc[n * v + u] == cnt && local_vc[n * u + v] == cnt);
				}
			}
		}
#endif
		ui *local_vc_backup = nullptr;
		std::map<std::pair<ui, ui>, std::vector<char>> all_s_connected_backup;

		bool need_update_local_vc = false;
		if (isSbundle_preliminary(R_end, true))
		{
			if (vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			{
				assert(R_end > best_solution_size);
				store_solution(R_end);
				return;
			}
			else
			{
				if (R_end == best_solution_size + 1)
					return;
				need_update_local_vc = true;
			}
		}

		if (need_update_local_vc)
		{
			update_local_vc(local_vc, local_vc_backup, all_s_connected, all_s_connected_backup, R_end);
		}

		// if (S_end >= 2)
		// {
		// 	ui vp_bound = get_vertex_pair_ub(S_end, R_end);
		// 	if (vp_bound <= best_solution_size)
		// 	{
		// 		if (need_update_local_vc)
		// 		{
		// 			restore_local_vc(local_vc, local_vc_backup, all_s_connected, all_s_connected_backup);
		// 		return;
		// 	}
		// }

		if (choose_zero && SR_rid[0] < R_end)
		{
			assert(S_end == 0 && SR_rid[0] < R_end);
			addtoS(0, S_end);
		}

		assert(R_end > best_solution_size);
		ui u = n;
		ui u2 = n;
		ui min_vc = R_end - s;
		int branching_case = 4;
		ui min_vc_in_S = n;

		find_branching_vertex_pair(u, u2, S_end, S_end, branching_case, 1, min_vc, min_vc_in_S);
		if (min_vc_in_S + s <= best_solution_size)
		{
			if (need_update_local_vc)
				restore_local_vc(local_vc, local_vc_backup, all_s_connected, all_s_connected_backup);
			return;
		}
		if (branching_case == 4)
		{
			find_branching_vertex_pair(u, u2, S_end, R_end, branching_case, 2, min_vc, min_vc_in_S);
		}
		if (branching_case == 4)
		{
			find_branching_vertex_pair(u, u2, S_end, R_end, branching_case, 3, min_vc, min_vc_in_S);
		}
		assert(branching_case != 4 && (u >= 0 && u < R_end) && (u2 >= 0 && u2 < R_end));

		// Specific branching cases
		if (branching_case == 1)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				assert(d1 == degree[SR[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
			LVC_based_branching(u, u2, S_end, R_end, level);
		}
		else if (branching_case == 2)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				assert(d1 == degree[SR[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
			simple_binary_branching_1(u, u2, S_end, R_end, level, min_vc);
		}
		else if (branching_case == 3)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				assert(d1 == degree[SR[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
			simple_binary_branching_2(u, u2, S_end, R_end, level);
		}

		if (need_update_local_vc)
		{
			restore_local_vc(local_vc, local_vc_backup, all_s_connected, all_s_connected_backup);
		}
	}

	// Find an s-unsatisfied vertex pair
	void find_branching_vertex_pair(ui &branch_vertex_1, ui &branch_vertex_2, ui pointer_1, ui pointer_2, int &branching_case, int branching_case_value, ui &min_vc, ui &min_vc_in_S)
	{
		if (branching_case_value == 1)
		{
			for (ui i = 0; i < pointer_1; i++)
			{
				ui *t_local_vc = local_vc + SR[i] * n;
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = i + 1; j < pointer_2; j++)
				{
					if (!t_matrix[SR[j]])
					{
						ui temp_vc = t_local_vc[SR[j]];
						if (temp_vc < min_vc)
						{
							min_vc = temp_vc;
							branch_vertex_1 = i;
							branch_vertex_2 = j;
							branching_case = branching_case_value;
						}
						if (temp_vc < min_vc_in_S)
						{
							min_vc_in_S = temp_vc;
						}
					}
				}
			}
		}
		else if (branching_case_value == 2)
		{
			for (ui i = 0; i < pointer_1; i++)
			{
				char *t_matrix = matrix + SR[i] * n;
				ui *t_local_vc = local_vc + SR[i] * n;
				for (ui j = pointer_1; j < pointer_2; j++)
				{
					if (!t_matrix[SR[j]])
					{
						ui temp_vc = t_local_vc[SR[j]];
						if (temp_vc < min_vc)
						{
							min_vc = temp_vc;
							branch_vertex_1 = i;
							branch_vertex_2 = j;
							branching_case = branching_case_value;
						}
					}
				}
			}
		}
		else if (branching_case_value == 3)
		{
			for (ui i = pointer_1; i < pointer_2; i++)
			{
				ui *t_local_vc = local_vc + SR[i] * n;
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = i + 1; j < pointer_2; j++)
				{
					if (!t_matrix[SR[j]])
					{
						ui temp_vc = t_local_vc[SR[j]];
						if (temp_vc < min_vc)
						{
							min_vc = temp_vc;
							branch_vertex_1 = i;
							branch_vertex_2 = j;
							branching_case = branching_case_value;
						}
					}
				}
			}
		}
	}

	void restore_deleted_vertices(ui &R_end, ui level, ui &count_restore)
	{
		while (!temp_deleted.empty())
		{
			ui last_deleted = temp_deleted.back();
			if (level_id[last_deleted] == level)
			{
				addtoC(last_deleted, R_end, level);
				temp_deleted.pop_back();
				count_restore++;
			}
			else
			{
				break;
			}
		}
	}

	void update_local_vc(ui *local_vc, ui *&local_vc_backup, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected_backup, ui R_end)
	{
		local_vc_backup = new ui[sizeof(ui) * ((long long)n) * n];
		std::copy(local_vc, local_vc + ((long long)n) * n, local_vc_backup);
		memset(local_vc, 0, sizeof(ui) * ((long long)n) * n);
		all_s_connected_backup = all_s_connected;
		all_s_connected.clear();
		vc->get_local_vertexConnectivity(R_end, SR, matrix, n, local_vc, all_s_connected, s);
		assert(local_vc != nullptr && local_vc_backup != nullptr);
	}

	void restore_local_vc(ui *local_vc, ui *&local_vc_backup, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected, std::map<std::pair<ui, ui>, std::vector<char>> all_s_connected_backup)
	{
		assert(local_vc != nullptr && local_vc_backup != nullptr);
		std::copy(local_vc_backup, local_vc_backup + ((long long)n) * n, local_vc);
		if (local_vc_backup != nullptr)
		{
			delete[] local_vc_backup;
			local_vc_backup = nullptr;
		}
		all_s_connected.clear();
		all_s_connected = all_s_connected_backup;
	}

	void LVC_based_branching(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u < S_end && u2 < S_end);

		find_s_connected(s_connected, u, u2);
		assert(!s_connected.empty());

		ui notadj_S = 0;
		std::vector<ui> notadj_C;
		get_notadj_S_and_notadj_C(notadj_S, notadj_C, S_end, R_end, level);
		assert(!notadj_C.empty() && notadj_S <= s);

		assert(s >= notadj_S);
		ui canselect = s - notadj_S;
		ui pos = n;
		assert(canselect >= 0 && canselect < notadj_C.size());

		bool all = true;
		generate_canselect_branches(all, canselect, S_end, R_end, level, pos, notadj_C);

		generate_the_last_branch(all, canselect, S_end, R_end, level, notadj_C);

		if (pos != n)
		{
			for (ui i = pos; i >= 0 && i <= notadj_C.size(); i--)
			{
				delfrS(notadj_C[i], S_end);
			}
		}
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
	}

	void simple_binary_branching_1(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level, ui min_vc)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u < S_end && u2 >= S_end);
		ui record_index = SR[u2];

		// The first branch with prune
		if (min_vc + s > best_solution_size && canadd(record_index, S_end, true))
		{
			assert(level_id[record_index] > level);
			addtoS(record_index, S_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C(S_end, R_end, level, count_deleted);
			SymBD_L(S_end, R_end, level + 1, false);
			restore_deleted_vertices(R_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS(record_index, S_end);
		}

		// The second branch
		assert(SR_rid[record_index] < R_end && SR_rid[record_index] >= S_end);
		delfrC(record_index, R_end, level);
		SymBD_L(S_end, R_end, level + 1, false);
		addtoC(record_index, R_end, level);
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
	}

	void simple_binary_branching_2(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u >= S_end && u2 >= S_end);
		ui record_index = SR[u2];
		assert(level_id[record_index] > level);
		if (canadd(record_index, S_end, true))
		{
			addtoS(record_index, S_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C(S_end, R_end, level, count_deleted);
			SymBD_L(S_end, R_end, level + 1, false);
			restore_deleted_vertices(R_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS(record_index, S_end);
		}

		// The second branch
		assert(SR_rid[record_index] < R_end && SR_rid[record_index] >= S_end);
		delfrC(record_index, R_end, level);
		SymBD_L(S_end, R_end, level + 1, false);
		addtoC(record_index, R_end, level);
		assert(pre_Rend == R_end);
	}

	void generate_canselect_branches(bool &all, ui canselect, ui &S_end, ui &R_end, ui level, ui &pos, std::vector<ui> notadj_C)
	{
		for (ui i = 0; i < canselect; ++i)
		{
			assert(SR_rid[notadj_C[i]] < R_end && SR_rid[notadj_C[i]] >= S_end);
			delfrC(notadj_C[i], R_end, level);
			if (i && !canadd(notadj_C[i - 1], S_end, true))
			{
				addtoC(notadj_C[i], R_end, level);
				all = false;
				break;
			}
			ui count_deleted = 0;
			if (i)
			{
				assert(canadd(notadj_C[i - 1], S_end, true));
				addtoS(notadj_C[i - 1], S_end);
				delete_vertices_from_C(S_end, R_end, level, count_deleted);
				pos = i - 1;
			}

			SymBD_L(S_end, R_end, level + 1, false);

			ui count_restore = 0;
			if (i)
			{
				restore_deleted_vertices(R_end, level, count_restore);
			}
			assert(count_restore == count_deleted);
			addtoC(notadj_C[i], R_end, level);
		}
	}

	void generate_the_last_branch(bool all, ui canselect, ui &S_end, ui &R_end, ui level, std::vector<ui> notadj_C)
	{
		if (all && (canselect == 0 || canadd(notadj_C[canselect - 1], S_end, true)))
		{
			bool t_canadd = false;
			if (canselect)
			{
				t_canadd = canadd(notadj_C[canselect - 1], S_end, true);
			}
			ui pre_dd = 0;
			for (ui i = canselect; i < (ui)notadj_C.size(); ++i)
			{
				assert(level_id[notadj_C[i]] > level);
				if (level_id[notadj_C[i]] > level)
				{
					assert(SR_rid[notadj_C[i]] < R_end && SR_rid[notadj_C[i]] >= S_end);
					delfrC(notadj_C[i], R_end, level);
					pre_dd++;
				}
			}

			ui count_deleted = 0;
			if (canselect)
			{
				assert(t_canadd);
				addtoS(notadj_C[canselect - 1], S_end);
				delete_vertices_from_C(S_end, R_end, level, count_deleted);
			}

			SymBD_L(S_end, R_end, level + 1, false);
			ui count_restore = 0;
			if (canselect)
			{
				restore_deleted_vertices(R_end, level, count_restore);
				delfrS(notadj_C[canselect - 1], S_end);
			}
			assert(count_restore == count_deleted);

			ui after_dd = 0;
			for (ui i = notadj_C.size() - 1; i >= canselect && i <= notadj_C.size(); i--)
			{
				if (level_id[notadj_C[i]] == level)
				{
					addtoC(notadj_C[i], R_end, level);
					after_dd++;
				}
			}
			assert(pre_dd == after_dd);
		}
	}

	void get_notadj_S_and_notadj_C(ui &notadj_S, std::vector<ui> &notadj_C, ui S_end, ui R_end, ui level)
	{
		assert(notadj_C.empty());
		for (ui i = 0; i < R_end; i++)
		{
			if (s_connected[SR[i]] == 0)
			{
				if (i < S_end)
				{
					notadj_S++;
				}
				else
				{
					assert(level_id[SR[i]] > level);
					notadj_C.push_back(SR[i]);
				}
			}
		}

#ifdef _BETTER_SEQEUNCE_
		ui t_n = notadj_C.size();
		std::vector<ui> scores(t_n);

		for (ui i = 0; i < t_n; ++i)
		{
			ui u = notadj_C[i];
			scores[i] = degree[u];
		}
		std::vector<ui> indices(t_n);
		for (ui i = 0; i < t_n; ++i)
		{
			indices[i] = i;
		}
		std::sort(indices.begin(), indices.end(), [&scores](ui i, ui j)
				  { return scores[i] < scores[j]; });

		std::vector<ui> sorted_arr1(t_n);
		for (ui i = 0; i < t_n; ++i)
		{
			sorted_arr1[i] = notadj_C[indices[i]];
		}

		notadj_C = sorted_arr1;
#endif

		assert(notadj_S <= s);
		assert(notadj_C.size() > s - notadj_S);
	}

	void delete_vertices_from_C(ui S_end, ui &R_end, ui level, ui &count_deleted)
	{
		for (ui i = S_end; i < R_end; i++)
		{
			ui temp_record = SR[i];
			if (level_id[temp_record] > level)
			{
				bool can_delete = false;
				if (degree[temp_record] < best_solution_size - s)
				{
					can_delete = true;
				}
				else if (!canadd(temp_record, S_end, false))
				{
					can_delete = true;
				}
				if (can_delete)
				{
					temp_deleted.push_back(temp_record);
					assert(level_id[temp_record] > level);
					assert(SR_rid[temp_record] < R_end && SR_rid[temp_record] >= S_end);
					delfrC(temp_record, R_end, level);
					count_deleted++;
				}
			}
		}
	}

	void find_s_connected(std::vector<char> &s_connected, ui vertex_1, ui vertex_2)
	{
		ui u = vertex_1;
		ui u2 = vertex_2;
		s_connected.clear();
		if (SR[u] > SR[u2])
		{
			s_connected = all_s_connected[std::make_pair(SR[u2], SR[u])];
		}
		else
		{
			s_connected = all_s_connected[std::make_pair(SR[u], SR[u2])];
		}
		assert(!s_connected.empty());
	}

	bool isSbundle_preliminary(ui end, bool is_R_end)
	{
		if (end <= s)
			return true;
		if (is_R_end)
		{
			for (ui i = 0; i < end; i++)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = i + 1; j < end; j++)
				{
					if (!t_matrix[SR[j]] && local_vc[SR[i] * n + SR[j]] < end - s)
						return false;
				}
			}
		}
		return true;
	}

	// Our proposed new upper bound
	ui get_vertex_pair_ub(ui S_end, ui R_end)
	{
		assert(S_end >= 2);
		std::vector<std::pair<ui, ui>> vp_in_S;

		for (ui i = 0; i < S_end; i++)
		{
			char *t_matrix = matrix + SR[i] * n;
			for (ui j = i + 1; j < S_end; j++)
			{
				if (!t_matrix[SR[j]])
				{
					if (SR[i] > SR[j])
					{
						vp_in_S.push_back(std::make_pair(SR[j], SR[i]));
					}
					else
					{
						vp_in_S.push_back(std::make_pair(SR[i], SR[j]));
					}
				}
			}
		}

		if (vp_in_S.size() == 0)
			return R_end;

		ui ret_ub = S_end;
		assert(!vp_in_S.empty());
		ui t_size = vp_in_S.size();
		std::vector<std::vector<ui>> lost_in_S(t_size, std::vector<ui>(2));

		for (ui i = 0; i < t_size; i++)
		{
			lost_in_S[i][0] = i;
		}
		std::vector<ui> partitions(t_size, 0);
		ui connect_all_vp_in_S = R_end - S_end;
		std::vector<ui> is_all_sconnect_with_S(R_end, 1);

		for (size_t i = 0; i < vp_in_S.size(); i++)
		{
			std::pair<ui, ui> cur_vp = vp_in_S[i];
			std::vector<char> cur_s_con = all_s_connected[cur_vp];
			assert(!cur_s_con.empty());
			for (ui j = 0; j < R_end; j++)
			{
				if (cur_s_con[SR[j]] == 0)
				{
					if (j < S_end)
					{
						if (lost_in_S[i][1] == s)
						{
							return 0;
						}
						lost_in_S[i][1]++;
					}
					else
					{
						if (is_all_sconnect_with_S[j] == 1)
						{
							is_all_sconnect_with_S[j] = 0;
							assert(connect_all_vp_in_S > 0);
							connect_all_vp_in_S--;
						}
					}
				}
			}
			assert(s >= lost_in_S[i][1]);
		}

		ret_ub = ret_ub + connect_all_vp_in_S;
		if (ret_ub > best_solution_size)
		{
			return R_end;
		}

		std::sort(lost_in_S.begin(), lost_in_S.end(), [](const std::vector<ui> &a, const std::vector<ui> &b)
				  { return a[1] > b[1]; });

		for (ui i = S_end; i < R_end; i++)
		{
			if (is_all_sconnect_with_S[i] == 0)
			{
				ui t_v = SR[i];
				for (ui j = 0; j < t_size; j++)
				{
					ui t_index = lost_in_S[j][0];
					std::vector<char> cur_s_con = all_s_connected[vp_in_S[t_index]];
					assert(!cur_s_con.empty());
					if (cur_s_con[t_v] == 0)
					{
						partitions[t_index]++;
						break;
					}
				}
			}
		}

#ifndef NDEBUG
		ui result = 0;
		for (ui i = 0; i < t_size; i++)
		{
			result += partitions[i];
		}
		result += connect_all_vp_in_S;
		assert(result == R_end - S_end);
#endif

		for (ui i = 0; i < t_size; i++)
		{
			ui t_index = lost_in_S[i][0];
			assert(lost_in_S[i][1] <= s);
			if (s - lost_in_S[i][1] > partitions[t_index])
			{
				ret_ub += partitions[t_index];
			}
			else
			{
				ret_ub += s - lost_in_S[i][1];
			}
			// early termination
			if (ret_ub > best_solution_size)
				return R_end;
		}

		assert(ret_ub <= R_end);
		return ret_ub;
	}

	bool canadd(ui u, ui S_end, bool precise)
	{
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < S_end; i++)
		{
			if (!t_matrix[SR[i]])
			{
				if (local_vc[u * n + SR[i]] < best_solution_size - s)
				{
					return false;
				}
			}
		}

		if (precise)
		{
			if (SR_rid[u] != S_end)
			{
				swap_pos(S_end, SR_rid[u]);
			}
			return vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);
		}

		return true;
	}

	bool removeByVpUB(ui u, ui S_end, ui R_end)
	{
		addtoS(u, S_end);
		ui vp_bound = get_vertex_pair_ub(S_end, R_end);
		if (vp_bound > best_solution_size)
		{
			delfrS(u, S_end);
			return false;
		}
		delfrS(u, S_end);
		return true;
	}

	void delfrC(ui u, ui &R_end, ui level)
	{
		assert(level_id[u] == n);
		level_id[u] = level;
		--R_end;
		swap_pos(R_end, SR_rid[u]);
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
			{
				assert(degree[SR[i]] > 0);
				--degree[SR[i]];
			}

		for (ui i = 0; i < R_end; i++)
		{
			char *t_matrix = matrix + SR[i] * n;
			for (ui j = i + 1; j < R_end; j++)
			{
				if (!t_matrix[SR[j]])
				{
					std::vector<char> t_s_connected;
					if (SR[i] > SR[j])
					{
						t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
					}
					assert(!t_s_connected.empty());
					if (t_s_connected[u] == 1)
					{
						assert((local_vc[SR[i] * n + SR[j]] > 0) && (local_vc[SR[j] * n + SR[i]] > 0));
						local_vc[SR[i] * n + SR[j]]--;
						local_vc[SR[j] * n + SR[i]]--;
					}
				}
			}
		}
	}

	void addtoC(ui u, ui &R_end, ui level)
	{
		assert(level_id[u] != n && level_id[u] == level);
		level_id[u] = n;
		{
			for (ui i = 0; i < R_end; i++)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = i + 1; j < R_end; j++)
				{
					if (!t_matrix[SR[j]])
					{
						std::vector<char> t_s_connected;
						if (SR[i] > SR[j])
						{
							t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
						}
						else
						{
							t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
						}
						assert(!t_s_connected.empty());
						if (t_s_connected[u] == 1)
						{
							assert((local_vc[SR[i] * n + SR[j]] < n) && (local_vc[SR[j] * n + SR[i]] < n));
							local_vc[SR[i] * n + SR[j]]++;
							local_vc[SR[j] * n + SR[i]]++;
						}
					}
				}
			}

			for (ui i = 0; i < R_end; i++)
			{
				char *t_matrix = matrix + SR[i] * n;
				if (!t_matrix[u])
				{
					ui v = SR[i];
					ui cnt = 0;
					std::vector<char> t_s_connected;
					if (u > v)
					{
						t_s_connected = all_s_connected[std::make_pair(v, u)];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(u, v)];
					}
					assert(!t_s_connected.empty());
					for (ui j = 0; j < R_end; j++)
					{
						int w = SR[j];
						if (t_s_connected[w])
							cnt++;
					}
					if (t_s_connected[u])
						cnt++;
					local_vc[n * v + u] = local_vc[n * u + v] = cnt;
				}
			}
		}

		assert(SR_rid[u] == R_end);
		SR_rid[u] = R_end;
		SR[R_end] = u;
		++R_end;

		char *t_matrix = matrix + u * n;
		degree[u] = 0;
		for (ui i = 0; i < R_end; i++)
		{
			if (t_matrix[SR[i]])
			{
				++degree[SR[i]];
				degree[u]++;
			}
		}
	}

	void addtoS(ui u, ui &S_end)
	{
		swap_pos(S_end, SR_rid[u]);
		++S_end;
	}

	void delfrS(ui u, ui &S_end)
	{
		--S_end;
		assert(S_end == SR_rid[u]);
	}

	void delete_vertices_based_on_lower_bound(ui S_end, ui R_end, std::vector<ui> &results)
	{
		for (ui i = S_end; i < R_end; i++)
		{
			ui u = SR[i];
			char *t_matrix = matrix + u * n;
			for (ui j = 0; j < S_end; j++)
			{
				ui v = SR[j];
				assert(u != v);
				if (!t_matrix[v] && local_vc[u * n + v] <= best_solution_size - s) // UB2
				{
					results.push_back(u);
					break;
				}
			}
		}
	}

	bool collect_removable_vertices_based_on_vp_bound(ui S_end, ui &R_end, ui level)
	{
		ui pre_S_end = S_end;
		ui pre_R_end = R_end;
		for (ui i = pre_S_end; i < pre_R_end; i++)
		{
			ui u = SR[i];
			if (level_id[u] > level)
			{
				addtoS(u, S_end);
				assert(S_end == pre_S_end + 1);
				ui t_ub = get_vertex_pair_ub(S_end, R_end); // UB3
				if (t_ub <= best_solution_size)
				{
					delfrS(u, S_end);
					level_id[u] = level;
					Qv.push(u);
					if (!remove_vertices_and_prune(S_end, R_end, level))
					{
						return false;
					}
					continue;
				}
				assert(SR_rid[u] + 1 == S_end);
				delfrS(u, S_end);
			}
			assert(S_end == pre_S_end);
		}
		assert(S_end == pre_S_end);
		return true;
	}

	void collect_removable_vertices_based_on_total_edges(ui S2_n, ui S_end, ui R_end, ui level)
	{
		vp.resize(R_end - S_end);
		ui max_nn = 0;
		for (ui i = S_end; i < R_end; i++)
		{
			ui nn = 0;
			if (S2_n != S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = 0; j < S2_n; j++)
					if (!t_matrix[S2[j]])
						++nn;
			}
			else
				nn = S_end - degree_in_S[SR[i]];
			if (nn > max_nn)
				max_nn = nn;
			vp[i - S_end].first = SR[i];
			vp[i - S_end].second = nn;
		}
		ui *cnt = neighbors;
		for (ui i = 0; i <= max_nn; i++)
			cnt[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
			++cnt[vp[i].second];
		for (ui i = 0; i < max_nn; i++)
			cnt[i + 1] += cnt[i];
		for (ui i = max_nn; i > 0; i--)
			cnt[i] = cnt[i - 1];
		cnt[0] = 0;
		ui *ids = nonneighbors;
		for (ui i = 0; i < vp.size(); i++)
			ids[cnt[vp[i].second]++] = i;

		ui total_support = 0;
		for (ui i = 0; i < S2_n; i++)
			total_support += s - S_end + degree_in_S[S2[i]];

		ui new_n = 0;
		while (!Qv.empty())
			Qv.pop();
		for (ui i = 0; i < vp.size(); i++)
		{
			ui idx = ids[i], v = vp[ids[i]].first;
			ui t_support = total_support - vp[idx].second;
			char *t_matrix = matrix + v * n;
			ui j = 0, v_support = s - 1 - S_end + degree_in_S[v], ub = S_end + 1;
			while (true)
			{
				if (j == new_n)
					j = i + 1;
				if (j >= vp.size() || ub > best_solution_size || ub + vp.size() - j <= best_solution_size)
					break;
				ui u = vp[ids[j]].first, nn = vp[ids[j]].second;
				if (t_support < nn)
					break;
				if (t_matrix[u])
				{
					t_support -= nn;
					++ub;
				}
				else if (v_support > 0)
				{
					--v_support;
					t_support -= nn;
					++ub;
				}
				++j;
			}
			if (ub <= best_solution_size)
			{
				level_id[v] = level;
				Qv.push(v);
			}
			else
				ids[new_n++] = ids[i];
		}
	}

	void get_neighbors_and_nonneighbors(ui u, ui R_end, ui &neighbors_n, ui &nonneighbors_n)
	{
		neighbors_n = 0;
		nonneighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (SR[i] != u)
			{
				if (t_matrix[SR[i]])
					neighbors[neighbors_n++] = SR[i];
				else
					nonneighbors[nonneighbors_n++] = SR[i];
			}
	}

	bool move_u_to_S_with_prune(ui u, ui &S_end, ui &R_end, ui level)
	{
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end && SR[SR_rid[u]] == u);
		assert(degree_in_S[u] + s > S_end);
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
		if (SR_rid[u] != S_end)
			swap_pos(S_end, SR_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_nonneighbors(u, R_end, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
#endif

		while (!Qv.empty())
			Qv.pop();
		// reduction rules based on the fact that each vertex can have at most s-1 nonneighbors
		for (ui i = 0; i < nonneighbors_n; i++)
		{
			ui v = nonneighbors[i];
			if (SR_rid[v] >= S_end)
			{
				if (level_id[v] == level)
					continue;
				if (S_end - degree_in_S[v] >= s || S_end - degree_in_S[u] == s)
				{
					level_id[v] = level;
					Qv.push(v);
				}
			}
			else if (S_end - degree_in_S[v] == s)
			{
				char *tt_matrix = matrix + v * n;
				for (ui j = S_end; j < R_end; j++)
					if (level_id[SR[j]] > level && !tt_matrix[SR[j]])
					{
						level_id[SR[j]] = level;
						Qv.push(SR[j]);
					}
			}
		}

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(level_id[SR[j]] == level || t_matrix[SR[j]]);
			}
#endif
		return remove_vertices_and_prune(S_end, R_end, level);
	}

	bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level)
	{
		while (true)
		{
			while (!Qv.empty())
			{
				ui u = Qv.front();
				Qv.pop(); // remove u
				assert(SR[SR_rid[u]] == u && SR_rid[u] >= S_end && SR_rid[u] < R_end);
				--R_end;
				swap_pos(SR_rid[u], R_end);

				// Dynamically maintain local vertex connectivity
				for (ui i = 0; i < R_end; i++)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = i + 1; j < R_end; j++)
					{
						if (!t_matrix[SR[j]])
						{
							std::vector<char> t_s_connected;
							if (SR[i] > SR[j])
							{
								t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
							}
							else
							{
								t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
							}
							assert(!t_s_connected.empty());
							if (t_s_connected[u] == 1)
							{
								assert((local_vc[SR[i] * n + SR[j]] > 0) && (local_vc[SR[j] * n + SR[i]] > 0));
								local_vc[SR[i] * n + SR[j]]--;
								local_vc[SR[j] * n + SR[i]]--;
							}
						}
					}
				}

				bool terminate = false;
				ui neighbors_n = 0;
				char *t_matrix = matrix + u * n;
				for (ui i = 0; i < R_end; i++)
					if (t_matrix[SR[i]])
					{
						ui w = SR[i];
						neighbors[neighbors_n++] = w;
						--degree[w];
						if (degree[w] + s <= best_solution_size)
						{
							if (i < S_end)
								terminate = true;
							else if (level_id[w] > level)
							{
								level_id[w] = level;
								Qv.push(w);
							}
						}
					}

				if (terminate)
				{
					for (ui i = 0; i < neighbors_n; i++)
						++degree[neighbors[i]];
					level_id[u] = n;
					{
						for (ui i = 0; i < R_end; i++)
						{
							char *t_matrix = matrix + SR[i] * n;
							for (ui j = i + 1; j < R_end; j++)
							{
								if (!t_matrix[SR[j]])
								{
									std::vector<char> t_s_connected;
									if (SR[i] > SR[j])
									{
										t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
									}
									else
									{
										t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
									}
									assert(!t_s_connected.empty());
									if (t_s_connected[u] == 1)
									{
										assert((local_vc[SR[i] * n + SR[j]] < n) && (local_vc[SR[j] * n + SR[i]] < n));
										local_vc[SR[i] * n + SR[j]]++;
										local_vc[SR[j] * n + SR[i]]++;
									}
								}
							}
						}

						for (ui i = 0; i < R_end; i++)
						{
							char *t_matrix = matrix + SR[i] * n;
							if (!t_matrix[u])
							{
								ui v = SR[i];
								ui cnt = 0;
								std::vector<char> t_s_connected;
								if (u > v)
								{
									t_s_connected = all_s_connected[std::make_pair(v, u)];
								}
								else
								{
									t_s_connected = all_s_connected[std::make_pair(u, v)];
								}
								assert(!t_s_connected.empty());
								for (ui j = 0; j < R_end; j++)
								{
									ui w = SR[j];
									if (t_s_connected[w])
										cnt++;
								}
								if (t_s_connected[u])
									cnt++;
								local_vc[n * v + u] = local_vc[n * u + v] = cnt;
							}
						}
					}
					++R_end;
					return false;
				}
			}

#ifdef _LOCAL_VC_DELETE_
			if (S_end >= 1)
			{
				std::vector<ui> deleted_vertices;
				delete_vertices_based_on_lower_bound(S_end, R_end, deleted_vertices);
				if (!deleted_vertices.empty())
				{
					for (ui i = 0; i < deleted_vertices.size(); i++)
					{
						ui t_v = deleted_vertices[i];
						assert(level_id[t_v] > level && SR_rid[t_v] >= S_end && SR_rid[t_v] < R_end);
						level_id[t_v] = level;
						Qv.push(t_v);
					}
				}
			}
#endif
			if (Qv.empty())
			{
				break;
			}
		}

		assert(Qv.empty());
		return true;
	}

	void restore_SR(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level)
	{
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
#endif
		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			if (level_id[SR[i]] <= level)
				printf("level_id[%u] = %u, level = %u, n = %u\n", SR[i], level_id[SR[i]], level, n);
			assert(level_id[SR[i]] > level);
		}
#endif
		for (; R_end < old_R_end; R_end++)
		{ // insert u back into R
			ui u = SR[R_end];
			assert(level_id[u] == level && SR_rid[u] == R_end);
			level_id[u] = n;
			{
				for (ui i = 0; i < R_end; i++)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = i + 1; j < R_end; j++)
					{
						if (!t_matrix[SR[j]])
						{
							std::vector<char> t_s_connected;
							if (SR[i] > SR[j])
							{
								t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
							}
							else
							{
								t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
							}
							assert(!t_s_connected.empty());
							if (t_s_connected[u] == 1)
							{
								assert((local_vc[SR[i] * n + SR[j]] < n) && (local_vc[SR[j] * n + SR[i]] < n));
								local_vc[SR[i] * n + SR[j]]++;
								local_vc[SR[j] * n + SR[i]]++;
							}
						}
					}
				}

				for (ui i = 0; i < R_end; i++)
				{
					char *t_matrix = matrix + SR[i] * n;
					if (!t_matrix[u])
					{
						int v = SR[i];
						int cnt = 0;
						std::vector<char> t_s_connected;
						if (u > v)
						{
							t_s_connected = all_s_connected[std::make_pair(v, u)];
						}
						else
						{
							t_s_connected = all_s_connected[std::make_pair(u, v)];
						}
						assert(!t_s_connected.empty());
						for (ui j = 0; j < R_end; j++)
						{
							int w = SR[j];
							if (t_s_connected[w])
								cnt++;
						}
						if (t_s_connected[u])
							cnt++;
						local_vc[n * v + u] = local_vc[n * u + v] = cnt;
					}
				}
			}

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					++degree[w];
					++degree[u];
					if (i < S_end)
						++degree_in_S[u];
				}
		}

		for (; S_end > old_S_end; S_end--) // move u from s to R
		{
			ui u = SR[S_end - 1];
			assert(SR_rid[u] == S_end - 1);

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					--degree_in_S[w];
				}
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(old_S_end == S_end || degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(old_S_end == S_end || t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
	}

	bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = SR[S_end - 1];
		--S_end;
		--R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		// Dynamically maintain local vertex connectivity
		for (ui i = 0; i < R_end; i++)
		{
			char *t_matrix = matrix + SR[i] * n;
			for (ui j = i + 1; j < R_end; j++)
			{
				if (!t_matrix[SR[j]])
				{
					std::vector<char> t_s_connected;
					if (SR[i] > SR[j])
					{
						t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
					}
					assert(!t_s_connected.empty());
					if (t_s_connected[u] == 1)
					{
						assert((local_vc[SR[i] * n + SR[j]] > 0) && (local_vc[SR[j] * n + SR[i]] > 0));
						local_vc[SR[i] * n + SR[j]]--;
						local_vc[SR[j] * n + SR[i]]--;
					}
				}
			}
		}

		bool terminate = false;
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
				neighbors[neighbors_n++] = SR[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
			--degree_in_S[v];
			--degree[v];
			if (degree[v] + s <= best_solution_size)
			{
				if (SR_rid[v] < S_end)
					terminate = true;
				else
				{
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if (terminate)
			return false;
		return true;
	}

	bool remove_u_from_C_with_prune(ui deleted_vertex, ui S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = deleted_vertex;
		--R_end;
		swap_pos(SR_rid[u], R_end);
		level_id[u] = level;

		// Dynamically updates local vertex connectivity
		for (ui i = 0; i < R_end; i++)
		{
			char *t_matrix = matrix + SR[i] * n;
			for (ui j = i + 1; j < R_end; j++)
			{
				if (!t_matrix[SR[j]])
				{
					std::vector<char> t_s_connected;
					if (SR[i] > SR[j])
					{
						t_s_connected = all_s_connected[std::make_pair(SR[j], SR[i])];
					}
					else
					{
						t_s_connected = all_s_connected[std::make_pair(SR[i], SR[j])];
					}
					assert(!t_s_connected.empty());
					if (t_s_connected[u] == 1)
					{
						assert((local_vc[SR[i] * n + SR[j]] > 0) && (local_vc[SR[j] * n + SR[i]] > 0));
						local_vc[SR[i] * n + SR[j]]--;
						local_vc[SR[j] * n + SR[i]]--;
					}
				}
			}
		}

		bool terminate = false;
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
				neighbors[neighbors_n++] = SR[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
			--degree[v];
			if (degree[v] + s <= best_solution_size)
			{
				if (SR_rid[v] < S_end)
					terminate = true;
				else
				{
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if (terminate)
			return false;
		return true;
	}

	bool collect_removable_vertices(ui S_end, ui R_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
			if (degree[SR[i]] + s <= best_solution_size)
				return false;

		for (ui i = S_end; i < R_end; i++)
			if (level_id[SR[i]] > level)
			{
				ui v = SR[i];
				if (S_end - degree_in_S[v] >= s || degree[v] + s <= best_solution_size)
				{
					level_id[v] = level;
					Qv.push(v);
					continue;
				}
				char *t_matrix = matrix + v * n;
				for (ui j = 0; j < S_end; j++)
				{
					if (S_end - degree_in_S[SR[j]] == s && !t_matrix[SR[j]])

					{
						level_id[v] = level;
						Qv.push(v);
						break;
					}
				}
			}
		return true;
	}

	void swap_pos(ui i, ui j)
	{
		std::swap(SR[i], SR[j]);
		SR_rid[SR[i]] = i;
		SR_rid[SR[j]] = j;
	}

	ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end, ui vertex_in_S)
	{
		assert(SR_rid[vertex_in_S] < S_end);
		ui u = n, min_degree_in_S = n;
		for (ui i = S_end; i < R_end; i++)
		{
			ui v = SR[i];
			if (!matrix[vertex_in_S * n + v])
			{
				if (degree_in_S[v] < min_degree_in_S)
				{
					u = v;
					min_degree_in_S = degree_in_S[v];
				}
				else if (degree_in_S[v] == min_degree_in_S)
				{
					if (degree[v] > degree[u])
					{
						u = v;
						min_degree_in_S = degree_in_S[v];
					}
				}
			}
		}
		assert(u != n && !matrix[vertex_in_S * n + u]);
		return u;
	}
};
#endif
