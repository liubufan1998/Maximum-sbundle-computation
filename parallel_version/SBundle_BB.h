#ifndef _SBUNDLE_BB_
#define _SBUNDLE_BB_

#include "Utility.h"
#include "Timer.h"
#include "sbundle_tool.h"

class SBundle_BB
{
private:
	ui n;
	ui temp_n;
	ept *pstart;
	ept *pend;
	ui *edges;
	ept edges_cap;

	vector<ui> degree;
	vector<ui> degree2;
	vector<ui> degree_in_S;

	ui s;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;

	vector<ui> neighbors;
	vector<ui> nonneighbors;

	vector<ui> S2;
	vector<ui> SR;
	vector<ui> SR_rid;
	vector<ui> SR_remap;
	vector<ui> SR_remap_rid;
	queue<ui> Qv;
	vector<ui> level_id;
	vector<ui> t_level_id;

	ui *buf;
	ui *buf1;
	ui *buf2;
	vector<char> vis;

	vector<pair<ui, ui>> vp;
	char *matrix;
	vector<ui> local_vc;
	vector<char> s_connected;
	vector<ui> temp_deleted;
	MaximumFlow *vc;
	long long branch_count;
	map<pair<ui, ui>, int> all_not_adj_inS;
	map<pair<ui, ui>, vector<char>> all_s_connected;
	vector<ui> temp_test;

public:
	SBundle_BB()
	{
		n = 0;
		edges_cap = 0;
		pstart = NULL;
		pend = NULL;
		edges = NULL;

		best_solution = NULL;
		s = best_solution_size = _UB_ = 0;

		buf = buf1 = buf2 = NULL;
	}

	~SBundle_BB()
	{
		if (pstart != NULL)
		{
			delete[] pstart;
			pstart = NULL;
		}
		if (pend != NULL)
		{
			delete[] pend;
			pend = NULL;
		}
		if (edges != NULL)
		{
			delete[] edges;
			edges = NULL;
		}
		if (best_solution != NULL)
		{
			delete[] best_solution;
			best_solution = NULL;
		}
		if (buf != NULL)
		{
			delete[] buf;
			buf = NULL;
		}
		if (buf1 != NULL)
		{
			delete[] buf1;
			buf1 = NULL;
		}
		if (buf2 != NULL)
		{
			delete[] buf2;
			buf2 = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;
		pstart = new ept[n + 1];
		pend = new ept[n + 1];
		edges = new ui[1];
		edges_cap = 1;

		degree.resize(n + 1);
		degree_in_S.resize(n + 1);
		best_solution = new ui[n + 1];
		S2.resize(n + 1);
		SR.resize(n + 1);
		SR_rid.resize(n + 1);
		neighbors.resize(n + 1);
		nonneighbors.resize(n + 1);
		level_id.resize(n + 1);

		buf = new ui[n + 1];
		buf1 = new ui[n + 1];
		buf2 = new ui[n + 1];
		vis.resize(n + 1);
	}

	void load_graph(ui _n, const vector<pair<ui, ui>> &vp)
	{
		n = _n;
		if (vp.size() * 2 > edges_cap)
		{
			do
			{
				edges_cap *= 2;
			} while (vp.size() * 2 > edges_cap);
			delete[] edges;
			edges = new ui[edges_cap];
		}
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ept i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first < n && vp[i].second < n);
			++degree[vp[i].first];
			++degree[vp[i].second];
		}
		for (ui i = 1; i < n; i++)
		{
			degree[i] += degree[i - 1];
			pstart[i] = degree[i - 1];
		}
		pstart[0] = 0;
		for (ept i = 0; i < vp.size(); i++)
		{
			ui a = vp[i].first, b = vp[i].second;
			edges[pstart[a]++] = b;
			edges[pstart[b]++] = a;
		}
		for (ui i = n; i > 0; i--)
			pstart[i] = pstart[i - 1];
		pstart[0] = 0;
	}

	void load_graph(ui n_, ui *_pstart, ui *_pend, ui *_edges)
	{
		ept m = 0;
		n = n_;
		for (ui i = 0; i < n; i++)
			m += _pend[i] - _pstart[i];
		if (m > edges_cap)
		{
			do
			{
				edges_cap *= 2;
			} while (m > edges_cap);
			delete[] edges;
			edges = new ui[edges_cap];
		}

		m = 0;
		for (ui i = 0; i < n; i++)
		{
			pstart[i] = m;
			for (ept j = _pstart[i]; j < _pend[i]; j++)
			{
				assert(_edges[j] < n);
				edges[m++] = _edges[j];
			}
		}
		pstart[n] = m;

		printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m / 2, double(m) / n / (n - 1));
	}

	void sBundle(ui S_, ui UB_, vector<ui> &sbundle, bool must_include_0, long long &branch_node)
	{
		s = S_;
		_UB_ = UB_;
		branch_count = 0;
		if (s <= 1)
		{
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return;
		}
		best_solution_size = sbundle.size();
		ui R_end;
		initialization(R_end);
		vc = new MaximumFlow();
		vc->prepare(n, edges_cap);

		if (R_end && best_solution_size < _UB_)
			splex_BB_search(0, R_end, 1, must_include_0);

		branch_node += branch_count;
		if (vc != nullptr)
		{
			delete vc;
			vc = nullptr;
		}
		if (best_solution_size > sbundle.size())
		{
			sbundle.clear();
			for (ui i = 0; i < best_solution_size; i++)
				sbundle.push_back(best_solution[i]);
		}
	}

	void pure_sbundle_MultiBB(ui S_end, ui R_end, ui level, bool choose_zero, vector<char> &sconnected)
	{

		branch_count++;

		if (level > 1)
			reorganize_edges(S_end, R_end, level);

		if (S_end > best_solution_size && vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end, pstart, pend, edges, s))
			store_solution(S_end);
		if (R_end > best_solution_size)
		{
			if (is_kplex(R_end))
			{
				if (vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
				{
					store_solution(R_end);
				}
			}
		}

		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		assert(choose_zero == false);
		assert(R_end > best_solution_size);

		ui u = n;
		ui u2 = n;
		ui branching_case = 4;
		ui min_vc = n;
		vector<char> s_connected2(n, 0);

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			assert(level_id[SR[i]] > level);
		}
#endif
		vc->get_branching_case(SR, u, u2, s_connected2, S_end, R_end, pstart, pend, edges, min_vc, s, n);
#ifndef NDEBUG
		ui cnt_vc = 0;
		for (ui i = 0; i < R_end; i++)
		{
			if (s_connected2[SR[i]] == 1)
			{
				cnt_vc++;
			}
		}
		assert(cnt_vc == min_vc);
#endif

		assert(u != n && u2 != n && !s_connected2.empty());
		assert(SR_rid[u] >= 0 && SR_rid[u] <= R_end && SR_rid[u2] >= 0 && SR_rid[u2] <= R_end);
		assert(level_id[u] > level && level_id[u2] > level);
		if (SR_rid[u] > SR_rid[u2])
		{
			ui temp = u;
			u = u2;
			u2 = temp;
		}
		assert(SR_rid[u] < SR_rid[u2]);

		if (SR_rid[u] < S_end && SR_rid[u2] < S_end)
		{
			branching_case = 1;

			if (min_vc + s <= best_solution_size)
			{
				return;
			}
		}
		else if (SR_rid[u] < S_end && SR_rid[u2] >= S_end)
		{
			branching_case = 2;
		}
		else
		{
			assert(SR_rid[u] >= S_end && SR_rid[u2] >= S_end);
			branching_case = 3;
		}
		assert(branching_case != 4);

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif

		if (branching_case == 1)
		{

			pure_multibranch_based_on_connetivity(u, u2, S_end, R_end, level, s_connected2, min_vc);

#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui u = SR[i], cnt = 0, cnt2 = 0;
				for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
					if (SR_rid[edges[j]] < R_end)
					{
						++cnt;
						if (SR_rid[edges[j]] < S_end)
							++cnt2;
					}
				assert(degree[u] == cnt);
			}
#endif
		}
		else if (branching_case == 2)
		{

			assert(SR_rid[u] < SR_rid[u2] && SR_rid[u2] >= S_end);
			pure_simple_binary_branching1(u, u2, S_end, R_end, level, min_vc);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui u = SR[i], cnt = 0, cnt2 = 0;
				for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
					if (SR_rid[edges[j]] < R_end)
					{
						++cnt;
						if (SR_rid[edges[j]] < S_end)
							++cnt2;
					}
				assert(degree[u] == cnt);
			}
#endif
		}
		else if (branching_case == 3)
		{

			pure_simple_binary_branching2(u, u2, S_end, R_end, level);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui u = SR[i], cnt = 0, cnt2 = 0;
				for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
					if (SR_rid[edges[j]] < R_end)
					{
						++cnt;
						if (SR_rid[edges[j]] < S_end)
							++cnt2;
					}
				assert(degree[u] == cnt);
			}
#endif
		}
	}

	void sbundle_MultiBB(ui S_end, ui R_end, ui level, bool choose_zero)
	{
		branch_count++;

		if (S_end > best_solution_size)
			store_solution2(S_end);
		if (R_end > best_solution_size && is_kplex2(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR_remap, R_end, matrix, temp_n, s, true))
			store_solution2(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		vector<ui> local_vc_backup;
		map<pair<ui, ui>, vector<char>> all_s_connected_backup;

		bool need_update_local_vc = false;
		if (isSbundle_preliminary(R_end, true))
		{
			if (vc->verify_SBundle_by_MaxFlowAlg0(SR_remap, R_end, matrix, temp_n, s, true))
			{
				assert(R_end > best_solution_size);
				store_solution2(R_end);
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

		assert(R_end > best_solution_size);

		ui u = temp_n;
		ui u2 = temp_n;
		ui min_vc = R_end - s;
		int branching_case = 4;
		ui min_vc_in_S = temp_n;

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

		if (branching_case == 1)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[SR_remap[i] * temp_n + SR_remap[j]])
						++d1;
				assert(d1 == degree2[SR_remap[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(t_level_id[SR_remap[i]] > level);
#endif
			multibranch_based_on_connetivity(u, u2, S_end, R_end, level);
		}
		else if (branching_case == 2)
		{
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[SR_remap[i] * temp_n + SR_remap[j]])
						++d1;
				assert(d1 == degree2[SR_remap[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(t_level_id[SR_remap[i]] > level);
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
					if (matrix[SR_remap[i] * temp_n + SR_remap[j]])
						++d1;
				assert(d1 == degree2[SR_remap[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(t_level_id[SR_remap[i]] > level);
#endif
			simple_binary_branching_2(u, u2, S_end, R_end, level);
		}

		if (need_update_local_vc)
		{
			restore_local_vc(local_vc, local_vc_backup, all_s_connected, all_s_connected_backup);
		}
	}

	void find_branching_vertex_pair(ui &branch_vertex_1, ui &branch_vertex_2, ui pointer_1, ui pointer_2, int &branching_case, int branching_case_value, ui &min_vc, ui &min_vc_in_S)
	{
		if (branching_case_value == 1)
		{
			for (ui i = 0; i < pointer_1; i++)
			{

				char *t_matrix = matrix + SR_remap[i] * temp_n;
				for (ui j = i + 1; j < pointer_2; j++)
				{
					if (!t_matrix[SR_remap[j]])
					{
						ui temp_vc = local_vc[SR_remap[i] * temp_n * SR_remap[j]];

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
				char *t_matrix = matrix + SR_remap[i] * temp_n;

				for (ui j = pointer_1; j < pointer_2; j++)
				{
					if (!t_matrix[SR_remap[j]])
					{
						ui temp_vc = local_vc[SR_remap[i] * temp_n + SR_remap[j]];

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

				char *t_matrix = matrix + SR_remap[i] * temp_n;
				for (ui j = i + 1; j < pointer_2; j++)
				{
					if (!t_matrix[SR_remap[j]])
					{
						ui temp_vc = local_vc[SR_remap[i] * temp_n + SR_remap[j]];

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
			if (t_level_id[last_deleted] == level)
			{
				addtoC(last_deleted, R_end, level);
				temp_deleted.pop_back();
				count_restore++;
			}
			else
				break;
		}
	}

	void restore_deleted_vertices2(ui &R_end, ui level, ui &count_restore)
	{
		while (!temp_deleted.empty())
		{
			ui last_deleted = temp_deleted.back();
			if (level_id[last_deleted] == level)
			{
				addtoC2(last_deleted, R_end, level);
				temp_deleted.pop_back();
				count_restore++;
			}
			else
			{
				break;
			}
		}
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif
	}

	void update_local_vc(vector<ui> &local_vc, vector<ui> &local_vc_backup, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected, std::map<std::pair<ui, ui>, std::vector<char>> &all_s_connected_backup, ui R_end)
	{
		local_vc_backup = local_vc;
		local_vc.clear();
		all_s_connected_backup = all_s_connected;
		all_s_connected.clear();
		vc->get_local_vertexConnectivity(R_end, SR_remap, matrix, temp_n, local_vc, all_s_connected, s);
	}

	void restore_local_vc(vector<ui> &local_vc, vector<ui> &local_vc_backup, map<pair<ui, ui>, vector<char>> &all_s_connected, map<pair<ui, ui>, vector<char>> all_s_connected_backup)
	{
		local_vc = local_vc_backup;
		local_vc_backup.clear();
		all_s_connected.clear();
		all_s_connected = all_s_connected_backup;
	}

	void pure_multibranch_based_on_connetivity(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level, vector<char> &s_connected2, ui min_vc)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert((SR_rid[u] < S_end && SR_rid[u2] < S_end));
		assert(!s_connected2.empty());

		ui notadj_S = 0;
		vector<ui> notadj_C;
		assert(min_vc + s > best_solution_size);
		bool should_terminate = false;
		get_notadj_S_and_notadj_C2(notadj_S, notadj_C, S_end, R_end, level, s_connected2, min_vc, should_terminate);
		if (should_terminate)
		{
			return;
		}
		assert(!notadj_C.empty() && 0 < notadj_S <= s);

		ui canselect = s - notadj_S;
		ui pos = n;

		bool test = canselect < notadj_C.size();
		if (!test)
		{
			printf("当前的canselect的大小是%u, 当前的notadj_S的大小是%u，当前的notadj_C的大小是%u\n", canselect, notadj_S, notadj_C.size());
			printf("当前的should_terminate的值是：%u\n", should_terminate);
		}
		assert(canselect >= 0 && canselect < notadj_C.size());

		bool all = true;
		generate_canselect_branches2(all, canselect, S_end, R_end, level, pos, notadj_C, u, u2, s_connected2);

		generate_the_last_branch2(all, canselect, S_end, R_end, level, notadj_C, u, u2, s_connected2);

		if (pos != n)
		{
			for (ui i = pos; i >= 0 && i <= notadj_C.size(); i--)
			{
				delfrS2(notadj_C[i], S_end);
			}
		}
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif
	}

	void pure_simple_binary_branching1(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level, ui min_vc)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;

		assert((SR_rid[u] < S_end && SR_rid[u2] >= S_end));
		ui record_index = u2;
		assert(level_id[record_index] > level);

		if (min_vc + s > best_solution_size)
		{
			swap_pos(S_end, SR_rid[record_index]);
			bool canadd = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
			swap_pos(S_end, SR_rid[record_index]);
			if (canadd)
			{
				addtoS2(record_index, S_end);
				ui count_deleted = 0, count_restore = 0;
				delete_vertices_from_C2(S_end, R_end, level, count_deleted);
				pure_sbundle_MultiBB(S_end, R_end, level + 1, false, s_connected);
				restore_deleted_vertices2(R_end, level, count_restore);
				assert(count_restore == count_deleted);
				delfrS2(record_index, S_end);
			}
		}

		assert(SR_rid[record_index] < R_end && SR_rid[record_index] >= S_end);
		delfrC2(record_index, R_end, level);
		pure_sbundle_MultiBB(S_end, R_end, level + 1, false, s_connected);
		addtoC2(record_index, R_end, level);
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif
	}

	void pure_simple_binary_branching2(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;

		assert((SR_rid[u] < S_end && SR_rid[u2] >= S_end) || (SR_rid[u] >= S_end && SR_rid[u2] >= S_end));
		ui record_index = u2;
		assert(level_id[record_index] > level);

		swap_pos(S_end, SR_rid[record_index]);
		bool canadd = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
		swap_pos(S_end, SR_rid[record_index]);
		if (canadd)
		{
			addtoS2(record_index, S_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C2(S_end, R_end, level, count_deleted);
			pure_sbundle_MultiBB(S_end, R_end, level + 1, false, s_connected);
			restore_deleted_vertices2(R_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS2(record_index, S_end);
		}

		assert(SR_rid[record_index] < R_end && SR_rid[record_index] >= S_end);
		delfrC2(record_index, R_end, level);
		pure_sbundle_MultiBB(S_end, R_end, level + 1, false, s_connected);
		addtoC2(record_index, R_end, level);
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif
	}
	void multibranch_based_on_connetivity(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u < S_end && u2 < S_end);

		find_s_connected(s_connected, u, u2);
		assert(!s_connected.empty());

		ui notadj_S = 0;
		vector<ui> notadj_C;
		get_notadj_S_and_notadj_C(notadj_S, notadj_C, S_end, R_end, level);
		assert(!notadj_C.empty() && notadj_S <= s);

		assert(s >= notadj_S);
		ui canselect = s - notadj_S;
		ui pos = temp_n;
		assert(canselect >= 0 && canselect < notadj_C.size());

		bool all = true;
		generate_canselect_branches(all, canselect, S_end, R_end, level, pos, notadj_C);

		generate_the_last_branch(all, canselect, S_end, R_end, level, notadj_C);

		if (pos != temp_n)
		{
			for (ui i = pos; i >= 0 && i <= notadj_C.size(); i--)
			{
				delfrS(notadj_C[i], S_end, R_end);
			}
		}
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
	}

	void simple_binary_branching_1(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level, ui min_vc)
	{
		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;

		assert(u < S_end && u2 >= S_end);
		ui record_index = SR_remap[u2];

		if (min_vc + s > best_solution_size && canadd(record_index, S_end, true))
		{
			assert(t_level_id[record_index] > level);
			addtoS(record_index, S_end, R_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C(S_end, R_end, level, count_deleted);
			sbundle_MultiBB(S_end, R_end, level + 1, false);
			restore_deleted_vertices(R_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS(record_index, S_end, R_end);
		}

		assert(SR_remap_rid[record_index] < R_end && SR_remap_rid[record_index] >= S_end);
		delfrC(record_index, R_end, level);
		sbundle_MultiBB(S_end, R_end, level + 1, false);
		addtoC(record_index, R_end, level);
		ui after_R_end = R_end;
		assert(pre_Rend == after_R_end);
	}

	void simple_binary_branching_2(ui branching_vertex_1, ui branching_vertex_2, ui S_end, ui R_end, ui level)
	{

		ui pre_Rend = R_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u >= S_end && u2 >= S_end);
		ui record_index = SR_remap[u2];
		assert(t_level_id[record_index] > level && canadd(record_index, S_end, true));
		if (canadd(record_index, S_end, true))
		{
			addtoS(record_index, S_end, R_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C(S_end, R_end, level, count_deleted);
			sbundle_MultiBB(S_end, R_end, level + 1, false);
			restore_deleted_vertices(R_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS(record_index, S_end, R_end);
		}

		assert(SR_remap_rid[record_index] < R_end && SR_remap_rid[record_index] >= S_end);
		delfrC(record_index, R_end, level);
		sbundle_MultiBB(S_end, R_end, level + 1, false);
		addtoC(record_index, R_end, level);
		assert(pre_Rend == R_end);
	}

	void generate_canselect_branches(bool &all, ui canselect, ui &S_end, ui &R_end, ui level, ui &pos, vector<ui> notadj_C)
	{
		for (ui i = 0; i < canselect; ++i)
		{
			assert(SR_remap_rid[notadj_C[i]] < R_end && SR_remap_rid[notadj_C[i]] >= S_end);
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
				addtoS(notadj_C[i - 1], S_end, R_end);
				delete_vertices_from_C(S_end, R_end, level, count_deleted);

				pos = i - 1;
			}

			sbundle_MultiBB(S_end, R_end, level + 1, false);

			ui count_restore = 0;
			if (i)
			{
				restore_deleted_vertices(R_end, level, count_restore);
			}
			assert(count_restore == count_deleted);
			addtoC(notadj_C[i], R_end, level);
		}
	}

	void generate_canselect_branches2(bool &all, ui canselect, ui &S_end, ui &R_end, ui level, ui &pos, vector<ui> &notadj_C, ui u1, ui u2, vector<char> &sconnected)
	{
		for (ui i = 0; i < canselect; ++i)
		{
			assert(SR_rid[notadj_C[i]] < R_end && SR_rid[notadj_C[i]] >= S_end);
			assert(level_id[notadj_C[i]] > level);

			delfrC2(notadj_C[i], R_end, level);

			bool canadd = false;
			if (i)
			{
				swap_pos(S_end, notadj_C[i - 1]);
				canadd = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
				swap_pos(S_end, notadj_C[i - 1]);
			}
			if (i && !canadd)
			{
				addtoC2(notadj_C[i], R_end, level);
				all = false;
				break;
			}
			ui count_deleted = 0;
			if (i)
			{
				assert(canadd);

				addtoS2(notadj_C[i - 1], S_end);
				delete_vertices_from_C2(S_end, R_end, level, count_deleted);

				pos = i - 1;
				assert(pos < n);
			}

			pure_sbundle_MultiBB(S_end, R_end, level + 1, false, sconnected);

			ui count_restore = 0;
			if (i)
			{
				restore_deleted_vertices2(R_end, level, count_restore);
			}
			assert(count_restore == count_deleted);

			addtoC2(notadj_C[i], R_end, level);
		}
	}

	void generate_the_last_branch(bool all, ui canselect, ui &S_end, ui &R_end, ui level, vector<ui> notadj_C)
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

				assert(t_level_id[notadj_C[i]] > level);
				if (t_level_id[notadj_C[i]] > level)
				{
					assert(SR_remap_rid[notadj_C[i]] < R_end && SR_remap_rid[notadj_C[i]] >= S_end);
					delfrC(notadj_C[i], R_end, level);
					pre_dd++;
				}
			}

			ui count_deleted = 0;
			if (canselect)
			{

				assert(t_canadd);
				addtoS(notadj_C[canselect - 1], S_end, R_end);
				delete_vertices_from_C(S_end, R_end, level, count_deleted);
			}

			sbundle_MultiBB(S_end, R_end, level + 1, false);

			ui count_restore = 0;
			if (canselect)
			{
				restore_deleted_vertices(R_end, level, count_restore);
				delfrS(notadj_C[canselect - 1], S_end, R_end);
			}
			assert(count_restore == count_deleted);

			/*20240306：发现的第二个问题就是下面的这一些代码应该是放在if(canselect)外面的，而我之前是放在了里面。
			20240321：我猜之前在候选集中更新顶点的时候可能出错的一个原因是我下面的添加回顶点的时候出了问题。因为前面我是从canselect开始删除的，而正常来讲在恢复的时候我应该从数组的末尾进行恢复。*/
			ui after_dd = 0;
			for (ui i = notadj_C.size() - 1; i >= canselect && i <= notadj_C.size(); i--)
			{
				if (t_level_id[notadj_C[i]] == level)
				{
					addtoC(notadj_C[i], R_end, level);
					after_dd++;
				}
			}
			assert(pre_dd == after_dd);
		}
	}

	void generate_the_last_branch2(bool all, ui canselect, ui &S_end, ui &R_end, ui level, vector<ui> notadj_C, ui u1, ui u2, vector<char> &sconnected)
	{

		bool canadd = false;
		if (canselect != 0)
		{
			assert(level_id[notadj_C[canselect - 1]] > level);
			assert(SR_rid[notadj_C[canselect - 1]] >= S_end && SR_rid[notadj_C[canselect - 1]] < R_end);
			swap_pos(S_end, SR_rid[notadj_C[canselect - 1]]);
			canadd = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
			swap_pos(S_end, SR_rid[notadj_C[canselect - 1]]);
		}

		if (all && (canselect == 0 || canadd))
		{
			ui pre_dd = 0;
			for (ui i = canselect; i < (ui)notadj_C.size(); ++i)
			{

				assert(level_id[notadj_C[i]] > level);
				if (level_id[notadj_C[i]] > level)
				{
					assert(SR_rid[notadj_C[i]] < R_end && SR_rid[notadj_C[i]] >= S_end);

					delfrC2(notadj_C[i], R_end, level);
					pre_dd++;
				}
			}

			ui count_deleted = 0;
			if (canselect)
			{
				assert(canadd);

				addtoS2(notadj_C[canselect - 1], S_end);
				delete_vertices_from_C2(S_end, R_end, level, count_deleted);
			}

			pure_sbundle_MultiBB(S_end, R_end, level + 1, false, sconnected);

			ui count_restore = 0;
			if (canselect)
			{
				restore_deleted_vertices2(R_end, level, count_restore);

				delfrS2(notadj_C[canselect - 1], S_end);
			}
			assert(count_restore == count_deleted);

			/*20240306：发现的第二个问题就是下面的这一些代码应该是放在if(canselect)外面的，而我之前是放在了里面。
			20240321：我猜之前在候选集中更新顶点的时候可能出错的一个原因是我下面的添加回顶点的时候出了问题。因为前面我是从canselect开始删除的，而正常来讲在恢复的时候我应该从数组的末尾进行恢复。*/
			ui after_dd = 0;
			for (ui i = notadj_C.size() - 1; i >= canselect && i <= notadj_C.size(); i--)
			{
				if (level_id[notadj_C[i]] == level)
				{
					addtoC2(notadj_C[i], R_end, level);
					after_dd++;
				}
			}
			assert(pre_dd == after_dd);
		}
	}

	void get_notadj_S_and_notadj_C(ui &notadj_S, vector<ui> &notadj_C, ui S_end, ui R_end, ui level)
	{
		assert(notadj_C.empty());
		for (ui i = 0; i < R_end; i++)
		{
			if (s_connected[SR_remap[i]] == 0)
			{
				if (i < S_end)
				{
					notadj_S++;
				}
				else
				{
					assert(t_level_id[SR_remap[i]] > level);
					notadj_C.push_back(SR_remap[i]);
				}
			}
		}
#ifdef _BETTER_SEQEUNCE_
		ui t_n = notadj_C.size();
		vector<ui> scores(t_n);

		for (ui i = 0; i < t_n; ++i)
		{
			ui u = notadj_C[i];
			scores[i] = degree2[u];
		}

		vector<ui> indices(t_n);
		for (ui i = 0; i < t_n; ++i)
		{
			indices[i] = i;
		}

		sort(indices.begin(), indices.end(), [&scores](ui i, ui j)
			 { return scores[i] < scores[j]; });

		vector<ui> sorted_arr1(t_n);
		for (ui i = 0; i < t_n; ++i)
		{
			sorted_arr1[i] = notadj_C[indices[i]];
		}

		notadj_C = sorted_arr1;
#endif
		assert(notadj_S <= s);
		assert(notadj_C.size() > s - notadj_S);
	}

	void get_notadj_S_and_notadj_C2(ui &notadj_S, vector<ui> &notadj_C, ui S_end, ui R_end, ui level, vector<char> &s_connected, ui min_vc, bool &terminate)
	{
		assert(notadj_C.empty());
		assert(s_connected.size() == n);
		for (ui i = 0; i < R_end; i++)
		{
			assert(level_id[SR[i]] > level);
			if (s_connected[SR[i]] == 0)
			{
				if (i < S_end)
				{
					notadj_S++;
				}
				else
				{
					notadj_C.push_back(SR[i]);
				}
			}
		}
		if (notadj_S > s)
		{
			terminate = true;
			return;
		}
		assert(notadj_S <= s);
		assert(notadj_C.size() > s - notadj_S);
	}

	void delete_vertices_from_C(ui S_end, ui &R_end, ui level, ui &count_deleted)
	{
		for (ui i = S_end; i < R_end; i++)
		{
			ui temp_record = SR_remap[i];
			if (t_level_id[temp_record] > level)
			{
				bool can_delete = false;
				if (degree2[temp_record] < best_solution_size - s)
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
					assert(t_level_id[temp_record] > level);
					assert(SR_remap_rid[temp_record] < R_end && SR_remap_rid[temp_record] >= S_end);
					delfrC(temp_record, R_end, level);
					count_deleted++;
				}
			}
		}
	}

	void delete_vertices_from_C2(ui S_end, ui &R_end, ui level, ui &count_deleted)
	{
		for (ui i = S_end; i < R_end; i++)
		{
			ui temp_record = SR[i];
			if (level_id[temp_record] > level)
			{
				bool can_delete = false;
				if (degree[temp_record] <= best_solution_size - s)
				{
					can_delete = true;
				}
				if (can_delete)
				{
					temp_deleted.push_back(temp_record);
					assert(level_id[temp_record] > level);
					assert(SR_rid[temp_record] < R_end && SR_rid[temp_record] >= S_end);
					delfrC2(temp_record, R_end, level);
					count_deleted++;
				}
			}
		}
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
		}
#endif
	}

	void find_s_connected(vector<char> &s_connected, ui vertex_1, ui vertex_2)
	{
		ui u = vertex_1;
		ui u2 = vertex_2;
		s_connected.clear();
		if (SR_remap[u] > SR_remap[u2])
		{
			s_connected = all_s_connected[make_pair(SR_remap[u2], SR_remap[u])];
		}
		else
		{
			s_connected = all_s_connected[make_pair(SR_remap[u], SR_remap[u2])];
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
				char *t_matrix = matrix + SR_remap[i] * temp_n;
				for (ui j = i + 1; j < end; j++)
				{
					if (!t_matrix[SR_remap[j]] && local_vc[SR_remap[i] * temp_n + SR_remap[j]] < end - s)
						return false;
				}
			}
		}
		return true;
	}

	ui get_vertex_pair_ub(ui S_end, ui R_end)
	{
		assert(S_end >= 2);
		all_not_adj_inS.clear();

		for (ui i = 0; i < S_end; i++)
		{
			char *t_matrix = matrix + SR_remap[i] * temp_n;
			for (ui j = i + 1; j < S_end; j++)
			{
				if (!t_matrix[SR_remap[j]])
				{
					if (SR_remap[i] > SR_remap[j])
					{
						all_not_adj_inS[make_pair(SR_remap[j], SR_remap[i])] = 1;
					}
					else
					{
						all_not_adj_inS[make_pair(SR_remap[i], SR_remap[j])] = 1;
					}
				}
			}
		}

		if (all_not_adj_inS.size() == 0)
			return R_end;

		ui ret_ub = S_end;

		assert(!all_not_adj_inS.empty());
		ui t_size = all_not_adj_inS.size();
		ui lost_in_S[t_size] = {};
		ui lost_in_C[t_size] = {};

		double dise[t_size] = {};
		ui cost[t_size] = {};
		vector<vector<ui>> vp_partition(t_size, vector<ui>(R_end));
		ui connect_all_vp_in_S = R_end - S_end;
		vector<ui> is_all_sconnect_with_S(R_end, 1);
		ui vp_cnt = 0;

		for (auto it = all_not_adj_inS.begin(); it != all_not_adj_inS.end(); ++it)
		{
			pair<ui, ui> cur_vp = it->first;
			vector<char> cur_s_con = all_s_connected[cur_vp];
			assert(!cur_s_con.empty());
			if (!cur_s_con.empty())
			{
				for (ui j = 0; j < R_end; j++)
				{
					if (cur_s_con[SR_remap[j]] == 0)
					{
						if (j < S_end)
						{
							lost_in_S[vp_cnt]++;
						}
						else
						{
							lost_in_C[vp_cnt]++;
							vp_partition[vp_cnt][j] = 1;
							if (is_all_sconnect_with_S[j] == 1)
							{
								is_all_sconnect_with_S[j] = 0;
								assert(connect_all_vp_in_S > 0);
								connect_all_vp_in_S--;
							}
						}
					}
				}

				if (s - lost_in_S[vp_cnt] > lost_in_C[vp_cnt])
				{
					cost[vp_cnt] = lost_in_C[vp_cnt];
				}
				else
				{
					cost[vp_cnt] = s - lost_in_S[vp_cnt];
				}

				if (cost[vp_cnt] == 0)
				{
					dise[vp_cnt] = 0;
				}
				else
				{
					dise[vp_cnt] = lost_in_C[vp_cnt] / static_cast<double>(cost[vp_cnt]);
				}
				vp_cnt++;
			}
		}
		if (vp_cnt != t_size)
		{
			printf("当前的t_size的大小是%u, 以及当前的vp_cnt的大小是%u ", t_size, vp_cnt);
		}
		assert(vp_cnt == t_size);

		ret_ub = ret_ub + connect_all_vp_in_S;
		if (ret_ub > best_solution_size)
		{
			return R_end;
		}

		vector<bool> visited(t_size, false);
		ui sum_visit = t_size;
		while (sum_visit > 0)
		{
			ui max_dise = 0;

			ui impossible_number = temp_n * temp_n * temp_n;
			ui max_dise_index = impossible_number;
			for (ui i = 0; i < t_size; i++)
			{
				if (!visited[i])
				{
					if (dise[i] > max_dise)
					{
						max_dise = dise[i];
						max_dise_index = i;
					}
					if (dise[i] == 0)
					{
						visited[i] = true;
						assert(sum_visit > 0);
						sum_visit--;
					}
				}
			}

			if (sum_visit == 0)
			{
				break;
			}
			assert(max_dise_index != impossible_number);
			visited[max_dise_index] = true;
			assert(sum_visit > 0);
			sum_visit--;

			ret_ub = ret_ub + cost[max_dise_index];
			if (ret_ub > best_solution_size)
			{
				return R_end;
			}

			for (ui i = S_end; i < R_end; i++)
			{
				if (vp_partition[max_dise_index][i] == 1)
				{
					for (ui j = 0; j < t_size; j++)
					{
						if (!visited[j])
						{
							if (vp_partition[j][i] == 1)
							{
								assert(lost_in_C[j] > 0);
								lost_in_C[j]--;
								if (lost_in_C[j] == 0)
								{
									visited[j] = true;
									assert(sum_visit > 0);
									sum_visit--;
									dise[j] = 0;
									continue;
								}

								if (s - lost_in_S[j] > lost_in_C[j])
								{
									cost[j] = lost_in_C[j];
								}
								else
								{
									cost[j] = s - lost_in_S[j];
								}

								if (cost[j] == 0)
								{
									dise[j] = 0;
								}
								else
								{
									dise[j] = lost_in_C[j] / static_cast<double>(cost[j]);
								}
								vp_partition[j][i] = 0;
							}
						}
					}
				}
			}
		}

		if (ret_ub > R_end)
		{
			printf("当前的ret_ub的大小是：%u\n", ret_ub);
			printf("当前的R_end的大小是：%u\n", R_end);
		}
		assert(ret_ub <= R_end);

		return ret_ub;
	}

	bool canadd(ui u, ui S_end, bool precise)
	{
		char *t_matrix = matrix + u * temp_n;
		for (ui i = 0; i < S_end; i++)
		{
			if (!t_matrix[SR_remap[i]])
			{

				if (local_vc[u * temp_n + SR_remap[i]] <= best_solution_size - s)
				{
					return false;
				}
			}
		}

		if (precise)
		{
			if (SR_remap_rid[u] != S_end)
			{
				swap_pos2(S_end, SR_remap_rid[u]);
			}
			return vc->verify_SBundle_by_MaxFlowAlg0(SR_remap, S_end + 1, matrix, temp_n, s, false);
		}

		return true;
	}

	bool removeByVpUB(ui u, ui S_end, ui R_end)
	{
		addtoS(u, S_end, R_end);
		ui vp_bound = get_vertex_pair_ub(S_end, R_end);
		if (vp_bound > best_solution_size)
		{
			delfrS(u, S_end, R_end);
			return false;
		}
		delfrS(u, S_end, R_end);
		return true;
	}

	void delfrC(ui u, ui &R_end, ui level)
	{
		assert(t_level_id[u] == temp_n);
		t_level_id[u] = level;
		--R_end;
		swap_pos2(R_end, SR_remap_rid[u]);

		char *t_matrix = matrix + u * temp_n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR_remap[i]])
			{
				assert(degree2[SR_remap[i]] > 0);
				--degree2[SR_remap[i]];
			}

		for (ui i = 0; i < R_end; i++)
		{
			char *t_matrix = matrix + SR_remap[i] * temp_n;
			for (ui j = i + 1; j < R_end; j++)
			{
				if (!t_matrix[SR_remap[j]])
				{
					vector<char> t_s_connected;
					if (SR_remap[i] > SR_remap[j])
					{
						t_s_connected = all_s_connected[make_pair(SR_remap[j], SR_remap[i])];
					}
					else
					{
						t_s_connected = all_s_connected[make_pair(SR_remap[i], SR_remap[j])];
					}
					assert(!t_s_connected.empty());
					if (t_s_connected[u] == 1)
					{
						assert((local_vc[SR_remap[i] * temp_n + SR_remap[j]] > 0) && (local_vc[SR_remap[j] * temp_n + SR_remap[i]] > 0));
						local_vc[SR_remap[i] * temp_n + SR_remap[j]]--;
						local_vc[SR_remap[j] * temp_n + SR_remap[i]]--;
					}
				}
			}
		}
	}

	void addtoC(ui u, ui &R_end, ui level)
	{
		assert(t_level_id[u] != temp_n && t_level_id[u] == level);
		t_level_id[u] = temp_n;
		{
			for (ui i = 0; i < R_end; i++)
			{
				char *t_matrix = matrix + SR_remap[i] * temp_n;
				for (ui j = i + 1; j < R_end; j++)
				{
					if (!t_matrix[SR_remap[j]])
					{
						vector<char> t_s_connected;
						if (SR_remap[i] > SR_remap[j])
						{
							t_s_connected = all_s_connected[make_pair(SR_remap[j], SR_remap[i])];
						}
						else
						{
							t_s_connected = all_s_connected[make_pair(SR_remap[i], SR_remap[j])];
						}
						assert(!t_s_connected.empty());
						if (t_s_connected[u] == 1)
						{
							assert((local_vc[SR_remap[i] * temp_n + SR_remap[j]] < temp_n) && (local_vc[SR_remap[j] * temp_n + SR_remap[i]] < temp_n));
							local_vc[SR_remap[i] * temp_n + SR_remap[j]]++;
							local_vc[SR_remap[j] * temp_n + SR_remap[i]]++;
						}
					}
				}
			}

			for (ui i = 0; i < R_end; i++)
			{
				char *t_matrix = matrix + SR_remap[i] * temp_n;
				if (!t_matrix[u])
				{
					int v = SR_remap[i];
					int cnt = 0;
					vector<char> t_s_connected;
					if (u > v)
					{
						t_s_connected = all_s_connected[make_pair(v, u)];
					}
					else
					{
						t_s_connected = all_s_connected[make_pair(u, v)];
					}
					assert(!t_s_connected.empty());
					for (ui j = 0; j < R_end; j++)
					{
						int w = SR_remap[j];
						if (t_s_connected[w])
							cnt++;
					}
					if (t_s_connected[u])
						cnt++;
					local_vc[temp_n * v + u] = local_vc[temp_n * u + v] = cnt;
				}
			}
		}

		SR_remap_rid[u] = R_end;
		SR_remap[R_end] = u;
		++R_end;

		char *t_matrix = matrix + u * temp_n;
		degree2[u] = 0;
		for (ui i = 0; i < R_end; i++)
		{
			if (t_matrix[SR_remap[i]])
			{
				++degree2[SR_remap[i]];
				degree2[u]++;
			}
		}
	}

	void addtoS(ui u, ui &S_end, ui R_end)
	{
		swap_pos2(S_end, SR_remap_rid[u]);
		++S_end;
	}

	void delfrS(ui u, ui &S_end, ui R_end)
	{
		--S_end;
		assert(S_end == SR_remap_rid[u]);
	}

	void delfrC2(ui u, ui &R_end, ui level)
	{
		assert(level_id[u] == n && SR[SR_rid[u]] == u);
		assert(SR_rid[u] < R_end);
		level_id[u] = level;

		--R_end;
		swap_pos(R_end, SR_rid[u]);

		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < R_end)
			{
				ui w = edges[i];
				assert(degree[w] > 0);
				--degree[w];
			}
	}

	void addtoC2(ui u, ui &R_end, ui level)
	{
		assert(level_id[u] == level && SR_rid[u] == R_end);
		level_id[u] = n;
		++R_end;
		degree[u] = 0;

		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
		{
			ui w = edges[i];
			assert(SR[SR_rid[w]] == w);
			if (SR_rid[w] < R_end)
			{
				assert(degree[w] < n && degree[u] < n);
				++degree[w];
				++degree[u];
			}
		}
	}

	void addtoS2(ui u, ui &S_end)
	{
		swap_pos(S_end, SR_rid[u]);
		++S_end;
	}

	void delfrS2(ui u, ui &S_end)
	{
		--S_end;
		assert(S_end == SR_rid[u]);
	}

private:
	void initialization(ui &R_end)
	{
		vis.clear();
		degree_in_S.clear();
		for (ui i = 0; i < n; i++)
			level_id[i] = n;
		for (ui i = 0; i < n; i++)
			SR[i] = SR_rid[i] = i;
		for (ui i = 0; i < n; i++)
		{
			pend[i] = pstart[i + 1];
		}

		while (!Qv.empty())
			Qv.pop();

		for (ui i = 0; i < n; i++)
		{
			degree[i] = pstart[i + 1] - pstart[i];
			if (degree[i] + s <= best_solution_size)
			{
				level_id[i] = 0;
				Qv.push(i);
			}
		}
		R_end = n;
		if (!remove_vertices_with_prune(0, R_end, 0))
			R_end = 0;
	}

	void reorganize_edges(ui S_end, ui R_end, ui level)
	{
		assert(level > 0);
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			assert(level_id[u] > level);
			ui non_neighbors_n = 0, end = pstart[u];
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level - 1; j++)
			{
				assert(level_id[edges[j]] == level - 1 || level_id[edges[j]] == n);
				if (level_id[edges[j]] >= level)
				{
					edges[end++] = edges[j];
				}
				else
					nonneighbors[non_neighbors_n++] = edges[j];
			}

			// if (degree[u] != end - pstart[u])
			// {
			// 	printf("Something wrong 333!\n");
			// 	printf("当前的顶点u是%u\n", u);
			// 	printf("当前的level值是%u，当前顶点u的level值是%u\n", level, level_id[u]);
			// 	printf("当前的degree[u]是%u, 当前的end是%u，当前的pstart[u]是%u，当前的end-pstart[u]是%u\n", degree[u], end, pstart[u], end - pstart[u]);
			// }
			assert(degree[u] == end - pstart[u]);
			for (ui j = 0; j < non_neighbors_n; j++)
				edges[end++] = nonneighbors[j];
			assert((end < pend[u] && level_id[edges[end]] < level - 1) || end == pend[u]);
#ifndef NDEBUG
			for (ui j = end; j < pend[u]; j++)
			{
				if (level_id[edges[j]] >= level)
					printf("removed_level[edges[j]]: %u, level: %u\n", level_id[edges[j]], level);
				assert(level_id[edges[j]] < level);
			}
#endif
		}
	}

	void compute_a_heuristic_solution_and_prune(ui &R_end, ui level)
	{

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			assert(degree[u] + s > best_solution_size);
			ui end = pstart[u];
			while (end < pend[u] && level_id[edges[end]] >= level)
				++end;
#ifndef NDEBUG
			for (ui j = end; j < pend[u]; j++)
			{
				if (level_id[edges[j]] >= level)
					printf("removed_level[edges[j]]: %u, level: %u\n", level_id[edges[j]], level);
				assert(level_id[edges[j]] < level);
			}
#endif
			if (degree[u] != end - pstart[u])
				printf("degree[u]: %u, %u\n", degree[u], end - pstart[u]);
			assert(degree[u] == end - pstart[u]);
		}
#endif
		vector<ui> core = neighbors;
		vector<ui> rid = nonneighbors;
		ui *id = buf;
		ui *t_degree = buf1;
		for (ui i = 0; i < R_end; i++)
		{
			id[i] = 0;
			t_degree[SR[i]] = degree[SR[i]];
		}
		for (ui i = 0; i < R_end; i++)
			++id[t_degree[SR[i]]];
		for (ui i = 1; i < R_end; i++)
			id[i] += id[i - 1];

		for (ui i = 0; i < R_end; i++)
			rid[SR[i]] = --id[t_degree[SR[i]]];
		for (ui i = 0; i < R_end; i++)
			id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for (ui i = 0, j = 0; i <= R_end; i++)
		{
			while (j < R_end && t_degree[id[j]] < i)
				++j;
			degree_start[i] = j;
		}

		ui max_core = 0, pre_solution_size = best_solution_size;
		for (ui i = 0; i < R_end; i++)
		{
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if (t_degree[u] > max_core)
				max_core = t_degree[u];
			core[u] = max_core;

			if (t_degree[u] + s >= R_end - i && R_end - i > best_solution_size)
			{
				best_solution_size = R_end - i;
				for (ui j = i; j < R_end; j++)
					best_solution[j - i] = id[j];
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			++degree_start[t_degree[u]];
			if (t_degree[u] == 0)
				continue;

			degree_start[t_degree[u] - 1] = degree_start[t_degree[u]];
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (rid[edges[j]] > i)
				{
					ui v = edges[j];
					ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
					swap(id[pos1], id[pos2]);
					rid[id[pos1]] = pos1;
					rid[id[pos2]] = pos2;
					++degree_start[t_degree[v]];
					--t_degree[v];
				}
		}

		assert(Qv.empty());
		for (ui i = 0; i < R_end; i++)
			if (core[SR[i]] + s <= best_solution_size)
			{
				assert(level_id[SR[i]] > level);
				level_id[SR[i]] = level;
				Qv.push(SR[i]);
			}
		remove_vertices_with_prune(0, R_end, level);
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

	void store_solution2(ui size)
	{
		if (size <= best_solution_size)
		{
			printf("!!! the solution to store is no larger than the current best solution!");
			return;
		}
		best_solution_size = size;
		for (ui i = 0; i < best_solution_size; i++)
			best_solution[i] = SR_remap[i];
	}

	bool is_kplex(ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
			if (degree[SR[i]] + s < R_end)
				return false;
		return true;
	}

	bool is_kplex2(ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
			if (degree2[SR_remap[i]] + s < R_end)
				return false;
		return true;
	}

	void splex_BB_search(ui S_end, ui R_end, ui level, bool choose_zero)
	{

		if (S_end > best_solution_size && vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end, pstart, pend, edges, s))
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		branch_count++;

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
		{
			assert(degree[SR[i]] + s > best_solution_size);
			assert(degree_in_S[SR[i]] + s >= S_end);
		}
#endif

		ui old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());

		if (level > 1)
			reorganize_edges(S_end, R_end, level);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
			assert(level_id[u] > level);
		}
#endif
		assert(choose_zero == false);

#ifdef _SUPPORT_BOUND_
		ui S2_n = 0;
		for (ui i = 0; i < S_end; i++)
			if (R_end - degree[SR[i]] > s)
				S2[S2_n++] = SR[i];
		if (S2_n >= 2)
		{
			collect_removable_vertices_based_on_total_edges(S2_n, S_end, R_end, level);
			if (!remove_vertices_with_prune(S_end, R_end, level))
			{
				restore_SR(S_end, R_end, old_S_end, old_R_end, level);
				return;
			}
		}
#endif
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		ui u = n;
		ui min_deg = n;
		for (ui i = 0; i < R_end; i++)
		{
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}
		}
		assert(u != n && SR_rid[u] < R_end);

#ifdef _SBUNDLE_BRANCHING_

		if (min_deg >= R_end - s)
		{

			temp_n = R_end;
			SR_remap.resize(temp_n);
			SR_remap_rid.resize(temp_n);
			for (ui i = 0; i < temp_n; ++i)
			{
				SR_remap[i] = i;
				SR_remap_rid[i] = i;
			}

			matrix = new char[temp_n * temp_n];
			memset(matrix, 0, sizeof(char) * temp_n * temp_n);

			degree2.assign(R_end, 0);
			for (ui i = 0; i < R_end; ++i)
			{
				ui t_u = SR[i];
				for (ui j = i + 1; j < R_end; ++j)
				{
					ui t_v = SR[j];
					assert(level_id[t_v] > level);
					for (ui k = pstart[t_u]; k < pend[t_u]; ++k)
					{
						if (edges[k] == t_v)
						{
							matrix[SR_remap[i] * temp_n + SR_remap[j]] = matrix[SR_remap[j] * temp_n + SR_remap[i]] = 1;
							degree2[SR_remap[i]]++;
							degree2[SR_remap[j]]++;
							break;
						}
					}
				}
			}

			assert(R_end == temp_n);
			t_level_id.resize(temp_n);
			for (ui i = 0; i < temp_n; i++)
			{
				t_level_id[SR_remap[i]] = temp_n;
			}

			local_vc.resize(temp_n * temp_n, 0);
			all_s_connected.clear();
			vc->get_local_vertexConnectivity(R_end, SR_remap, matrix, temp_n, local_vc, all_s_connected, s);

			ui pre_R_end = R_end, pre_S_end = S_end;
			sbundle_MultiBB(S_end, R_end, 1, false);
			assert(R_end == pre_R_end && S_end == pre_S_end);

			degree2.clear();

			if (matrix != nullptr)
			{
				delete[] matrix;
				matrix = nullptr;
			}
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}
#else

		if (min_deg >= R_end - s)
		{

			ui pre_S_end = S_end, pre_R_end = R_end;
			sbundle_binary_BB(S_end, R_end, level + 1);

			assert(S_end == pre_S_end && R_end == pre_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}
#endif

		assert(min_deg < R_end - s);

		if (SR_rid[u] < S_end)
		{
			u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end, u);
		}
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
		assert(degree[u] + s > best_solution_size && degree[u] + s > S_end);

		bool can_add_u = false;
		if (S_end + 1 <= s)
		{
			can_add_u = true;
		}
		else
		{
			swap_pos(S_end, SR_rid[u]);
			can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
		}

		if (can_add_u)
		{

			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			if (move_u_to_S_with_prune(u, S_end, R_end, level))
				splex_BB_search(S_end, R_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;

			assert(S_end == t_old_S_end + 1 && SR[S_end - 1] == u);
			restore_SR(S_end, R_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
			for (ui i = 0; i < n; i++)
				assert(!vis[i]);
#endif

			assert(Qv.empty());
			bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_with_prune(S_end, R_end, level);
			if (succeed)
				splex_BB_search(S_end, R_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
		else
		{
			ui pre_best_solution_size = best_solution_size;

			assert(Qv.empty());
			bool succeed = remove_u_from_C_with_prune(u, S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_with_prune(S_end, R_end, level);
			if (succeed)
			{
				splex_BB_search(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
	}

	void sbundle_binary_BB(ui S_end, ui R_end, ui level)
	{
		branch_count++;

		if (S_end > best_solution_size && vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end, pstart, pend, edges, s))
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		if (level > 1)
			reorganize_edges(S_end, R_end, level);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
			assert(level_id[u] > level);
		}
#endif

		ui u = n;
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

		assert(SR_rid[u] < R_end && SR_rid[u] >= S_end);

		assert(level_id[u] == n);
		level_id[u] = level;
		--R_end;
		swap_pos(R_end, SR_rid[u]);

		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
		{
			if (SR_rid[edges[i]] < R_end)
			{
				ui w = edges[i];
				--degree[w];
			}
		}
		sbundle_binary_BB(S_end, R_end, level + 1);

		assert(level_id[u] == level);
		assert(SR_rid[u] == R_end);
		++R_end;
		level_id[u] = n;
		degree[u] = 0;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
		{
			ui w = edges[i];
			assert(SR[SR_rid[w]] == w);
			if (SR_rid[w] < R_end)
			{
				++degree[w];
				++degree[u];
			}
		}

		assert(level_id[u] > level);
		if (SR_rid[u] != S_end)
		{
			swap_pos(S_end, SR_rid[u]);
		}
		bool can_add = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
		if (can_add)
		{
			assert(SR_rid[u] == S_end);
			++S_end;
			sbundle_binary_BB(S_end, R_end, level + 1);
			--S_end;
		}
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
				for (ept j = pstart[SR[i]]; j < pend[SR[i]] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 1;
				for (ui j = 0; j < S2_n; j++)
					if (!vis[S2[j]])
						++nn;
				for (ept j = pstart[SR[i]]; j < pend[SR[i]] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 0;
			}
			else
				nn = S_end - degree_in_S[SR[i]];
			if (nn > max_nn)
				max_nn = nn;
			vp[i - S_end].second = nn;
		}
		vector<ui> cnt(max_nn + 1, 0);
		for (ui i = 0; i <= max_nn; i++)
			cnt[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
			++cnt[vp[i].second];
		for (ui i = 0; i < max_nn; i++)
			cnt[i + 1] += cnt[i];
		for (ui i = max_nn; i > 0; i--)
			cnt[i] = cnt[i - 1];
		cnt[0] = 0;
		vector<ui> ids = nonneighbors;
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
			ui idx = ids[i], v = SR[S_end + ids[i]];
			ui t_support = total_support - vp[idx].second;
			for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
				vis[edges[j]] = 1;
			ui j = 0, v_support = s - 1 - S_end + degree_in_S[v], ub = S_end + 1;
			while (true)
			{
				if (j == new_n)
					j = i + 1;
				if (j >= vp.size() || ub > best_solution_size || ub + vp.size() - j <= best_solution_size)
					break;
				ui u = SR[S_end + ids[j]], nn = vp[ids[j]].second;
				if (t_support < nn || (nn != 0 && ub + (t_support / nn) <= best_solution_size))
					break;
				if (vis[u])
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
			for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
				vis[edges[j]] = 0;
			if (ub <= best_solution_size)
			{
				level_id[v] = level;
				Qv.push(v);
			}
			else
				ids[new_n++] = ids[i];
		}
	}

	bool greedily_add_vertices_to_S(ui &S_end, ui &R_end, ui level)
	{
		while (true)
		{
			vector<ui> candidates;
			ui candidates_n = 0;
			for (ui i = S_end; i < R_end; i++)
			{
				ui u = SR[i];
				if (R_end - degree[u] > s)
					continue;

				ui neighbors_n = 0, non_neighbors_n = 0;
				get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, non_neighbors_n);
				assert(non_neighbors_n < s);
				bool OK = true;
				for (ui j = 0; j < non_neighbors_n; j++)
					if (R_end - degree[nonneighbors[j]] > s)
					{
						OK = false;
						break;
					}
				if (OK)
					candidates[candidates_n++] = u;
			}

			if (!candidates_n)
				break;

			while (candidates_n)
			{
				ui u = candidates[--candidates_n];
				assert(SR_rid[u] >= S_end);
				if (SR_rid[u] >= R_end)
					return false;

				bool can_add_u = false;
				if (S_end + 1 <= s)
				{
					can_add_u = true;
				}
				else
				{
					swap_pos(S_end, SR_rid[u]);
					can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
				}
				if (!can_add_u || !move_u_to_S_with_prune(u, S_end, R_end, level))
					return false;
			}
		}
		return true;
		/*该函数的目标是在保证集合 S 仍为 k-plex 的前提下，贪婪地尝试将 R 集合中的顶点添加到 S 集合中。如果某个顶点的添加导致 S 不满足 k-plex 的条件，则该顶点不会被添加到 S 中。
		这个过程会不断重复，直到没有更多合适的顶点可以添加为止。如果在某次尝试中发现顶点无法被添加而不违反 k-plex 的条件，则函数返回 false，否则，如果所有顶点都成功添加，返回 true。*/
	}

	bool greedily_add_nonneighbors(ui *candidates, ui candidates_n, ui &S_end, ui &R_end, ui level)
	{
		while (candidates_n)
		{
			ui u = candidates[--candidates_n];
			assert(SR_rid[u] >= S_end);

			bool can_add_u = false;
			if (S_end + 1 <= s)
			{
				can_add_u = true;
			}
			else
			{
				swap_pos(S_end, SR_rid[u]);
				can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
			}
			if (SR_rid[u] >= R_end || !can_add_u || !move_u_to_S_with_prune(u, S_end, R_end, level))
				return false;
		}
		return true;
	}

	void get_neighbors_and_non_neighbors(ui u, ui idx_start, ui idx_end, ui level, ui &neighbors_n, ui &non_neighbors_n)
	{
		neighbors_n = non_neighbors_n = 0;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < idx_end)
				vis[edges[i]] = 1;
		for (ui i = idx_start; i < idx_end; i++)
			if (SR[i] != u)
			{
				if (vis[SR[i]])
					neighbors[neighbors_n++] = SR[i];
				else
					nonneighbors[non_neighbors_n++] = SR[i];
			}
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			vis[edges[i]] = 0;
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
#endif
		if (SR_rid[u] != S_end)
			swap_pos(S_end, SR_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];

		while (!Qv.empty())
			Qv.pop();

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
				for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 1;
				for (ui j = S_end; j < R_end; j++)
					if (level_id[SR[j]] > level && !vis[SR[j]])
					{
						level_id[SR[j]] = level;
						Qv.push(SR[j]);
					}
				for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 0;
			}
		}
		return remove_vertices_with_prune(S_end, R_end, level);
	}

	bool remove_vertices_with_prune(ui S_end, ui &R_end, ui level)
	{
		while (true)
		{
			while (!Qv.empty())
			{
				ui u = Qv.front();
				Qv.pop();
				assert(SR[SR_rid[u]] == u);
				assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
				--R_end;
				swap_pos(SR_rid[u], R_end);

				bool terminate = false;
				ui neighbors_n = 0;

				for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
					if (SR_rid[edges[i]] < R_end)
					{

						ui w = edges[i];
						neighbors[neighbors_n++] = w;
						--degree[w];
						if (degree[w] + s <= best_solution_size)
						{
							if (SR_rid[w] < S_end)
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
					++R_end;
					return false;
				}
			}

#ifndef NDEBUG

			for (ui i = 0; i < R_end; i++)
			{
				assert(level_id[SR[i]] > level);
			}
#endif

			assert(Qv.empty());

			if (S_end >= 1)
			{
				vector<ui> deleted_vertices;
				delete_vertices_based_on_dist(S_end, R_end, s, deleted_vertices, level);
				if (!deleted_vertices.empty())
				{
					for (ui i = 0; i < deleted_vertices.size(); i++)
					{
						ui t_v = deleted_vertices[i];
						level_id[t_v] = level;
						assert(SR_rid[t_v] >= S_end && SR_rid[t_v] < R_end);
						Qv.push(t_v);
					}
				}
			}
			if (Qv.empty())
			{
				break;
			}
		}

		return true;
	}

	void restore_SR(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level)
	{
		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

		for (; R_end < old_R_end; R_end++)
		{
			ui u = SR[R_end];
			assert(level_id[u] == level && SR_rid[u] == R_end);
			level_id[u] = n;

			degree[u] = degree_in_S[u] = 0;
			for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			{
				ui w = edges[i];
				assert(SR[SR_rid[w]] == w);
				if (SR_rid[w] < R_end)
				{
					++degree[w];
					++degree[u];
				}
				if (SR_rid[w] < S_end)
					++degree_in_S[u];
			}
		}

		for (; S_end > old_S_end; S_end--)
		{
			ui u = SR[S_end - 1];
			assert(SR_rid[u] == S_end - 1);

			for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			{
				ui w = edges[i];
				if (SR_rid[w] < R_end)
					--degree_in_S[w];
			}
		}
	}

	void delete_vertices_based_on_dist(ui S_end, ui R_end, ui s, vector<ui> &results, ui level)
	{

		vector<ui> dist(R_end, n);
		queue<ui> qq;
		assert(S_end > 0 && s >= 2);
		ui tt_u = SR[S_end - 1];
		assert(SR_rid[tt_u] < S_end);
		dist[S_end - 1] = 0;
		qq.push(tt_u);
		while (!qq.empty())
		{
			ui current = qq.front();
			qq.pop();
			for (ui i = pstart[current]; i < pend[current] && level_id[edges[i]] >= level; ++i)
			{
				ui neighbor = edges[i];
				ui neighbor_id = SR_rid[neighbor];
				if (neighbor_id < R_end && dist[neighbor_id] == n)
				{
					dist[neighbor_id] = dist[SR_rid[current]] + 1;
					qq.push(neighbor);
				}
			}
		}

		for (ui i = S_end; i < R_end; i++)
		{
			assert(s >= 2);

			if (dist[i] > 2 + (s - 2) / (best_solution_size + 1 - s))
			{
				results.push_back(SR[i]);
			}
		}
	}

	bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = SR[S_end - 1];
		--S_end;
		--R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool terminate = false;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < R_end)
			{
				ui v = edges[i];
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

		bool terminate = false;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < R_end)
			{
				ui v = edges[i];
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
			}

		return true;
	}

	void swap_pos(ui i, ui j)
	{
		swap(SR[i], SR[j]);
		SR_rid[SR[i]] = i;
		SR_rid[SR[j]] = j;
	}

	void swap_pos2(ui i, ui j)
	{
		swap(SR_remap[i], SR_remap[j]);
		SR_remap_rid[SR_remap[i]] = i;
		SR_remap_rid[SR_remap[j]] = j;
	}

	ui choose_branch_vertex(ui S_end, ui R_end, ui level)
	{
		ui *D = buf;
		ui D_n = 0;
		for (ui i = 0; i < R_end; i++)
			if (R_end - degree[SR[i]] > s)
				D[D_n++] = SR[i];
		assert(D_n != 0);

		ui min_degree_in_S = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] < min_degree_in_S)
				min_degree_in_S = degree_in_S[D[i]];
		ui u = n, min_degree = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] == min_degree_in_S && degree[D[i]] < min_degree)
			{
				min_degree = degree[D[i]];
				u = D[i];
			}
		assert(u != n);

		if (SR_rid[u] < S_end)
		{
			ui max_degree = 0, b = n;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(u, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for (ui i = 0; i < nonneighbors_n; i++)
				if (degree[nonneighbors[i]] > max_degree)
				{
					max_degree = degree[nonneighbors[i]];
					b = nonneighbors[i];
				}
			return b;
		}
		else if (degree_in_S[u] < S_end || R_end - degree[u] > s + 1)
			return u;
		else
		{
			ui max_degree = 0, w = n;
			for (ui i = S_end; i < R_end; i++)
				if (degree[SR[i]] > max_degree)
				{
					max_degree = degree[SR[i]];
					w = SR[i];
				}
			if (degree[w] + 1 >= R_end)
			{
				printf("!!! WA degree[w]: %u, R_end: %u\n", degree[w], R_end);
			}
			assert(degree[w] + 1 < R_end);
			if (R_end - degree[w] == 2)
				return w;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(w, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for (ui i = 0; i < nonneighbors_n; i++)
				if (R_end - degree[nonneighbors[i]] == s + 1)
					return nonneighbors[i];
		}

		printf("!!! WA in choose_branch_vertex\n");
		return n;
	}

	ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end, ui vertex_in_S)
	{
		assert(SR_rid[vertex_in_S] < S_end);
		ui u = n, min_degree_in_S = n;
		for (ui i = S_end; i < R_end; i++)
		{
			ui v = SR[i];
			bool is_adj = false;
			for (ui j = pstart[vertex_in_S]; j < pend[vertex_in_S]; j++)
			{
				if (edges[j] == v)
				{
					is_adj = true;
					break;
				}
			}

			if (!is_adj)
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

#ifndef NDEBUG
		bool is_adj = false;
		for (ui j = pstart[vertex_in_S]; j < pend[vertex_in_S]; j++)
		{
			if (edges[j] == u)
			{
				is_adj = true;
				break;
			}
		}
		assert(u != n && !is_adj);
#endif
		return u;
	}

	bool binary_find(ui u, ui w, ept begin, ept end, ui *edges)
	{
		if (begin >= end)
			return 0;

		ui idx;
		while (begin + 1 < end)
		{
			idx = begin + (end - begin) / 2;
			if (edges[idx] > w)
				end = idx;
			else
				begin = idx;
		}

		if (edges[begin] == w)
		{
			return 1;
		}

		return 0;
	}
};
#endif
