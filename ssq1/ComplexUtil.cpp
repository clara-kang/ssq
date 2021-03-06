#include <memory>
#include <vector>
#include "complexUtil.h"
#include "PathFinding.h"
#include <algorithm>
#include <map>
#include <cmath>
#include <queue>
#include <Eigen/Dense>
#include <Eigen/SparseExtra>

//#include <tuple>

using namespace ComplexUtil;
using namespace std;

double PI = 3.1415926;


std::map<std::pair<int, int>, trnsfr_func_t> line_match_to_func = {
{{0, 2}, {false, {1, 1}, {0, 1}}},
{{2, 0}, {false, {1, 1}, {0, -1}}},
{{0, 3}, {true, {-1, 1}, {0, 0}}},
{{3, 0}, {true, {1, -1}, {0, 0}}},
{{0, 0}, {false, {-1, -1}, {1, 0}}},
{{0, 1}, {true, {1, -1}, {1, 1}}},
{{1, 0}, {true, {-1, 1}, {1, -1}}},
{{1, 3}, {false, {1, 1}, {-1, 0}}},
{{3, 1}, {false, {1, 1}, {1, 0}}},
{{1, 1}, {false, {-1, -1}, {2, 1}}},
{{1, 2}, {true, {1, -1}, {0, 2}}},
{{2, 1}, {true, {-1, 1}, {2, 0}}},
{{2, 2}, {false, {-1, -1}, {1, 2}}},
{{2, 3}, {true, {1, -1}, {-1, 1}}},
{{3, 2}, {true, {-1, 1}, {1, 1}}},
{{3, 3}, {false, {-1, -1}, {0, 1}}}
};

void ComplexUtil::findVertexTypes(HalfEdge::vert_t *saddles, HalfEdge::vert_t *maxs,
	HalfEdge::vert_t *mins, vert_type_t *vert_types,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U) {

	//find type of each vertex
	auto done_v_types = vector<VERT_TYPE>(HE_verts->size());
	auto done_saddles = vector<int>();
	auto done_maxs = vector<int>();
	auto done_mins = vector<int>();

	for (int v_idx = 0; v_idx < HE_verts->size(); v_idx++) {
		// get vertex value
		double v_val = U(v_idx);

		// get its neighbors
		auto neighbors = HalfEdge::getNeighbors(v_idx, HE_verts, HE_edges);

		//// get the values of all neighbors
		//std::vector<double> values(neighbors->size());
		//for (int nb_idx = 1; nb_idx < neighbors->size(); nb_idx++) {
		// values[nb_idx] = U(nb_idx);
		//}

		// get signs for all neighbors, +1 for >, 0 for =, -1 for <
		std::vector<int> signs(neighbors->size());
		for (int nb_idx = 0; nb_idx < neighbors->size(); nb_idx++) {
			double nb_val = U(neighbors->at(nb_idx));
			if (nb_val > v_val) {
				signs[nb_idx] = 1;
			}
			else if (nb_val == v_val) {
				signs[nb_idx] = 0;
			}
			else {
				signs[nb_idx] = -1;
			}
		}

		// we want all sign=0 to equal to previous sign that is not zero
		int nb_idx = 0;
		// find first non_zero sign
		while (signs[nb_idx] == 0) {
			nb_idx++;
		}
		nb_idx++;
		// make all sign=0 to equal to previous sign that is not zero
		for (int i = 0; i < neighbors->size(); i++) {
			if (signs[nb_idx % neighbors->size()] == 0) {
				signs[nb_idx % neighbors->size()] = signs[(nb_idx - 1) % neighbors->size()];
			}
			nb_idx++;
		}

		// get sign changes
		int sign_change_times = 0;
		for (int nb_idx = 0; nb_idx < neighbors->size(); nb_idx++) {
			if (signs[nb_idx] != signs[(nb_idx + neighbors->size() - 1) % neighbors->size()]) {
				sign_change_times++;
			}
		}

		//bool first_sign = (U(neighbors->at(0)) >= v_val);
		//bool last_sign = first_sign; // false for less than
		//bool sign = false;
		//// iterate over neighbor values
		//for (int nb_idx = 1; nb_idx < neighbors->size(); nb_idx++) {
		// double nb_val = U(neighbors->at(nb_idx));
		// sign = (nb_val > v_val);
		// if (sign != last_sign) {
		// sign_change_times++;
		// }
		// last_sign = sign;
		//}
		//// last one
		//if (last_sign != first_sign) {
		// sign_change_times++;
		//}

		if (sign_change_times == 0) {
			if (signs[0] > 0) {
				done_v_types.at(v_idx) = VERT_TYPE::MIN;
				done_mins.push_back(v_idx);
			}
			else {
				done_v_types.at(v_idx) = VERT_TYPE::MAX;
				done_maxs.push_back(v_idx);
			}
		}
		else if (sign_change_times == 2) {
			done_v_types.at(v_idx) = VERT_TYPE::REG;
		}
		else {
			//cout << "sign_change_times: " << sign_change_times << endl;
			done_v_types.at(v_idx) = VERT_TYPE::SADDLE;
			done_saddles.push_back(v_idx);
		}
	}
	vert_type_t done_v_types_ptr = std::make_shared<vector<VERT_TYPE>>(done_v_types);
	HalfEdge::vert_t done_saddles_ptr = std::make_shared<vector<int>>(done_saddles);
	HalfEdge::vert_t done_maxs_ptr = std::make_shared<vector<int>>(done_maxs);
	HalfEdge::vert_t done_mins_ptr = std::make_shared<vector<int>>(done_mins);
	swap(done_v_types_ptr, *vert_types);
	swap(done_saddles_ptr, *saddles);
	swap(done_maxs_ptr, *maxs);
	swap(done_mins_ptr, *mins);
}

int findNext(HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U,
	int v_indx, bool asc) {
	// get neighbors of vertex
	std::shared_ptr<vector<int>> neighbors = HalfEdge::getNeighbors(v_indx, HE_verts, HE_edges);

	// print out neighbor U values
	//for (auto it = neighbors->begin(); it != neighbors->end(); ++it) {
	// cout << "nb " << *it << ": " << U(*it) << endl;
	//}

	if (asc) {
		// find the max/min among neighbors
		auto comp = [&U](const int &lhs, const int &rhs) { if (U(lhs) == U(rhs)) { return false; } return U(lhs) < U(rhs); };
		auto max_nb_it = std::max_element(neighbors->begin(), neighbors->end(), comp);
		//cout << "next elem is " << U(*max_nb_it) << endl;
		return *max_nb_it;
	}
	else {
		auto comp = [&U](const int &lhs, const int &rhs) {
			/*cout << "U(lhs): " << U(lhs) << endl;
			cout << "U(rhs): " << U(rhs) << endl;*/
			if (U(lhs) == U(rhs)) {
				return false;
			}
			return U(lhs) < U(rhs); };
		auto min_nb_it = std::min_element(neighbors->begin(), neighbors->end(), comp);
		//cout << "next elem is " << U(*min_nb_it) << endl;
		return *min_nb_it;
	}

}

steep_lines_t ComplexUtil::findSteepLines(HalfEdge::vert_t saddles, vert_type_t vert_types,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U) {

	// map to record all steep lines, start vertex and end vertex used as keys
	std::map<std::pair<int, int>, std::vector<int>> steeplines;

	int sl_per_sdl = 0;
	for (auto it = saddles->begin(); it != saddles->end(); ++it) {

		sl_per_sdl = 0;

		// find steep lines going out from each saddle
		// get all neighbors of saddle
		HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(*it, HE_verts, HE_edges);

		// print U values for neighbors
		//for (int nb_id = 0; nb_id < neighbors->size(); nb_id++) {
		// cout << "nb: " << neighbors->at(nb_id) << ", U: " << U(neighbors->at(nb_id)) << endl;
		//}
		for (int nb_id = 0; nb_id < neighbors->size(); nb_id++) {
			int nb = neighbors->at(nb_id);
			int next_nb = neighbors->at((nb_id + 1) % neighbors->size());
			int prev_nb = neighbors->at((nb_id - 1 + neighbors->size()) % neighbors->size());

			bool ascending = true;
			// check if it's a local max
			if (U(nb) >= U(*it) && U(nb) >= U(prev_nb) && U(nb) >= U(next_nb)) {
				ascending = true;
				sl_per_sdl++;
			}
			// check if it's a local min
			else if (U(nb) < U(*it) && U(nb) < U(prev_nb) && U(nb) < U(next_nb)) {
				ascending = false;
				sl_per_sdl++;
			}
			// not a local min/max
			else {
				continue;
			}

			// find the steep line
			vector<int> steep_line;
			int node_in_line = nb;
			steep_line.push_back(*it);
			steep_line.push_back(node_in_line);
			// while not reaching max/min, keep getting next one
			VERT_TYPE extrm_type = ascending ? VERT_TYPE::MAX : VERT_TYPE::MIN;
			while (vert_types->at(node_in_line) != extrm_type) {
				node_in_line = findNext(HE_verts, HE_edges, U, node_in_line, ascending);
				steep_line.push_back(node_in_line);
			}
			// record steep line
			steeplines.insert(steeplines.begin(), { std::make_pair(steep_line[0], steep_line.back()), steep_line });
		}
		cout << "sl_per_sdl: " << sl_per_sdl << endl;
	}
	steep_lines_t steep_lines_ptr = std::make_shared<std::map<std::pair<int, int>, std::vector<int>>>(steeplines);
	return steep_lines_ptr;
}

std::shared_ptr<std::vector<std::vector<int>>> ComplexUtil::buildMsPatches(steep_lines_t steeplines, std::shared_ptr<Eigen::MatrixXd> vertices_ptr,
	std::shared_ptr<Eigen::MatrixXd> vns_ptr) {

	// the projections of start to end vectors of sls, perpendicular to vertex normals 
	std::map<std::pair<int, int>, Eigen::Vector3d> s2e_angl;

	// the half steeplines going out from vertex
	std::map<int, std::vector<int>> v2_sl;

	// iterate over steeplines
	for (auto it = steeplines->begin(); it != steeplines->end(); ++it) {
		// get vertex position of start and end
		Eigen::Vector3d start = vertices_ptr->row(it->first.first);
		Eigen::Vector3d end = vertices_ptr->row(it->first.second);
		// get vertex normal of start and end
		Eigen::Vector3d start_norm = vns_ptr->row(it->first.first);
		Eigen::Vector3d end_norm = vns_ptr->row(it->first.second);

		// vector from start to end, projected perpendicular to start normal
		Eigen::Vector3d s2e = end - start;
		Eigen::Vector3d s2e_prpndclr = s2e - s2e.dot(start_norm) * start_norm;

		// vector from end to start, projected perpendicular to end normal
		Eigen::Vector3d e2s = start - end;
		Eigen::Vector3d e2s_prpndclr = e2s - e2s.dot(end_norm) * end_norm;

		s2e_angl.insert(s2e_angl.begin(), { it->first, s2e_prpndclr });
		s2e_angl.insert(s2e_angl.begin(),
			{ std::make_pair(it->first.second, it->first.first), e2s_prpndclr });
	}

	// build v2_sl, the adjacency list
	for (auto it = s2e_angl.begin(); it != s2e_angl.end(); ++it) {
		if (v2_sl.find(it->first.first) == v2_sl.end()) {
			v2_sl.insert(v2_sl.begin(), { it->first.first, std::vector<int>({it->first.second}) });
		}
		else {
			v2_sl[it->first.first].push_back(it->first.second);
		}
	}

	// patches
	std::vector<std::vector<int>> ms_patches;

	// filter? indicating whether the sl belongs to a patch
	std::map<std::pair<int, int>, bool> sl_cover;
	for (auto it = s2e_angl.begin(); it != s2e_angl.end(); ++it) {
		sl_cover.insert(sl_cover.begin(), { it->first, false });
	}

	for (auto v2_sl_it = v2_sl.begin(); v2_sl_it != v2_sl.end(); ++v2_sl_it) {
		for (auto ends_it = v2_sl_it->second.begin(); ends_it != v2_sl_it->second.end(); ends_it++) {

			// if half sl does not belong to any patch
			if (!sl_cover[{v2_sl_it->first, *ends_it}]) {

				int end_vert = *ends_it;
				int start_vert = v2_sl_it->first;
				std::vector<int> patch({ start_vert, end_vert });

				for (int l_num = 0; l_num < 3; l_num++) {

					// mark sl as dealt with
					sl_cover[{start_vert, end_vert}] = true;

					// look for the sl that makes the smallest counter clockwise angle with current sl
					std::map<int, double> end_angles;
					std::vector<int> &end_nbs = v2_sl[end_vert];
					for (auto end_nb_it = end_nbs.begin(); end_nb_it != end_nbs.end(); ++end_nb_it) {
						if (*end_nb_it == start_vert) {
							continue;
						}
						// get ccw angle
						else {
							Eigen::Vector3d v1 = s2e_angl[{end_vert, start_vert}];
							Eigen::Vector3d v2 = s2e_angl[{end_vert, *end_nb_it}];
							Eigen::Vector3d end_norm = vns_ptr->row(end_vert);

							// get angle
							double angle = acos(v2.dot(v1) / (v1.norm() * v2.norm()));
							// if angle should be > 180
							if (end_norm.dot(v2.cross(v1)) < 0) {
								angle = 2.0 * PI - angle;
							}
							// diff between angle and 90 degrees
							double angle_diff = abs(PI / 2.0 - angle);
							end_angles.insert(end_angles.begin(), { *end_nb_it, angle });
						}
					}
					// find the nb with smallest end_angle
					double min_angle = INFINITY;
					int min_nb;
					for (auto end_angles_it = end_angles.begin(); end_angles_it != end_angles.end(); ++end_angles_it) {
						if (end_angles_it->second < min_angle) {
							min_angle = end_angles_it->second;
							min_nb = end_angles_it->first;
						}
					}

					// min_nb is the next vertex in the patch
					start_vert = end_vert;
					end_vert = min_nb;
					patch.push_back(min_nb);
				}
				sl_cover[{start_vert, end_vert}] = true;
				ms_patches.push_back(patch);
			}
		}
	}
	std::shared_ptr<vector<vector<int>>> ms_patch_ptr = std::make_shared<vector<vector<int>>>(ms_patches);
	return ms_patch_ptr;
}

vector<int> &find_sl(steep_lines_t steeplines, int start_vert, int end_vert, bool *reversed) {
	auto sl_it = steeplines->find({ start_vert, end_vert });
	if (sl_it != steeplines->end()) {
		vector<int> &sl = steeplines->at({ start_vert, end_vert });
		*reversed = false;
		return sl;
	}
	else {
		vector<int> &sl = steeplines->at({ end_vert, start_vert });
		*reversed = true;
		return sl;
	}
}

// given a specific patch, return the index (the ith) of the steepline in the patch
int find_sl_id(int cur_vert, const steep_lines_t steeplines,
	const vector<int> &patch_nodes) {
	bool reversed;
	for (int i = 0; i < 4; i++) {
		vector<int> &sl = find_sl(steeplines, patch_nodes[i], patch_nodes[(i + 1) % 4], &reversed);
		if (std::find(sl.begin(), sl.end(), cur_vert) != sl.end()) {
			return i;
		}
	}
	return -1;
}

void graphSearchRec(const steep_lines_t steeplines, int start_vert,
	HalfEdge::vert_t HE_vert, HalfEdge::edge_t HE_edges,
	const vector<int> &patch_nodes, vector<int> &patch_verts, vector<bool> &visited) {

	std::queue<int> to_visit({ start_vert });

	while (!to_visit.empty()) {
		int cur_vert = to_visit.front();
		to_visit.pop();

		if (!visited[cur_vert]) {
			visited[cur_vert] = true;
			patch_verts.push_back(cur_vert);
			int cur_line_id = find_sl_id(cur_vert, steeplines, patch_nodes);
			// if reach corner, stop
			if (cur_line_id >= 0 && std::find(patch_nodes.begin(), patch_nodes.end(), cur_vert) != patch_nodes.end()) {
				continue;
			}
			// only allow to progress along boundary
			else if (cur_line_id >= 0) {
				bool rev;
				// push next and prec on sl
				vector<int> &sl = find_sl(steeplines, patch_nodes[cur_line_id],
					patch_nodes[(cur_line_id + 1) % 4], &rev);
				auto cur_id_in_sl = std::find(sl.begin(), sl.end(), cur_vert);
				to_visit.push(*(cur_id_in_sl - 1));
				to_visit.push(*(cur_id_in_sl + 1));
			}
			else {
				// get neighbors
				HalfEdge::vert_t nbs = HalfEdge::getNeighbors(cur_vert, HE_vert, HE_edges);

				for (auto nb_it = nbs->begin(); nb_it != nbs->end(); ++nb_it) {
					to_visit.push(*nb_it);
				}
			}
		}
	}
}

void ComplexUtil::fillMsPatches(steep_lines_t steeplines, patch_t ms_patches,
	HalfEdge::vert_t HE_vert, HalfEdge::edge_t HE_edges,
	std::shared_ptr<Eigen::MatrixXd> vertices_ptr,
	std::shared_ptr<Eigen::MatrixXd> vns_ptr,
	std::shared_ptr<vector<vector<int>>> *patch_verts) {

	vector<vector<int>> patches_vertices(ms_patches->size());
	int cnt = 0;
	for (auto patch_it = ms_patches->begin(); patch_it != ms_patches->end(); patch_it++) {
		// get vertex on the left side of the first steepline as start
		bool reversed; // whether the sl should be reversed
		vector<int> &first_sl = find_sl(steeplines, (*patch_it)[0], (*patch_it)[1], &reversed);
		int mid_vert = first_sl[first_sl.size() / 2];
		// next vert after mid vert
		int next_vert;
		if (reversed) {
			next_vert = first_sl[first_sl.size() / 2 - 1];
		}
		else {
			next_vert = first_sl[first_sl.size() / 2 + 1];
		}

		// get neighbors of mid vertex
		HalfEdge::vert_t mid_vert_nbs = HalfEdge::getNeighbors(mid_vert, HE_vert, HE_edges);

		// find index of next_vert in mid_vert_nbs
		auto next_v_it = std::find(mid_vert_nbs->begin(), mid_vert_nbs->end(), next_vert);
		int next_v_indx = std::distance(mid_vert_nbs->begin(), next_v_it);
		// get next vertex ccw, it should be a vertex in the patch
		int inside_vert = mid_vert_nbs->at((next_v_indx - 1 + mid_vert_nbs->size()) % mid_vert_nbs->size());

		vector<int> patch_verts;
		vector<bool> visited(vertices_ptr->rows(), false);
		cout << "inside_vert: " << inside_vert << endl;
		graphSearchRec(steeplines, inside_vert, HE_vert, HE_edges,
			*patch_it, patch_verts, visited);
		patches_vertices[cnt] = patch_verts;
		cnt++;
	}
	shared_ptr<vector<vector<int>>> patches_vertices_ptr = std::make_shared<vector<vector<int>>>(patches_vertices);
	swap(*patch_verts, patches_vertices_ptr);
}

bool ComplexUtil::complexValidityCheck(patch_t ms_patches) {
	for (auto patch_it = ms_patches->begin(); patch_it != ms_patches->end(); ++patch_it) {
		if (patch_it->at(4) != patch_it->at(0)) {
			return false;
		}
	}
	return true;
}

bool ComplexUtil::patchValidityCheck(std::shared_ptr<vector<vector<int>>> patch_verts,
	int N, const vector<int> &verts_on_sls) {
	vector<int> vert_patch_ids(N, -1);
	int patch_cnt = 0;
	bool valid = true;
	int dup = 0;
	for (auto patch_it = patch_verts->begin(); patch_it != patch_verts->end(); ++patch_it) {
		for (auto vert_it = patch_it->begin(); vert_it != patch_it->end(); ++vert_it) {
			// if vertex is not on any steepline, it can only appear in one patch
			if (std::find(verts_on_sls.begin(), verts_on_sls.end(), *vert_it) == verts_on_sls.end()) {
				// it already belongs to a patch
				if (vert_patch_ids[*vert_it] >= 0) {
					valid = false;
					cout << "collision patch " << vert_patch_ids[*vert_it] << ", and patch " << patch_cnt << endl;
					dup++;
				}
			}
			vert_patch_ids[*vert_it] = patch_cnt;
		}
		patch_cnt++;
	}

	return valid;
}

std::pair<int, int> getSlKey(steep_lines_t steeplines, std::pair<int, int> sl_id) {
	auto sl_it = steeplines->find(sl_id);
	if (sl_it != steeplines->end()) {
		return sl_id;
	}
	else {
		return { sl_id.second, sl_id.first };
	}
}

steep_lines_t ComplexUtil::buildSLtoPatchMap(steep_lines_t steeplines, patch_t ms_patches) {
	std::map<std::pair<int, int>, std::vector<int>> sl_patch_map;
	for (auto sl_it = steeplines->begin(); sl_it != steeplines->end(); ++sl_it) {
		sl_patch_map.insert(sl_patch_map.begin(), { sl_it->first, std::vector<int>({}) });
	}

	for (int patch_cnt = 0; patch_cnt < ms_patches->size(); patch_cnt++) {
		for (int i = 0; i < 4; i++) {
			int s = ms_patches->at(patch_cnt)[i];
			int e = ms_patches->at(patch_cnt)[(i + 1) % 4];
			sl_patch_map[getSlKey(steeplines, { s, e })].push_back(patch_cnt);
		}
	}
	std::shared_ptr <std::map<std::pair<int, int>, std::vector<int>>> ptr
		= std::make_shared<std::map<std::pair<int, int>, std::vector<int>>>(sl_patch_map);
	steep_lines_t sl_patch_map_ptr = ptr;
	return sl_patch_map_ptr;
}

int which_line(std::vector<int> &patch, int start, int end) {
	for (int i = 0; i < 4; i++) {
		if (patch[i] == start && patch[(i + 1) % 4 == end]) {
			return i;
		}
	}
	return -1;
}

int ComplexUtil::which_node(std::vector<int> &patch, int vert) {
	auto it = std::find(patch.begin(), patch.end(), vert);
	return std::distance(patch.begin(), it);
}

int getNbPatchId(steep_lines_t steep_lines, steep_lines_t sl_patch_map, int self_id, int start, int end) {
	vector<int> &ids = sl_patch_map->at(getSlKey(steep_lines, { start, end }));
	if (ids[0] == self_id) {
		return ids[1];
	}
	return ids[0];
}

// patch corner 0,1,2,3 corresponds to (0,0), (1,0), (1,1), (0,1)
trnsfr_funcs_map_t ComplexUtil::buildTrnsfrFuncsAndPatchGraph(steep_lines_t steeplines,
	steep_lines_t sl_patch_map, patch_t ms_patches, patch_t *patch_graph_in) {
	std::map<std::pair<int, int>, trnsfr_func_t> trnsfr_funcs_map;
	std::vector<std::vector<int>> patch_graph(ms_patches->size(), std::vector<int>({ 0, 0, 0, 0 }));
	int patch_cnt = 0;
	for (auto patch_it = ms_patches->begin(); patch_it != ms_patches->end(); patch_it++) {
		for (int line_id = 0; line_id < 4; line_id++) {
			int nb_start = (*patch_it)[(line_id + 1) % 4];
			int nb_end = (*patch_it)[line_id];
			int nb_patch_id = getNbPatchId(steeplines, sl_patch_map, patch_cnt, nb_start, nb_end);
			int nb_line_id = which_line(ms_patches->at(nb_patch_id), nb_start, nb_end);
			trnsfr_funcs_map.insert(trnsfr_funcs_map.begin(),
				{ { patch_cnt, nb_patch_id }, line_match_to_func[{line_id, nb_line_id}] });
			patch_graph[patch_cnt][line_id] = nb_patch_id;
		}
		patch_cnt++;
	}
	std::shared_ptr<std::map<std::pair<int, int>, trnsfr_func_t>> ptr =
		std::make_shared<std::map<std::pair<int, int>, trnsfr_func_t>>(trnsfr_funcs_map);
	trnsfr_funcs_map_t trnsfr_funcs_map_ptr = ptr;
	std::shared_ptr<std::vector<std::vector<int>>> ptr2 = std::make_shared<std::vector<std::vector<int>>>(patch_graph);
	patch_t patch_graph_ptr = ptr2;
	std::swap(*patch_graph_in, patch_graph_ptr);
	return trnsfr_funcs_map_ptr;
}

// return f2(f1)
trnsfr_func_t composeTrnsfrFunc(trnsfr_func_t &f1, trnsfr_func_t &f2) {
	trnsfr_func_t composed_f;
	std::get<0>(composed_f) = std::get<0>(f1) != std::get<0>(f2);
	if (std::get<0>(f2)) {
		std::get<1>(composed_f) = { std::get<1>(f1).second * std::get<1>(f2).first, std::get<1>(f1).first * std::get<1>(f2).second };
		std::get<2>(composed_f) = { std::get<2>(f1).second * std::get<1>(f2).first, std::get<2>(f1).first * std::get<1>(f2).second };
	}
	else {
		std::get<1>(composed_f) = { std::get<1>(f1).first * std::get<1>(f2).first, std::get<1>(f1).second * std::get<1>(f2).second };
		std::get<2>(composed_f) = { std::get<2>(f1).first * std::get<1>(f2).first, std::get<2>(f1).second * std::get<1>(f2).second };
	}
	std::get<2>(composed_f).first += std::get<2>(f2).first;
	std::get<2>(composed_f).second += std::get<2>(f2).second;
	return composed_f;
}

trnsfr_func_t inverseTrnsfrFunc(trnsfr_func_t &f) {
	trnsfr_func_t inv_func;
	std::get<0>(inv_func) = std::get<0>(f);
	std::pair<int, int> constants;

	if (std::get<0>(inv_func)) {
		std::get<1>(inv_func) = { std::get<1>(f).second, std::get<1>(f).first };
		constants = { std::get<2>(f).second, std::get<2>(f).first };

	}
	else {
		std::get<1>(inv_func) = std::get<1>(f);
		constants = std::get<2>(f);
	}
	constants.first *= -std::get<1>(inv_func).first;
	constants.second *= -std::get<1>(inv_func).second;
	std::get<2>(inv_func) = constants;
	return inv_func;
}

bool shareNodes(vector<int> &nodes1, vector<int> &nodes2) {
	for (int i = 0; i < 4; i++) {
		if (std::find(nodes2.begin(), nodes2.end(), nodes1[i]) != nodes2.end()) {
			return true;
		}
	}
	return false;
}

void addTrnsfrFunc(int patch_from, int patch_to, patch_t patch_graph, trnsfr_funcs_map_t trnsfr_funcs_map) {
	std::shared_ptr<std::vector<int>> path_to_nb = PathFinding::getPath(patch_from, patch_to, *patch_graph);
	trnsfr_func_t cmpsd_func = trnsfr_funcs_map->at({ patch_from, path_to_nb->at(1) });
	// compose trnsfr function along path
	for (int i = 1; i < path_to_nb->size() - 1; i++) {
		trnsfr_func_t next_func = trnsfr_funcs_map->at({ path_to_nb->at(i),path_to_nb->at(i + 1) });
		cmpsd_func = composeTrnsfrFunc(cmpsd_func, next_func);
	}
	// update trnsfr_funcs_map
	trnsfr_funcs_map->insert({ { patch_from, patch_to }, cmpsd_func });
	trnsfr_funcs_map->insert({ { patch_to, patch_from }, inverseTrnsfrFunc(cmpsd_func) });
	// update patch_graph
	patch_graph->at(patch_from).push_back(patch_to);
	patch_graph->at(patch_to).push_back(patch_from);
}

void ComplexUtil::buildExtendedTrnsfrfuncs(patch_t ms_patches, patch_t patch_graph, trnsfr_funcs_map_t trnsfr_funcs_map) {
	for (int patch_id = 0; patch_id < ms_patches->size(); patch_id++) {
		vector<int> &nb_patch_ids = patch_graph->at(patch_id);
		vector<int> extended_nbs;
		for (auto nb_it = nb_patch_ids.begin(); nb_it != nb_patch_ids.end(); ++nb_it) {
			// get nb's nbs
			vector<int> &nb_nb_patch_ids = patch_graph->at(*nb_it);
			for (auto nb_nb_it = nb_nb_patch_ids.begin(); nb_nb_it != nb_nb_patch_ids.end(); ++nb_nb_it) {
				// add to extended nb sets if not already in, if not already a nb, if share common vertices
				if (*nb_nb_it != patch_id &&
					std::find(nb_patch_ids.begin(), nb_patch_ids.end(), *nb_nb_it) == nb_patch_ids.end() &&
					std::find(extended_nbs.begin(), extended_nbs.end(), *nb_nb_it) == extended_nbs.end() &&
					shareNodes(ms_patches->at(patch_id), ms_patches->at(*nb_nb_it))) {

					extended_nbs.push_back(*nb_nb_it);
				}
			}
		}
		for (auto ex_nb_it = extended_nbs.begin(); ex_nb_it != extended_nbs.end(); ++ex_nb_it) {
			addTrnsfrFunc(patch_id, *ex_nb_it, patch_graph, trnsfr_funcs_map);
		}
	}
}

bool ComplexUtil::isPatchNode(patch_t ms_patches, int vert_idx, int patch_id) {
	vector<int> &patch_nodes = ms_patches->at(patch_id);
	auto it = std::find(patch_nodes.begin(), patch_nodes.end(), vert_idx);
	if (it != patch_nodes.end()) {
		return true;
	}
	return false;
}

bool isNode(vector<int> &patch_nodes_v, int v_idx) {
	return std::find(patch_nodes_v.begin(), patch_nodes_v.end(), v_idx) != patch_nodes_v.end();
}

void ComplexUtil::applyTrnsfrFunc(Eigen::Vector2d &initial_uv,
	Eigen::Vector2d &res_uv, trnsfr_func_t trnsfr_func) {

	if (std::get<0>(trnsfr_func)) {
		res_uv << initial_uv(1), initial_uv(0);
	}
	else {
		res_uv << initial_uv(0), initial_uv(1);
	}
	res_uv(0) *= std::get<1>(trnsfr_func).first;
	res_uv(1) *= std::get<1>(trnsfr_func).second;
	res_uv(0) += std::get<2>(trnsfr_func).first;
	res_uv(1) += std::get<2>(trnsfr_func).second;
}

std::shared_ptr<Eigen::MatrixXd> ComplexUtil::solveForCoords(Eigen::SparseMatrix<double> &L, std::vector<int> &vert_patch_ids,
	ComplexUtil::trnsfr_funcs_map_t trnsfr_funcs_map, int node_num, patch_t ms_patches,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {

	int unknown_size = vert_patch_ids.size() - node_num;
	Eigen::SparseMatrix<double> M(unknown_size * 2, unknown_size * 2);
	Eigen::VectorXd knowns(unknown_size * 2);

	for (int i = 0; i < M.rows(); i++) {
		knowns(i) = 0;
	}

	std::vector<int> verts_matrx_map(vert_patch_ids.size(), -1);

	// build relationship between vertex index to matrix row
	int row_indx = 0;
	for (int v_indx = 0; v_indx < vert_patch_ids.size(); v_indx++) {
		//if (!isPatchNode(ms_patches, v_indx, vert_patch_ids[v_indx])) {
		if (vert_patch_ids[v_indx] >= 0) {
			verts_matrx_map[v_indx] = row_indx;
			row_indx++;
		}
	}

	for (int vert_indx = 0; vert_indx < vert_patch_ids.size(); vert_indx++) {
		if (verts_matrx_map[vert_indx] < 0) {
			continue;
		}

		int patch_id = vert_patch_ids[vert_indx];
		int v_mmap_indx_u = verts_matrx_map[vert_indx];
		int v_mmap_indx_v = v_mmap_indx_u + unknown_size;

		M.insert(v_mmap_indx_u, v_mmap_indx_u) = 1;
		M.insert(v_mmap_indx_v, v_mmap_indx_v) = 1;
		// get its neighbors
		auto neighbors = HalfEdge::getNeighbors(vert_indx, HE_verts, HE_edges);

		double node_weight = 0;

		for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {
			node_weight -= L.coeff(vert_indx, *nb_it); // L is negative
		}
		for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {

			double nb_weight = -L.coeff(vert_indx, *nb_it) / node_weight;
			int nb_patch_id = vert_patch_ids[*nb_it];

			if (verts_matrx_map[*nb_it] >= 0) {
				int nb_mmap_indx_x = verts_matrx_map[*nb_it];
				int nb_mmap_indx_y = verts_matrx_map[*nb_it] + unknown_size;

				if (nb_patch_id != patch_id) {
					// try to find transfer function
					if (trnsfr_funcs_map->find({ nb_patch_id, patch_id }) != trnsfr_funcs_map->end()) {
						trnsfr_func_t &trnsfer_func = trnsfr_funcs_map->at({ nb_patch_id, patch_id });
						// x, y inverted
						if (std::get<0>(trnsfer_func)) {
							std::swap(nb_mmap_indx_x, nb_mmap_indx_y);
						}
						knowns(v_mmap_indx_u) += nb_weight * (std::get<2>(trnsfer_func)).first;
						knowns(v_mmap_indx_v) += nb_weight * (std::get<2>(trnsfer_func)).second;
						M.insert(v_mmap_indx_u, nb_mmap_indx_x) = -nb_weight * (std::get<1>(trnsfer_func)).first;
						M.insert(v_mmap_indx_v, nb_mmap_indx_y) = -nb_weight * (std::get<1>(trnsfer_func)).second;
					}
					else {
						cout << "cannot find trnsfr func between " << nb_patch_id << " and " << patch_id << endl;
					}
				}
				else {
					M.insert(v_mmap_indx_u, nb_mmap_indx_x) = -nb_weight;
					M.insert(v_mmap_indx_v, nb_mmap_indx_y) = -nb_weight;
				}
			}
			else {
				int nb_node_indx = which_node(ms_patches->at(patch_id), *nb_it);
				if (nb_node_indx == 1 || nb_node_indx == 2) {
					knowns(v_mmap_indx_u) += nb_weight;
				}
				if (nb_node_indx == 2 || nb_node_indx == 3) {
					knowns(v_mmap_indx_u) += nb_weight;
				}
			}
		}
	}

	//auto M_dense = Eigen::MatrixXd(M);
	//for (int i = 0; i < M_dense.rows(); i++) {
	//	cout << M_dense(i, i) << endl;
	//}

	//std::cout << "knowns: " << knowns << endl;
	M.makeCompressed();

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solverA;
	solverA.compute(M);

	if (solverA.info() != Eigen::Success) {
		cout << "Oh: Very bad" << endl;
	}

	Eigen::VectorXd X = Eigen::VectorXd(unknown_size * 2);
	X = solverA.solve(knowns);
	//for (int indx = 0; indx < X.rows(); indx++) {
	// std::cout << X(indx) << std::endl;
	//}

	Eigen::MatrixXd uv_coords(vert_patch_ids.size(), 2);

	for (int v_indx = 0; v_indx < vert_patch_ids.size(); v_indx++) {
		int v_mmap_indx_u = verts_matrx_map[v_indx];
		if (v_mmap_indx_u >= 0) {
			int v_mmap_indx_v = v_mmap_indx_u + unknown_size;
			uv_coords.row(v_indx) << X(v_mmap_indx_u), X(v_mmap_indx_v);
		}
		else {
			// it's a node
			uv_coords.row(v_indx) << -100, -100;
		}
	}

	std::shared_ptr<Eigen::MatrixXd> uv_coords_ptr = std::make_shared<Eigen::MatrixXd>(uv_coords);
	return uv_coords_ptr;
}

bool inRange(Eigen::Vector2d uv) {
	return (uv(0) >= 0 && uv(0) <= 1 && uv(1) >= 0 && uv(1) <= 1);
}

void swapBtwnPatches(int v_idx, std::vector<int> &vert_patch_ids,
	std::shared_ptr<vector<vector<int>>> patch_verts,
	int patch_from, int patch_to) {
	// remove from current patch
	vector<int> &from_patch_verts = patch_verts->at(patch_from);
	if (v_idx == 336) {
		cout << "336" << endl;
	}
	auto v_it = std::find(from_patch_verts.begin(), from_patch_verts.end(), v_idx);
	from_patch_verts.erase(v_it);
	// add to the other patch
	patch_verts->at(patch_to).push_back(v_idx);
	// update vert_patch_ids
	vert_patch_ids[v_idx] = patch_to;
}

bool inRangeInPath(std::shared_ptr<Eigen::MatrixXd> uvs, Eigen::Vector2d &uv, int cur_patch_id, int target_patch_id,
	patch_t patch_graph, trnsfr_funcs_map_t trnsfr_functions) {
	Eigen::Vector2d uv_in_other;
	if (trnsfr_functions->find({ cur_patch_id, target_patch_id }) == trnsfr_functions->end()) {
		cout << cur_patch_id << ", " << target_patch_id << endl;
		// transfer function not constructed, append it now, should be rare
		addTrnsfrFunc(cur_patch_id, target_patch_id, patch_graph, trnsfr_functions);

	}
	trnsfr_func_t cur_to_other = trnsfr_functions->at({ cur_patch_id, target_patch_id });
	applyTrnsfrFunc(uv, uv_in_other, cur_to_other);
	return inRange(uv_in_other);
}

void addToMap(map<int, vector<int>> &m, int key, int val) {
	if (m.find(key) != m.end()) {
		m[key].push_back(val);
	}
	else {
		m[key] = std::vector<int>(1, val);
	}
}

void ComplexUtil::adjustBndrys(std::vector<int> &vert_patch_ids, steep_lines_t steeplines,
	patch_t patch_graph, trnsfr_funcs_map_t trnsfr_functions,
	std::shared_ptr<Eigen::MatrixXd> uvs, std::shared_ptr<vector<vector<int>>> patch_verts,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {

	// checked map, true if vertex is in range in current patch
	std::vector<bool> checked(vert_patch_ids.size(), false);

	for (auto sl_it = steeplines->begin(); sl_it != steeplines->end(); ++sl_it) {
		vector<int> &sl_verts = sl_it->second;
		map<int, vector<int>> fronts;

		for (auto sl_v = sl_verts.begin(); sl_v != sl_verts.end(); ++sl_v) {
			// the vertex might not be in its original patch anymore
			if (checked[*sl_v]) {
				continue;
			}

			int cur_patch_id = vert_patch_ids[*sl_v];
			// is node
			if (cur_patch_id < 0) {
				continue;
			}
			Eigen::Vector2d sl_v_uv = uvs->row(*sl_v);

			// not in range in current patch
			if (!inRange(sl_v_uv)) {
				// get the potential patches that sl_v might be in
				vector<int> &cur_patch_nbs = patch_graph->at(cur_patch_id);
				for (auto patch_nb_it = cur_patch_nbs.begin(); patch_nb_it != cur_patch_nbs.end(); ++patch_nb_it) {
					if (inRangeInPath(uvs, sl_v_uv, cur_patch_id, *patch_nb_it, patch_graph, trnsfr_functions)) {
						// swap vert on sl to other patch
						swapBtwnPatches(*sl_v, vert_patch_ids, patch_verts, cur_patch_id, *patch_nb_it);
						// add to other patches's front
						addToMap(fronts, *patch_nb_it, *sl_v);
						break;
					}
				}
			}
			// in range in current patch
			else {
				addToMap(fronts, cur_patch_id, *sl_v);
			}
			checked[*sl_v] = true;
		}

		// advance the fronts
		for (auto fronts_it = fronts.begin(); fronts_it != fronts.end(); ++fronts_it) {
			std::vector<int> &front = fronts_it->second;
			while (front.size() > 0) {
				// get first element
				int v_idx = front[0];
				front.erase(front.begin());
				int front_patch_id = fronts_it->first;

				HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(v_idx, HE_verts, HE_edges);

				for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {
					if (checked[*nb_it]) {
						continue;
					}
					int nb_patch_id = vert_patch_ids[*nb_it];
					// nb is node
					if (nb_patch_id < 0) {
						continue;
					}
					// nb is not in the same patch as front, see if can be assimilated
					if (nb_patch_id != front_patch_id) {
						Eigen::Vector2d cur_uv = uvs->row(*nb_it);
						if (!inRange(cur_uv) && inRangeInPath(uvs, cur_uv, nb_patch_id, front_patch_id, patch_graph, trnsfr_functions)) {
							swapBtwnPatches(*nb_it, vert_patch_ids, patch_verts, nb_patch_id, front_patch_id);
							front.push_back(*nb_it);
							checked[*nb_it] = true;
						}
					}
				}
			}
		}
	}
}



int findVertOnBoundaryRec(std::vector<int> &vert_patch_ids, int start_v, std::vector<int> &visited,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {
	if (std::find(visited.begin(), visited.end(), start_v) != visited.end()) {
		return -1;
	}
	int patch_id = vert_patch_ids[start_v];
	HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(start_v, HE_verts, HE_edges); 
	for (int nb_idx = 0; nb_idx < neighbors->size(); nb_idx++) {
		// found vertex on boundary
		if (vert_patch_ids[neighbors->at(nb_idx)] != patch_id) {
			return start_v;
		}
	}
	// did not find vertex on boundary
	visited.push_back(start_v);
	for (int nb_idx = 0; nb_idx < neighbors->size(); nb_idx++) {
		int v_on_bndry = findVertOnBoundaryRec(vert_patch_ids, neighbors->at(nb_idx),
			visited, HE_verts, HE_edges);
		if (v_on_bndry != -1) {
			return v_on_bndry;
		}
	}
	return -1;
}

int findClosestBoundary(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, 
	int node_v, vector<int> &bndry_verts,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {
	vector<double> dist2node(bndry_verts.size());
	for (int i = 0; i < bndry_verts.size(); i++) {
		dist2node[i] = (vertices_ptr->row(bndry_verts[i]) - vertices_ptr->row(node_v)).norm();
	}
	auto min_bndry_v = std::min_element(dist2node.begin(), dist2node.end());
	return bndry_verts[std::distance(dist2node.begin(), min_bndry_v)];
}

void updatePatchNode(patch_t ms_patches, int old_node, int new_node, 
	std::vector<int> &vert_patch_ids, int assign_patch_id) {
	int any_patch_id = -1;
	for (int patch_id = 0; patch_id < ms_patches->size(); patch_id++) {
		for (int i = 0; i < 4; i++) {
			if (ms_patches->at(patch_id)[i] == old_node) {
				any_patch_id = patch_id;
				ms_patches->at(patch_id)[i] = new_node;
			}
		}
	}
	if (assign_patch_id >= 0) {
		vert_patch_ids[old_node] = assign_patch_id;
	}
	else {
		vert_patch_ids[old_node] = any_patch_id;
	}
}

bool onBoundary(patch_t ms_patches, int vert, int patch_id, std::vector<int> &vert_patch_ids,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {
	HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(vert, HE_verts, HE_edges);
	std::vector<int> &patch_nodes = ms_patches->at(patch_id);
	for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {
		int nb_patch_id = vert_patch_ids[*nb_it];
		if (patch_id != nb_patch_id && 
			std::find(patch_nodes.begin(), patch_nodes.end(), *nb_it) == patch_nodes.end()) {
			return true;
		}
	}
	return false;
}

void ComplexUtil::retraceSteepLines(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, 
	std::vector<int> &vert_patch_ids, patch_t ms_patches, steep_lines_t *steeplines, 
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {

	std::map<std::pair<int, int>, std::vector<int>> new_sls;
	Eigen::SparseMatrix<int> conn_graph(vertices_ptr->rows(), vertices_ptr->rows());
	
	vector<int> patches_bndry;

	for (int patch_id = 0; patch_id < ms_patches->size(); patch_id++) {
		vector<int> patch_nodes = ms_patches->at(patch_id);
		vector<int> visited;
		vector<int> patch_bndry;

		//int boundary_node = findVertOnBoundaryRec(vert_patch_ids, start_v, visited, HE_verts, HE_edges);
		// find patch_bndry

		for (int v_idx = 0; v_idx < vert_patch_ids.size(); v_idx++) {
			if (vert_patch_ids[v_idx] == patch_id ){
				if (onBoundary(ms_patches, v_idx, patch_id, vert_patch_ids, HE_verts, HE_edges)) {
					patch_bndry.push_back(v_idx);
				}
			}
		}
		for (auto n_it = patch_nodes.begin(); n_it != patch_nodes.end(); ++n_it) {
			if (onBoundary(ms_patches, *n_it, patch_id, vert_patch_ids, HE_verts, HE_edges)) {
				patch_bndry.push_back(*n_it);
			}
		}
		patches_bndry.insert(patches_bndry.end(), patch_bndry.begin(), patch_bndry.end());

		// check if patch node is surrounded, if yes, project onto boundary
		for (int i = 0; i < 4; i++) {
			int node_v = ms_patches->at(patch_id)[i];
			HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(node_v, HE_verts, HE_edges);
			bool node_inside = !onBoundary(ms_patches, node_v, patch_id, vert_patch_ids, HE_verts, HE_edges);

			// project node on boundary
			if (node_inside) {
				int close_bndry_v = findClosestBoundary(vertices_ptr, node_v,
					patch_bndry, HE_verts, HE_edges);
				updatePatchNode(ms_patches, node_v, close_bndry_v, vert_patch_ids, patch_id);
			}
		}
	}

	// build connectivity graph
	for (auto v = patches_bndry.begin(); v != patches_bndry.end(); ++v) {
		HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(*v, HE_verts, HE_edges);
		int num_nb_conn = 0;
		for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {
			if (std::find(patches_bndry.begin(), patches_bndry.end(), *nb_it) != patches_bndry.end()) {
				num_nb_conn++;
				conn_graph.coeffRef(*v, *nb_it) = 1;
				conn_graph.coeffRef(*nb_it, *v) = 1;
			}
		}
	}

	// find new sls
	for (int patch_id = 0; patch_id < ms_patches->size(); patch_id++) {
		for (int i = 0; i < 4; i++) {
			int node_v = ms_patches->at(patch_id)[i];
			int next_node_v = ms_patches->at(patch_id)[(i+1)%4];
			if (new_sls.find({ node_v, next_node_v }) == new_sls.end() &&
				new_sls.find({ next_node_v, node_v }) == new_sls.end()) {
				std::shared_ptr<std::vector<int>> path = PathFinding::getPath(node_v, next_node_v, conn_graph);
				new_sls[{node_v, next_node_v}] = *path;
			}
			
		}
		
	}

	steep_lines_t ptr = std::make_shared<std::map<std::pair<int, int>, std::vector<int>>>(new_sls);
	std::swap(ptr, *steeplines);
}

void insertSl2v(std::map<int, std::vector<int>> &sl_out_of_v, int start, int end) {
	if (sl_out_of_v.find(start) == sl_out_of_v.end()) {
		sl_out_of_v[start] = std::vector<int>(1, end);
	}
	else {
		sl_out_of_v[start].push_back(end);
	}
}

bool is_sl_rev(steep_lines_t steeplines, int start_vert, int end_vert) {
	auto sl_it = steeplines->find({ start_vert, end_vert });
	if (sl_it != steeplines->end()) {
		return false;
	}
	else {
		return true;
	}
}

void goToNextElem(std::vector<int>::iterator &it, bool rev) {
	if (rev) {
		it--;
	}
	else {
		it++;
	}
}

std::vector<int>::iterator getStartIt(std::vector<int> *sl, bool rev) {
	if (rev) {
		return (*sl).end() - 2;
	}
	else {
		// start at second elem, first always same
		return (*sl).begin() + 1;
	}
}

void relocateNode(std::vector<std::vector<int> *> sl_ptrs, std::vector<bool> revs, 
	int node_old_loc, int node_new_loc) {

}

void ComplexUtil::relocateNodes(steep_lines_t steeplines, patch_t ms_patches, 
	std::vector<int> &nodes, std::vector<int> &vert_patch_ids) {
	// get steeplines going out of a node
	std::map<int, std::vector<int>> sl_out_of_v;
	// map between current node and new node
	std::map<int, int> old_new_node;

	for (auto sl_it = steeplines->begin(); sl_it != steeplines->end(); ++sl_it) {
		int start = sl_it->first.first;
		int end = sl_it->first.second;
		insertSl2v(sl_out_of_v, start, end);
		insertSl2v(sl_out_of_v, end, start);
	}

	// check for overlap steeplines going out of a node
	for (auto node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
		vector<int> &ends = sl_out_of_v[*node_it];
		std::vector<std::vector<int> *> sl_ptrs(ends.size());
		std::vector<bool> revs(ends.size());
		// indices correspond to ends, min steps
		std::vector<int> min_steps_on_sl(ends.size(), 0);
		// indices correspond to ends, contens are the end ids that it overlaps with
		vector<vector<int>> overlap_ids(ends.size(), std::vector<int>({}));

		int e_id = 0;
		for (auto e_it = ends.begin(); e_it != ends.end(); ++e_it) {
			// end has been relocated
			if (old_new_node.find(*e_it) != old_new_node.end()) {
				*e_it = old_new_node[*e_it];
			}
			bool rev = is_sl_rev(steeplines, *node_it, *e_it);
			if (rev) {
				sl_ptrs[e_id] = &steeplines->at({*e_it, *node_it});
			}
			else {
				sl_ptrs[e_id] = &steeplines->at({ *node_it, *e_it });
			}
			revs[e_id] = rev;
			e_id++;
		}
		
		// compare every pair of sls going out of node
		for (int i = 0; i < ends.size() - 1; i++) {
			bool rev_i = revs[i];
			std::vector<int>::iterator start_i_it = getStartIt(sl_ptrs[i], rev_i);
			for (int j = i + 1; j < ends.size(); j++) {
				bool rev_j = revs[j];
				std::vector<int>::iterator start_j_it = getStartIt(sl_ptrs[j], rev_j);

				int steps_along_sl = 0;
				// start comparing, ASSUME two sls never overlap completely
				while (*start_i_it == *start_j_it) {
					steps_along_sl++;
					// get to next elem
					goToNextElem(start_i_it, rev_i);
					goToNextElem(start_j_it, rev_j);
				}
				// if there is overlap
				if (steps_along_sl > 0) {
					if (min_steps_on_sl[i] == 0 || steps_along_sl < min_steps_on_sl[i]) {
						min_steps_on_sl[i] = steps_along_sl;
					}
					if (min_steps_on_sl[j] == 0 || steps_along_sl < min_steps_on_sl[j]) {
						min_steps_on_sl[j] = steps_along_sl;
					}
					overlap_ids[i].push_back(j);
					overlap_ids[j].push_back(i);
				}
			}
		}
		
		int min_steps = 0;
		int move_sl_end_id = 0;
		for (int i = 0; i < ends.size(); i++) {
			if (min_steps_on_sl[i] > 0 && (min_steps == 0 || min_steps_on_sl[i] < min_steps)) {
				min_steps = min_steps_on_sl[i];
				move_sl_end_id = i;
			}
		}

		if (min_steps > 0) {
			bool move_sl_rev = revs[move_sl_end_id];
			auto it = getStartIt(sl_ptrs[move_sl_end_id], move_sl_rev);
			vector<int> overlapped(min_steps);
			for (int i = 0; i < min_steps; i++) {
				overlapped[i] = *it;
				goToNextElem(it, move_sl_rev);
			}
			// relocate node
			for (int i = 0; i < ends.size(); i++) {
				bool rev = revs[i];
				// if sl contains overlapped portion
				if (i == move_sl_end_id || 
					std::find(overlap_ids[i].begin(), overlap_ids[i].end(), move_sl_end_id) != overlap_ids[i].end()) {
					if (!rev) {
						sl_ptrs[i]->erase(sl_ptrs[i]->begin(), sl_ptrs[i]->begin() + min_steps);
					}
					else {
						sl_ptrs[i]->erase(sl_ptrs[i]->end() - min_steps, sl_ptrs[i]->end());
					}
					
				} 
				// if sl does not contain overlapped portion
				else {
					vector<int> rev_overlapped = overlapped;
					std::reverse(rev_overlapped.begin(), rev_overlapped.end());
					if (!rev) {
						sl_ptrs[i]->insert(sl_ptrs[i]->begin(), rev_overlapped.begin(), rev_overlapped.end());
					}
					else {
						sl_ptrs[i]->insert(sl_ptrs[i]->end(), overlapped.begin(), overlapped.end());
					}
				}
			}
			// update ms_patches
			int new_node = *(overlapped.end() - 1);
			updatePatchNode(ms_patches, *node_it, new_node, vert_patch_ids, -1);
			// update steeplines
			int e_id = 0;
			for (int i = 0; i < ends.size(); i++) {
				int end = ends[i];
				bool rev = is_sl_rev(steeplines, *node_it, end);
				if (rev) {
					(*steeplines)[{end, new_node}] = steeplines->at({ end, *node_it });
					steeplines->erase({ end, *node_it });
				}
				else {
					(*steeplines)[{new_node, end}] = steeplines->at({ *node_it, end });
					steeplines->erase({ *node_it, end });
				}
			}
			// update nodes
			old_new_node[*node_it] = new_node;
			//*node_it = new_node;
		}
	}
	for (auto n_it = old_new_node.begin(); n_it != old_new_node.end(); ++n_it) {
		auto node_it = std::find(nodes.begin(), nodes.end(), n_it->first);
		*node_it = n_it->second;
	}
}