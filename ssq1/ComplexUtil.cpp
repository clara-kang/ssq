#include <memory>
#include <vector>
#include "complexUtil.h"
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
		//	values[nb_idx] = U(nb_idx);
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
		//	double nb_val = U(neighbors->at(nb_idx));
		//	sign = (nb_val > v_val);
		//	if (sign != last_sign) {
		//		sign_change_times++;
		//	}
		//	last_sign = sign;
		//}
		//// last one
		//if (last_sign != first_sign) {
		//	sign_change_times++;
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
	//	cout << "nb " << *it << ": " << U(*it) << endl;
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
		//	cout << "nb: " << neighbors->at(nb_id) << ", U: " << U(neighbors->at(nb_id)) << endl;
		//}
		for (int nb_id = 0; nb_id < neighbors->size(); nb_id ++) {
			int nb = neighbors->at(nb_id);
			int next_nb = neighbors->at((nb_id+1)% neighbors->size());
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

		s2e_angl.insert(s2e_angl.begin(), { it->first, s2e_prpndclr});
		s2e_angl.insert(s2e_angl.begin(), 
			{ std::make_pair(it->first.second, it->first.first), e2s_prpndclr });
	}

	// build v2_sl, the adjacency list
	for (auto it = s2e_angl.begin(); it != s2e_angl.end(); ++it) {
		if (v2_sl.find(it->first.first) == v2_sl.end()) {
			v2_sl.insert(v2_sl.begin(), { it->first.first, std::vector<int>({it->first.second})});
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
			if ( !sl_cover[{v2_sl_it->first, *ends_it}]) {

				int end_vert = *ends_it;
				int start_vert = v2_sl_it->first;
				std::vector<int> patch({start_vert, end_vert});

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
							double angle_diff = abs(PI/2.0 - angle);
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
	
	for (int patch_cnt = 0; patch_cnt < ms_patches->size(); patch_cnt ++) {
		for (int i = 0; i < 4; i++) {
			int s = ms_patches->at(patch_cnt)[i];
			int e = ms_patches->at(patch_cnt)[(i+1)%4];
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
	vector<int> &ids = sl_patch_map->at(getSlKey(steep_lines, {start, end}));
	if (ids[0] == self_id) {
		return ids[1];
	}
	return ids[0];
}

// patch corner 0,1,2,3 corresponds to (0,0), (1,0), (1,1), (0,1)
trnsfr_funcs_map_t ComplexUtil::buildTrnsfrFuncsAndPatchGraph(steep_lines_t steeplines, 
	steep_lines_t sl_patch_map, patch_t ms_patches, patch_t *patch_graph_in) {
	std::map<std::pair<int, int>, trnsfr_func_t> trnsfr_funcs_map;
	std::vector<std::vector<int>> patch_graph(ms_patches->size(), std::vector<int>({0, 0, 0, 0}));
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

bool ComplexUtil::isPatchNode(patch_t ms_patches, int vert_idx, int patch_id) {
	vector<int> &patch_nodes = ms_patches->at(patch_id);
	auto it = std::find(patch_nodes.begin(), patch_nodes.end(), vert_idx);
	if (it != patch_nodes.end()) {
		return true;
	}
	return false;
}

void ComplexUtil::applyTrnsfrFunc(Eigen::Vector2d &initial_uv, 
	Eigen::Vector2d &res_uv, trnsfr_func_t trnsfr_func) {

	if (std::get<0>(trnsfr_func)) {
		res_uv << initial_uv(1), initial_uv(0);
	}
	else {
		res_uv << initial_uv(0), initial_uv(1);
	}
	res_uv(0) *= std::get<2>(trnsfr_func).first;
	res_uv(1) *= std::get<2>(trnsfr_func).second;
	res_uv(0) *= std::get<2>(trnsfr_func).first;
	res_uv(1) *= std::get<2>(trnsfr_func).second;
}

std::shared_ptr<Eigen::MatrixXd> ComplexUtil::solveForCoords(Eigen::SparseMatrix<double> &L, std::vector<int> &vert_patch_ids,
	ComplexUtil::trnsfr_funcs_map_t trnsfr_funcs_map, int node_num, patch_t ms_patches,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges) {

	int unknown_size = vert_patch_ids.size() - node_num;
	Eigen::SparseMatrix<double> M(unknown_size*2, unknown_size * 2);
	Eigen::VectorXd knowns(unknown_size * 2);

	for (int i = 0; i < M.rows(); i++) {
		knowns(i) = 0;
	}

	std::vector<int> verts_matrx_map(vert_patch_ids.size(), -1);

	// build relationship between vertex index to matrix row
	int row_indx = 0;
	for (int v_indx = 0; v_indx < vert_patch_ids.size(); v_indx++) {
		if (!isPatchNode(ms_patches, v_indx, vert_patch_ids[v_indx])) {
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

			double nb_weight = - L.coeff(vert_indx, *nb_it) / node_weight;
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
				}
				else {
					M.insert(v_mmap_indx_u, nb_mmap_indx_x) = -nb_weight;
					M.insert(v_mmap_indx_v, nb_mmap_indx_y) = -nb_weight;
				}
			}
			else {
				int nb_node_indx = which_node(ms_patches->at(nb_patch_id), *nb_it);
				if (nb_node_indx == 1 || nb_node_indx == 2) {
					knowns(v_mmap_indx_u) += nb_weight;
				}
				if (nb_node_indx == 2 || nb_node_indx == 3) {
					knowns(v_mmap_indx_u) += nb_weight;
				}
			}
		}
		//auto M_dense = Eigen::MatrixXd(M);
		//std::cout << M_dense.row(0) << std::endl;

	}

	//std::cout << "knowns: " << knowns << endl;
	M.makeCompressed();

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solverA;
	solverA.compute(M);

	if (solverA.info() != Eigen::Success) {
		cout << "Oh: Very bad" << endl;
	}

	Eigen::VectorXd X = Eigen::VectorXd(unknown_size * 2);
	X = solverA.solve(knowns);
	for (int indx = 0; indx < X.rows(); indx++) {
		std::cout << X(indx) << std::endl;
	}

	Eigen::MatrixXd uv_coords(vert_patch_ids.size(), 2);

	for (int v_indx = 0; v_indx < vert_patch_ids.size(); v_indx++) {
		int v_mmap_indx_u = verts_matrx_map[v_indx];
		if (v_mmap_indx_u >= 0) {
			int v_mmap_indx_v = v_mmap_indx_u + unknown_size;
			uv_coords.row(v_indx) << X(v_mmap_indx_u), X(v_mmap_indx_v);
		}
		else {
			uv_coords.row(v_indx) << -1, -1;
		}
	}

	std::shared_ptr<Eigen::MatrixXd> uv_coords_ptr = std::make_shared<Eigen::MatrixXd>(uv_coords);
	return uv_coords_ptr;
}

