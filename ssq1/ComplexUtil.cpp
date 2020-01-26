#include <memory>
#include <vector>
#include "complexUtil.h"
#include <algorithm>
#include <map>
#include <cmath>

//#include <tuple>

using namespace ComplexUtil;
using namespace std;

double PI = 3.1415926;

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

std::vector<std::vector<int>> ComplexUtil::buildMsPatches(steep_lines_t steeplines, std::shared_ptr<Eigen::MatrixXd> vertices_ptr,
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
		Eigen::Vector3d e2s = -s2e;
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
							Eigen::Vector3d v1 = -s2e_angl[{end_vert, start_vert}];
							Eigen::Vector3d v2 = s2e_angl[{end_vert, *end_nb_it}];
							Eigen::Vector3d end_norm = vns_ptr->row(end_vert);
						
							// get angle
							double angle = acos(v2.dot(v1) / (v1.norm() * v2.norm()));
							// if angle should be > 180
							if (end_norm.dot(v1.cross(v2)) < 0) {
								angle = 2.0 * PI - angle;
							}
							// diff between angle and 90 degrees
							double angle_diff = abs(PI/2.0 - angle);
							end_angles.insert(end_angles.begin(), { *end_nb_it, angle_diff });
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
				ms_patches.push_back(patch);
			}
		}
	}
	return ms_patches;

}
