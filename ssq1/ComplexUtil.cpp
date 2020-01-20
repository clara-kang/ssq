#include <memory>
#include <vector>
#include "complexUtil.h"
#include <algorithm>
#include <map>
//#include <tuple>

using namespace ComplexUtil;
using namespace std;

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
