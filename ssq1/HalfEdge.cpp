#include "HalfEdge.h"

void HalfEdge::createHEs(Eigen::MatrixXi &faces, int N, std::shared_ptr<std::vector<HalfEdge::HE>> *HE_edges,
	std::shared_ptr<std::vector<int>> *HE_verts, std::shared_ptr<std::vector<int>> *HE_faces) {
	std::map<std::pair<int, int>, int> edges_map;
	std::vector< HalfEdge::HE> he_edges(faces.rows() * 3);
	std::vector<int> he_faces(faces.rows()); // points to an edge
	std::vector<int> he_verts(N); // points to an edge that starts with v

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int v = faces(i, j);
			int v_next = faces(i, (j + 1) % 3);
			int edges_index = i * 3 + j;
			HE he;
			he.vert = v_next;
			he.face = i;
			he.next = i * 3 + (j + 1) % 3;
			edges_map[std::pair<int, int>(v, v_next)] = edges_index;
			he_edges[edges_index] = he;
			he_verts[v] = edges_index;
		}
		he_faces[i] = i * 3;
	}
	// get pairs
	std::map<std::pair<int, int>, int>::iterator it;
	for (it = edges_map.begin(); it != edges_map.end(); ++it) {
		int edge_index = it->second;
		// not yet have a pair
		if (he_edges[edge_index].pair < 0) {
			int start_v = it->first.first;
			int end_v = it->first.second;
			auto pair_index_it = edges_map.find(std::pair<int, int>(end_v, start_v));
			if (pair_index_it != edges_map.end()) {
				he_edges[edge_index].pair = pair_index_it->second;
			}
		}
	}
	auto done_HE_edges = std::make_shared<std::vector<HalfEdge::HE>>(he_edges);
	auto done_HE_verts = std::make_shared < std::vector <int>>(he_verts);
	auto done_HE_faces = std::make_shared < std::vector <int>>(he_faces);

	swap(done_HE_edges, *HE_edges);
	swap(done_HE_verts, *HE_verts);
	swap(done_HE_faces, *HE_faces);

}

HalfEdge::vert_t HalfEdge::getNeighbors(int v_index, HalfEdge::vert_t HE_vert, HalfEdge::edge_t HE_edges) {
	vert_t neighbors = std::make_shared<std::vector<int>>();
	int he_index = HE_vert->at(v_index);
	// append first neighbor
	neighbors->push_back(HE_edges->at(he_index).vert);

	int next_he_indx = HE_edges->at(HE_edges->at(he_index).pair).next;
	int nb_indx;
	while (next_he_indx != he_index) {
		nb_indx = HE_edges->at(next_he_indx).vert;
		neighbors->push_back(nb_indx);
		next_he_indx = HE_edges->at(HE_edges->at(next_he_indx).pair).next;
	}
	return neighbors;
}

