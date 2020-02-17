#pragma once
#include <memory>
#include <Eigen/SparseExtra>
#include <vector>

namespace HalfEdge {
	struct HE {
		int vert; // end vert
		int pair;
		int face;
		int next;
	};

	typedef std::shared_ptr<std::vector<int>> vert_t;
	typedef std::shared_ptr<std::vector<HE>> edge_t;

	//void hi();
	void createHEs(Eigen::MatrixXi &face, int N, std::shared_ptr<std::vector<HE>> *HE_edges,  
		std::shared_ptr<std::vector<int>> *HE_verts, std::shared_ptr<std::vector<int>> *HE_faces);

	// return neighbors in clockwise order??
	vert_t getNeighbors(int v_index, vert_t HE_vert, edge_t HE_edges);

	vert_t getFacesOfV(int v_index, vert_t HE_vert, edge_t HE_edges);
}
