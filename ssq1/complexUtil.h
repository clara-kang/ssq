#include <memory>
#include <vector>

#include "HalfEdge.h"
#pragma once
namespace ComplexUtil {
	enum VERT_TYPE { SADDLE, MAX, MIN, REG };

	typedef std::shared_ptr<std::vector< VERT_TYPE>> vert_type_t;
	typedef std::shared_ptr<std::map<std::pair<int, int>, std::vector<int>>> steep_lines_t;

	void findVertexTypes(HalfEdge::vert_t *saddles, HalfEdge::vert_t *maxs,
		HalfEdge::vert_t *mins, vert_type_t *vert_types,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U);

	steep_lines_t findSteepLines(HalfEdge::vert_t saddles, vert_type_t vert_types,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U);
}