#include <memory>
#include <vector>

#include "HalfEdge.h"
#pragma once

using namespace std;

namespace ComplexUtil {
	enum VERT_TYPE { SADDLE, MAX, MIN, REG };

	typedef std::shared_ptr<std::vector< VERT_TYPE>> vert_type_t;
	typedef std::shared_ptr<std::map<std::pair<int, int>, std::vector<int>>> steep_lines_t;
	typedef std::shared_ptr<std::vector<std::vector<int>>> patch_t;

	void findVertexTypes(HalfEdge::vert_t *saddles, HalfEdge::vert_t *maxs,
		HalfEdge::vert_t *mins, vert_type_t *vert_types,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U);

	steep_lines_t findSteepLines(HalfEdge::vert_t saddles, vert_type_t vert_types,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, Eigen::MatrixXd &U);

	patch_t buildMsPatches(steep_lines_t steeplines, std::shared_ptr<Eigen::MatrixXd> vertices_ptr,
		std::shared_ptr<Eigen::MatrixXd> vns_ptr);


	void fillMsPatches(steep_lines_t steeplines, patch_t ms_patches,
		HalfEdge::vert_t HE_vert, HalfEdge::edge_t HE_edges,
		std::shared_ptr<Eigen::MatrixXd> vertices_ptr,
		std::shared_ptr<Eigen::MatrixXd> vns_ptr,
		std::shared_ptr<vector<vector<int>>> *patch_verts);

	bool patchValidityCheck(std::shared_ptr<vector<vector<int>>> patch_verts,
		int N, const vector<int> &sl_indices);
}