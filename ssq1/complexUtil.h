#include <memory>
#include <vector>
#include <tuple>
#include "HalfEdge.h"
#pragma once

using namespace std;

namespace ComplexUtil {
	enum VERT_TYPE { SADDLE, MAX, MIN, REG };

	typedef std::shared_ptr<std::vector< VERT_TYPE>> vert_type_t;
	typedef std::shared_ptr<std::map<std::pair<int, int>, std::vector<int>>> steep_lines_t;
	typedef std::shared_ptr<std::vector<std::vector<int>>> patch_t; // the 4 verts
	typedef std::tuple<bool, std::pair<int, int>, std::pair<int, int>> trnsfr_func_t;
	typedef std::shared_ptr<std::map<std::pair<int, int>, trnsfr_func_t>> trnsfr_funcs_map_t;


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

	bool complexValidityCheck(patch_t ms_patches);

	bool patchValidityCheck(std::shared_ptr<vector<vector<int>>> patch_verts,
		int N, const vector<int> &sl_indices);

	// return map, key is {sl_start_vert, sl_end_vert}, value is {patch1, patch2}
	steep_lines_t buildSLtoPatchMap(steep_lines_t steeplines, patch_t ms_patches);

	trnsfr_funcs_map_t buildTrnsfrFuncsAndPatchGraph(steep_lines_t steeplines,
		steep_lines_t sl_patch_map, patch_t ms_patches, patch_t *patch_graph);

	void buildExtendedTrnsfrfuncs(patch_t ms_patches, patch_t patch_graph, trnsfr_funcs_map_t trnsfr_funcs_map);

	//// return vectors, each is [patch1, patch2, patch3, patch3]
	//patch_t buildPatchGraph(steep_lines_t sl_patch_map, patch_t ms_patches);

	std::shared_ptr<Eigen::MatrixXd> solveForCoords(Eigen::SparseMatrix<double> &L, std::vector<int> &vert_patch_ids,
		ComplexUtil::trnsfr_funcs_map_t trnsfr_funcs_map, int node_num, patch_t ms_patches,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges);

	bool isPatchNode(patch_t ms_patches, int vert_idx, int patch_id);

	int which_node(std::vector<int> &patch, int vert);

	void applyTrnsfrFunc(Eigen::Vector2d &initial_uv, Eigen::Vector2d &res_uv, trnsfr_func_t trnsfr_func);

	void adjustBndrys(std::vector<int> &vert_patch_ids, steep_lines_t steeplines,
		patch_t patch_graph, trnsfr_funcs_map_t trnsfr_functions,
		std::shared_ptr<Eigen::MatrixXd> uvs, std::shared_ptr<vector<vector<int>>> patch_verts,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges);

	void retraceSteepLines(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, 
		std::vector<int> &vert_patch_ids, patch_t ms_patches,
		steep_lines_t steeplines, std::shared_ptr<vector<vector<int>>> patch_verts,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges);
}
