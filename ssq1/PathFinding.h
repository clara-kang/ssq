#pragma once
#pragma once

#include <vector>
#include <memory>
#include <map>
#include <Eigen/SparseExtra>
#include "HalfEdge.h"

std::shared_ptr<std::vector<int>> trace_path(std::map<int, int> &prevs, int src_v, int dest_v) {
	std::vector<int> path;
	int v = dest_v;
	while (v != src_v) {
		path.push_back(v);
		v = prevs[v];
	}
	path.push_back(v);
	std::reverse(path.begin(), path.end());
	std::shared_ptr<std::vector<int>> path_ptr = std::make_shared<std::vector<int>>(path);
	return path_ptr;
}

std::shared_ptr<std::vector<int>> getSubPath(int src_v, int dest_v, std::vector<std::vector<int>> &conn_map,
	std::map<int, int> &prevs, std::vector<int> &work_list) {
	int from_v = work_list[0];
	work_list.erase(work_list.begin());

	std::vector<int> &nbs = conn_map.at(from_v);

	for (auto nb_it = nbs.begin(); nb_it != nbs.end(); ++nb_it) {
		if (*nb_it == dest_v) {
			// found
			prevs[*nb_it] = from_v;
			return trace_path(prevs, src_v, dest_v);
		}
		if (*nb_it != dest_v && prevs.find(*nb_it) == prevs.end()) {
			prevs[*nb_it] = from_v;
			work_list.push_back(*nb_it);
		}
	}

	return nullptr;
}

std::shared_ptr<std::vector<int>> getSubPath(int src_v, int dest_v, Eigen::SparseMatrix<int> &conn_map,
	std::map<int, int> &prevs, std::vector<int> &work_list) {
	int from_v = work_list[0];
	work_list.erase(work_list.begin());

	std::vector<int> nbs;
	for (Eigen::SparseMatrix<int>::InnerIterator it(conn_map, from_v); it; ++it) {
		nbs.push_back(it.row());
	}

	for (auto nb_it = nbs.begin(); nb_it != nbs.end(); ++nb_it) {
		if (*nb_it == dest_v) {
			// found
			prevs[*nb_it] = from_v;
			return trace_path(prevs, src_v, dest_v);
		}
		if (*nb_it != dest_v && prevs.find(*nb_it) == prevs.end()) {
			prevs[*nb_it] = from_v;
			work_list.push_back(*nb_it);
		}
	}

	return nullptr;
}

bool getPathToPatchRec(int start, int &dst, std::vector<int> &vert_patch_ids, int dest_patch_id,
	HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges,
	vector<int> &work_list, map<int, int> &prevs) {

	HalfEdge::vert_t neighbors = HalfEdge::getNeighbors(start, HE_verts, HE_edges);

	for (auto nb_it = neighbors->begin(); nb_it != neighbors->end(); ++nb_it) {
		if (vert_patch_ids[*nb_it] == dest_patch_id) {
			prevs[*nb_it] = start;
			dst = *nb_it;
			return true;
		}
		if (prevs.find(*nb_it) == prevs.end()) {
			prevs[*nb_it] = start;
			work_list.push_back(*nb_it);
		}
	}
	return false;
}

namespace PathFinding {
	std::shared_ptr<std::vector<int>> getPath(int src_v, int dest_v, std::vector<std::vector<int>> &conn_map) {

		std::vector<int> work_list(1, src_v);
		std::map<int, int> prevs;
		prevs[src_v] = src_v;

		while (work_list.size() != 0) {
			std::shared_ptr<std::vector<int>> path_ptr = getSubPath(src_v, dest_v, conn_map, prevs, work_list);
			if (path_ptr != nullptr) {
				return path_ptr;
			}
		}

		return nullptr;

	}

	std::shared_ptr<std::vector<int>> getPath(int src_v, int dest_v, Eigen::SparseMatrix<int> &conn_map) {

		std::vector<int> work_list(1, src_v);
		std::map<int, int> prevs;
		prevs[src_v] = src_v;

		while (work_list.size() != 0) {
			std::shared_ptr<std::vector<int>> path_ptr = getSubPath(src_v, dest_v, conn_map, prevs, work_list);
			if (path_ptr != nullptr) {
				return path_ptr;
			}
		}

		return nullptr;

	}

	void getPathToPatch(int start, std::vector<int> &vert_patch_ids, int dest_patch_id,
		HalfEdge::vert_t HE_verts, HalfEdge::edge_t HE_edges, vector<int> &path) {

		int dest;
		map<int, int> prevs;
		vector<int> work_list;

		work_list.push_back(start);
		while (work_list.size() > 0) {
			int next_start = work_list[0];
			work_list.erase(work_list.begin());
			if (getPathToPatchRec(next_start, dest, vert_patch_ids, dest_patch_id,
				HE_verts, HE_edges, work_list, prevs)) {

				path = *(trace_path(prevs, start, dest));
				return;
			}
		}
	}
}

