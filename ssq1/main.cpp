#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseExtra>
#include <set>
#include <string.h>

#include "HalfEdge.h"
#include "complexUtil.h"

using namespace std;

//int EIG_NUM = 10;
int EIG_NUM = 10;

bool LOAD_U = true;
bool LOAD_UV = false;

string u_file_name = "../U_sphere.txt";

void clamp(double &x, double low, double high) {
	if (x <= low) {
		x = low;
	}
	else if (x >= high) {
		x = high;
	}
}

vector<int> getPatchNodes(ComplexUtil::patch_t ms_patches) {
	// gather nodes
	set<int> patch_nodes;
	for (auto p_it = ms_patches->begin(); p_it != ms_patches->end(); ++p_it) {
		for (int i = 0; i < 4; i++) {
			int node_id = p_it->at(i);
			patch_nodes.insert(node_id);
		}
	}

	vector<int> patch_nodes_v(patch_nodes.begin(), patch_nodes.end());
	return patch_nodes_v;
}

void setNodesToNeg(vector<int> &patch_nodes_v, std::vector<int> &vert_patch_ids) {
	for (auto node_it = patch_nodes_v.begin(); node_it != patch_nodes_v.end(); ++node_it) {
		vert_patch_ids[*node_it] = -1;
	}
}

void writeUout(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, Eigen::MatrixXd &U) {

	// write U out as colors
	std::ofstream col_out("../cols.txt");
	col_out << "cols = [None] * " + std::to_string(vertices_ptr->rows()) << endl;
	// assign a color to each vertex based on U, blue -> green -> red
	double red[3] = { 1.0, 0.0, 0.0 };
	double green[3] = { 0.0, 1.0, 0.0 };
	double blue[3] = { 0.0, 0.0, 1.0 };
	for (int v_indx = 0; v_indx < vertices_ptr->rows(); v_indx++) {
		double v_col[3];
		// interpolate from blue to green
		if (U(v_indx) < 0.5) {
			for (int i = 0; i < 3; i++) {
				v_col[i] = 2.0 * (U(v_indx) * green[i] + (0.5 - U(v_indx)) * blue[i]);
			}
		}
		// interpolate from green to red
		else {
			for (int i = 0; i < 3; i++) {
				double val = U(v_indx) - 0.5;
				v_col[i] = 2.0 * ((0.5 - val) * green[i] + val * red[i]);
			}
		}
		col_out << "cols[" << v_indx << "]=(" << v_col[0] << ", " << v_col[1] << ", " << v_col[2] << ")" << endl;
	}
}

void writeMaxMinSdlOut(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, HalfEdge::vert_t maxs, 
	HalfEdge::vert_t mins, HalfEdge::vert_t saddles) {
	// iterate over points
	std::ofstream pts_out("../points_verts.txt");
	int pts_cnt = 0;
	pts_out << "max_locs = [None] * " + std::to_string(maxs->size()) + "\n";
	pts_out << "sdl_locs = [None] * " + std::to_string(saddles->size()) + "\n";
	pts_out << "min_locs = [None] * " + std::to_string(mins->size()) + "\n";
	for (int pt_type = 0; pt_type < 3; pt_type++) {
		string array_name;
		pts_cnt = 0;
		std::shared_ptr<vector<int>> pt_list;
		if (pt_type == 0) {
			pt_list = maxs;
			array_name = "max_locs";
		}
		else if (pt_type == 1) {
			pt_list = saddles;
			array_name = "sdl_locs";
		}
		else {
			pt_list = mins;
			array_name = "min_locs";
		}
		for (auto v_it = pt_list->begin(); v_it != pt_list->end(); ++v_it) {
			pts_out << array_name + "[" + std::to_string(pts_cnt) + "] = ("
				+ std::to_string((*vertices_ptr)(*v_it, 0)) + ", "
				+ std::to_string((*vertices_ptr)(*v_it, 1)) + ", "
				+ std::to_string((*vertices_ptr)(*v_it, 2)) + ")\n";
			pts_cnt++;
		}
	}

	pts_out.close();
}

void writeLineOUt(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, ComplexUtil::steep_lines_t steeplines) {
	// iterate over steep lines
	std::ofstream lines_out("../line_verts.txt");
	int lines_cnt = 0;
	lines_out << "verts = [None] * " + std::to_string(steeplines->size()) << endl;
	for (auto sl_it = steeplines->begin(); sl_it != steeplines->end(); ++sl_it) {
		vector<int> &sl = sl_it->second;
		lines_out << "verts[" + std::to_string(lines_cnt) + "] = [";
		for (auto v_it = sl.begin(); v_it != sl.end(); ++v_it) {
			auto v_loc = (*vertices_ptr).row(*v_it);
			lines_out << "(" + std::to_string((*vertices_ptr)(*v_it, 0)) + ", "
				+ std::to_string((*vertices_ptr)(*v_it, 1)) + ", "
				+ std::to_string((*vertices_ptr)(*v_it, 2)) + "),";
		}
		lines_out << +"]\n";
		lines_cnt++;
	}
	lines_out.close();
}

void writePatchOut(std::vector<int> &vert_patch_ids) {
	// write patch verts out
	std::ofstream patches_out("../patches.txt");
	patches_out << "patch_ids = [";
	for (int i = 0; i < vert_patch_ids.size(); i++) {
		patches_out << vert_patch_ids[i] << ", ";
	}
	patches_out << "]" << endl;
}

void writeNodesOut(std::shared_ptr<Eigen::MatrixXd> vertices_ptr, vector<int> &patch_nodes_v) {
	// write out patch nodes
	std::ofstream nodes_out("../patch_nodes.txt");
	nodes_out << "nodes_locs = [";
	for (auto node_it = patch_nodes_v.begin(); node_it != patch_nodes_v.end(); ++node_it) {
		Eigen::Vector3d loc = vertices_ptr->row(*node_it);
		nodes_out << "(" << loc(0) << ", " << loc(1) << ", " << loc(2) << "), ";
	}
	nodes_out << "]" << endl;
	nodes_out.close();
}

int main(int argc, char *argv[])
{

	//auto vertices_ptr = std::make_shared<Eigen::MatrixXd>();
	//auto faces_ptr = std::make_shared<Eigen::MatrixXi>();

	//igl::readOBJ("../models/sphere3.obj", *vertices_ptr, *faces_ptr);

	Eigen::MatrixXd tcs;
	Eigen::MatrixXi ftcs;

	auto vertices_ptr = std::make_shared<Eigen::MatrixXd>();
	auto vns_ptr = std::make_shared<Eigen::MatrixXd>();
	auto N = std::make_shared<Eigen::MatrixXd>();
	auto faces_ptr = std::make_shared<Eigen::MatrixXi>();
	auto fns_ptr = std::make_shared<Eigen::MatrixXi>();

	// //Load a mesh in OFF format
	igl::readOBJ("../models/high_sphere.obj", *vertices_ptr, tcs, *N, *faces_ptr, ftcs, *fns_ptr);

	// todo: fix later
	vns_ptr = vertices_ptr;
	cout << "vertices_ptr size: " << vertices_ptr->rows() << endl;
	cout << "vns_ptr size: " << vns_ptr->rows() << endl;
	//for (int i = 0; i < vns_ptr->rows(); i++) {
	// vns_ptr->row(i) = vertices_ptr->row(i).normalized();
	// cout << "vn[" << i << "]: " << vns_ptr->row(i) << endl;;
	//}

	Eigen::SparseMatrix<double> L;
	Eigen::SparseMatrix<double> M(vertices_ptr->rows(), vertices_ptr->rows());
	igl::cotmatrix(*vertices_ptr, *faces_ptr, L);
	L = (-L).eval();

	//igl::massmatrix(*vertices_ptr, *faces_ptr, igl::MASSMATRIX_TYPE_DEFAULT, M);
	M.setIdentity();
	

	Eigen::MatrixXd U;
	Eigen::VectorXd D;
	//igl::eigs(L, M, 5, igl::EIGS_TYPE_SM, U, D);
	if (LOAD_U) {
		U = Eigen::MatrixXd(vertices_ptr->rows(), 1);
		std::string line;
		std::ifstream rfile;
		rfile.open(u_file_name);
		if (rfile.is_open()) {
			int line_num = 0;
			while (std::getline(rfile, line)) {
				U(line_num) = stod(line);
				//cout << "U(line_num): " << U(line_num) << endl;
				line_num++;
			}
			rfile.close();
		}
	}
	else {
		if (!igl::eigs(L, M, EIG_NUM, igl::EIGS_TYPE_SM, U, D))
		{
			cout << "failed." << endl;
			exit(1);
		}
		std::ofstream u_file(u_file_name);
		U = ((U.array() - U.minCoeff()) / (U.maxCoeff() - U.minCoeff())).eval();
		for (int i = 0; i < U.rows(); i++) {
			u_file << U(i) << endl;
		}
		//std::cout << "U: " << U << endl;
		//u_file << U << endl;
		u_file.close();
	}

	//writeUout(vertices_ptr, U);

	std::shared_ptr< std::vector<HalfEdge::HE>> HE_edges = nullptr;
	std::shared_ptr<std::vector<int>> HE_verts = nullptr, HE_faces = nullptr;
	HalfEdge::createHEs(*faces_ptr, vertices_ptr->rows(), &HE_edges, &HE_verts, &HE_faces);

	// smooth U
	//Eigen::MatrixXd U_smooth = Eigen::MatrixXd(vertices_ptr->rows(), 1);
	//for (int i = 0; i < 1; i++) {
	//	ComplexUtil::smoothU(HE_verts, HE_edges, U, U_smooth);
	//	U = U_smooth;
	//}
	writeUout(vertices_ptr, U);

	HalfEdge::vert_t saddles = nullptr;
	HalfEdge::vert_t maxs = nullptr;
	HalfEdge::vert_t mins = nullptr;
	ComplexUtil::vert_type_t vert_types = nullptr;
	ComplexUtil::findVertexTypes(&saddles, &maxs, &mins, &vert_types, HE_verts, HE_edges, U);

	ComplexUtil::steep_lines_t steeplines;
	steeplines = ComplexUtil::findSteepLines(saddles, vert_types, HE_verts, HE_edges, U);
	// filter by persistence
	ComplexUtil::filterSteepLines(steeplines, U, vert_types, maxs, mins, saddles);

	writeMaxMinSdlOut(vertices_ptr, maxs, mins, saddles);
	writeLineOUt(vertices_ptr, steeplines);

	ComplexUtil::patch_t ms_patches = ComplexUtil::buildMsPatches(steeplines, vertices_ptr, vns_ptr);
	bool complex_valid = ComplexUtil::complexValidityCheck(ms_patches);

	std::shared_ptr<vector<vector<int>>> patch_verts = nullptr;
	std::vector<int> vert_patch_ids(vertices_ptr->rows(), -1);

	ComplexUtil::fillMsPatches(steeplines, ms_patches, HE_verts, HE_edges, vertices_ptr, vns_ptr, &patch_verts);

	// all vertices on steeplines, including duplicates
	vector<int> verts_on_sls;
	for (auto sl_it = steeplines->begin(); sl_it != steeplines->end(); ++sl_it) {
		vector<int> &sl = sl_it->second;
		for (auto v_it = sl.begin(); v_it != sl.end(); ++v_it) {
			verts_on_sls.push_back(*v_it);
		}
	}

	bool patches_valid = ComplexUtil::patchValidityCheck(patch_verts, vertices_ptr->rows(), verts_on_sls);

	ComplexUtil::steep_lines_t sl_patch_map = ComplexUtil::buildSLtoPatchMap(steeplines, ms_patches);

	ComplexUtil::patch_t patch_graph = nullptr;
	ComplexUtil::trnsfr_funcs_map_t trnsfr_funcs_map =
		ComplexUtil::buildTrnsfrFuncsAndPatchGraph(steeplines, sl_patch_map, ms_patches, &patch_graph);

	ComplexUtil::buildExtendedTrnsfrfuncs(ms_patches, patch_graph, trnsfr_funcs_map);

	int node_num = maxs->size() + mins->size() + saddles->size();
	vector<int> patch_nodes_v = getPatchNodes(ms_patches);

	//write patch verts
	int patch_cnt = 0;
	for (auto patch_it = patch_verts->begin(); patch_it != patch_verts->end(); ++patch_it) {
		for (auto vert_it = patch_it->begin(); vert_it != patch_it->end(); ++vert_it) {
			vert_patch_ids[*vert_it] = patch_cnt;
		}
		patch_cnt++;
	}

	setNodesToNeg(patch_nodes_v, vert_patch_ids);
	writePatchOut(vert_patch_ids);

	std::shared_ptr<Eigen::MatrixXd> uv_coords;
	//if (!LOAD_UV) {
	//	uv_coords = ComplexUtil::solveForCoords(L, vert_patch_ids, trnsfr_funcs_map, node_num, ms_patches,
	//		HE_verts, HE_edges);
	//	// write UVs out, for debugging only! 
	//	std::ofstream uv_file("../UV.txt");
	//	for (int i = 0; i < U.rows(); i++) {
	//		uv_file << (*uv_coords)(i, 0) << " " << (*uv_coords)(i, 1) << endl;
	//	}
	//	uv_file.close();
	//}
	//else {
	//	Eigen::MatrixXd uv = Eigen::MatrixXd(vertices_ptr->rows(), 2);
	//	std::string line;
	//	std::ifstream rfile;
	//	rfile.open("../UV.txt");
	//	if (rfile.is_open()) {
	//		int line_num = 0;
	//		while (std::getline(rfile, line)) {
	//			size_t pos = line.find(" ");
	//			std::string u_str = line.substr(0, pos);
	//			std::string v_str = line.substr(pos + 1, line.length());
	//			uv.row(line_num) << stod(u_str), stod(v_str);
	//			line_num++;
	//		}
	//		rfile.close();
	//	}
	//	uv_coords = std::make_shared<Eigen::MatrixXd>(uv);
	//}

	for (int i = 0; i < 1; i++) {

		uv_coords = ComplexUtil::solveForCoords(L, vert_patch_ids, trnsfr_funcs_map, node_num, ms_patches,
			HE_verts, HE_edges);


		//vector<vector<double>> uv_vec(vertices_ptr->rows());
		for (int i = 0; i < vertices_ptr->rows(); i++) {
			//uv_vec[i] = vector<double>({ (*uv_coords)(i, 0), (*uv_coords)(i, 1) });
			if (vert_patch_ids[i] != -1) {
				if (abs((*uv_coords)(i, 0)) > 2 || abs((*uv_coords)(i, 1)) > 2) {
					cout << "error: " << (*uv_coords)(i, 0) << ", " << (*uv_coords)(i, 1) << endl;
				}
			}
		}
		ComplexUtil::adjustBndrys(vert_patch_ids, ms_patches, patch_graph, trnsfr_funcs_map,
			uv_coords, patch_verts, HE_verts, HE_edges);

		//map<int, set<int>> node2patch;
		//ComplexUtil::buildNode2PatchesMap(ms_patches, node2patch);

		//ComplexUtil::relocatePatchNodes(vertices_ptr, vert_patch_ids, ms_patches,
		//	node2patch, uv_coords, HE_verts, HE_edges);

		//// gather nodes
		//patch_nodes_v = getPatchNodes(ms_patches);
		//setNodesToNeg(patch_nodes_v, vert_patch_ids);
		//uv_coords = ComplexUtil::solveForCoords(L, vert_patch_ids, trnsfr_funcs_map, node_num, ms_patches,
		//	HE_verts, HE_edges);
	}


	//------------------------------------------------------------------//
	writePatchOut(vert_patch_ids);
	writeNodesOut(vertices_ptr, patch_nodes_v);

	// write indices of verts on steeplines
	std::ofstream sl_vert_id_out("../sl_vert_ids.txt");
	sl_vert_id_out << "sl_verts = [";

	for (auto v_it = verts_on_sls.begin(); v_it != verts_on_sls.end(); ++v_it) {
		sl_vert_id_out << *v_it << ", ";
	}

	sl_vert_id_out << "]" << endl;




	// write texture coords out
	Eigen::MatrixXd tcs_out(faces_ptr->size() * 3, 2);
	Eigen::MatrixXi ftcs_out(faces_ptr->size(), 3);

	//vert_patch_ids, verts_on_sls
	// iterate over faces
	for (int face_id = 0; face_id < faces_ptr->rows(); face_id++) {

		int face_patch_id = -1;
		for (int i = 0; i < 3; i++) {
			int vert_idx = (*faces_ptr)(face_id, i);
			// do not use vertices on nodes to find patch
			face_patch_id = vert_patch_ids[vert_idx];
			if (face_patch_id >= 0) {
				break;
			}
			//if (std::find(verts_on_sls.begin(), verts_on_sls.end(), vert_idx) == verts_on_sls.end()) {
			//	face_patch_id = vert_patch_ids[vert_idx];
			//	break;
			//}
		}
	
		//// TODO :: deal with the case with 3 verts on sls
		//if (face_patch_id == -1) {
		//	for (int i = 0; i < 3; i++) {
		//		tcs_out.row(face_id * 3 + i) << 0, 0;
		//	}
		//	continue;
		//}
		//uv_coords
		for (int i = 0; i < 3; i++) {
			Eigen::Vector2d uv(-100, -100);
			int vert_idx = (*faces_ptr)(face_id, i);
			int vert_patch_id = vert_patch_ids[vert_idx];

			// if vertex is a node
			if (vert_patch_id < 0) {
				vector<int> &patch_nodes = ms_patches->at(face_patch_id);
				int which_node = ComplexUtil::which_node(patch_nodes, vert_idx);
				switch (which_node)
				{
				case 0:
					uv << 0, 0;
					break;
				case 1:
					uv << 1, 0;
					break;
				case 2:
					uv << 1, 1;
					break;
				case 3:
					uv << 0, 1;
					break;
				}
			}
			// not a node, test if belongs to a different patch
			else if (vert_patch_id != face_patch_id) {
				// find uv with trnsfr function
				// TODO:: find a path between disconnected patches
				if (trnsfr_funcs_map->find({ vert_patch_id, face_patch_id }) == trnsfr_funcs_map->end()) {
					cout << "error" << endl;
				}
				else {
					ComplexUtil::trnsfr_func_t trnsfr_func = trnsfr_funcs_map->at({ vert_patch_id, face_patch_id });
					Eigen::Vector2d initial_uv = uv_coords->row(vert_idx);
					ComplexUtil::applyTrnsfrFunc(initial_uv, uv, trnsfr_func);
				}
			}
			// vertex belongs to the patch
			else {
				uv = uv_coords->row(vert_idx);
				//if (uv(0) == -1 && uv(1) == -1) {
				//	cout << "error! " << endl;
				//}
			}
			//if (uv(0) == -100) {
			//	cout << "error!" << endl;
			//}
			tcs_out.row(face_id * 3 + i) = uv;
		}
	}
	for (int face_id = 0; face_id < faces_ptr->rows(); face_id++) {
		ftcs_out.row(face_id) << face_id * 3, face_id * 3 + 1, face_id * 3 + 2;
	}

	igl::writeOBJ("../models/high_sphere_uv.obj", *vertices_ptr, *faces_ptr, *N, *fns_ptr, tcs_out, ftcs_out);

}