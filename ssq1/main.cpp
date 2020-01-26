#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseExtra>

#include <string.h>

#include "HalfEdge.h"
#include "complexUtil.h"

using namespace std;

int EIG_NUM = 10;


int main(int argc, char *argv[])
{

	//auto vertices_ptr = std::make_shared<Eigen::MatrixXd>();
	//auto faces_ptr = std::make_shared<Eigen::MatrixXi>();

	//igl::readOBJ("../models/sphere3.obj", *vertices_ptr, *faces_ptr);

	Eigen::MatrixXd tcs;
	Eigen::MatrixXi ftcs;

	auto vertices_ptr = std::make_shared<Eigen::MatrixXd>();
	auto vns_ptr = std::make_shared<Eigen::MatrixXd>();
	auto faces_ptr = std::make_shared<Eigen::MatrixXi>();
	auto fns_ptr = std::make_shared<Eigen::MatrixXi>();

	// //Load a mesh in OFF format
	igl::readOBJ("../models/sphere3.obj", *vertices_ptr, tcs, *vns_ptr, *faces_ptr, ftcs, *fns_ptr);


	Eigen::SparseMatrix<double> L, M;
	igl::cotmatrix(*vertices_ptr, *faces_ptr, L);
	L = (-L).eval();

	igl::massmatrix(*vertices_ptr, *faces_ptr, igl::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::MatrixXd U;
	Eigen::VectorXd D;
	//igl::eigs(L, M, 5, igl::EIGS_TYPE_SM, U, D);
	if (!igl::eigs(L, M, EIG_NUM, igl::EIGS_TYPE_SM, U, D))
	{
		cout << "failed." << endl;
		exit(1);
	}

	std::ofstream u_file("../U.txt");
	U = ((U.array() - U.minCoeff()) / (U.maxCoeff() - U.minCoeff())).eval();
	for (int i = 0; i < U.rows(); i++) {
		u_file << i << ": (" << vertices_ptr->row(i) << "), U: " << U(i) << endl;
	}
	std::cout << "U: " << U << endl;
	u_file << U << endl;
	u_file.close();

	std::shared_ptr< std::vector<HalfEdge::HE>> HE_edges = nullptr;
	std::shared_ptr<std::vector<int>> HE_verts = nullptr, HE_faces = nullptr;
	HalfEdge::createHEs(*faces_ptr, vertices_ptr->rows(), &HE_edges, &HE_verts, &HE_faces);

	HalfEdge::vert_t saddles = nullptr;
	HalfEdge::vert_t maxs = nullptr;
	HalfEdge::vert_t mins = nullptr;
	ComplexUtil::vert_type_t vert_types = nullptr;
	ComplexUtil::findVertexTypes(&saddles, &maxs, &mins, &vert_types, HE_verts, HE_edges, U);

	ComplexUtil::steep_lines_t steeplines;
	steeplines = ComplexUtil::findSteepLines(saddles, vert_types, HE_verts, HE_edges, U);

	std::vector<std::vector<int>> ms_patches = ComplexUtil::buildMsPatches(steeplines, vertices_ptr, vns_ptr);

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
		lines_out << + "]\n";
		lines_cnt++;
	}
	lines_out.close();

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

	//// write out vertices faces
	//std::ofstream faces_out("../vert_faces.txt");
	//faces_out << "faces = [None] * " + std::to_string(vertices_ptr->rows()) << endl;
	//for (int v_indx = 0; v_indx < vertices_ptr->rows(); v_indx++) {
	//	std::shared_ptr<vector<int>> faces = HalfEdge::getFacesOfV(v_indx, HE_verts, HE_edges);
	//
	//	faces_out << "faces[" << v_indx << "] = [";
	//	for (auto f_it = faces->begin(); f_it != faces->end(); ++f_it) {
	//		faces_out << *f_it << ", ";
	//	}
	//	faces_out << "];" << endl;
	//}
	//faces_out.close();

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

	