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

int EIG_NUM = 1;


int main(int argc, char *argv[])
{

	auto vertices_ptr = std::make_shared<Eigen::MatrixXd>();
	auto faces_ptr = std::make_shared<Eigen::MatrixXi>();

	igl::readOBJ("../models/sphere1.obj", *vertices_ptr, *faces_ptr);


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

	U = ((U.array() - U.minCoeff()) / (U.maxCoeff() - U.minCoeff())).eval();
	std::cout << "U: " << U << endl;

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

	// iterate over steep lines
	std::ofstream lines_out("../line_verts.txt");
	int lines_cnt = 0;
	lines_out << "verts = [None] * " + std::to_string(steeplines->size());
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
}

	