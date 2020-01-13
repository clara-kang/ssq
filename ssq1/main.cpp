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
	ComplexUtil::vert_type_t vert_types = nullptr;
	ComplexUtil::findVertexTypes(&saddles, &vert_types, HE_verts, HE_edges, U);
	
	ComplexUtil::steep_lines_t steeplines;
	steeplines = ComplexUtil::findSteepLines(saddles, vert_types, HE_verts, HE_edges, U);
	
	// write all steeplines to obj files
	int file_cnt = 0;
	for (auto it = steeplines->begin(); it != steeplines->end(); ++it) {
		vector<int> steepline = it->second;
		Eigen::MatrixXd V(steepline.size(), 3);
		Eigen::MatrixXi F(steepline.size() - 2, 3);
		int row_id = 0;
		for (auto sl_it = steepline.begin(); sl_it != steepline.end(); ++sl_it) {
			V.row(row_id) = vertices_ptr->row(*sl_it);
			row_id++;
		}
		for (int row_id = 0; row_id < steepline.size() - 2; row_id++) {
			F.row(row_id) << row_id, row_id+1, row_id+2;
		}
		string file_name = "../models/lines/strand" + std::to_string(file_cnt) + ".obj";
		igl::writeOBJ(file_name, V, F);
		file_cnt++;
	}

	//Eigen::MatrixXd V(3, 3);
	//Eigen::MatrixXi F(1, 3);
	//V.row(0) << 0.0, 0.0, 0.0;
	//V.row(1) << 0.0, 0.2, 0.0;
	//V.row(2) << 0.0, 0.4, 0.0;

	//F.row(0) << 0, 1, 2;
	//igl::writeOBJ("strand.obj", V, F);
}