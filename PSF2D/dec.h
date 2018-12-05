#ifndef DEC_H
#define DEC_H
#include <Eigen/Eigen>
#include "mesh2d.h"

Eigen::SparseMatrix<float> hodge2(Mesh2D& mesh,bool inverse=false);
Eigen::SparseMatrix<float> hodge1(Mesh2D& mesh,bool inverse=false);
Eigen::SparseMatrix<float> hodge0(Mesh2D& mesh);

Eigen::SparseMatrix<float> derivative0(Mesh2D& mesh,bool dual=false);
Eigen::SparseMatrix<float> derivative1(Mesh2D& mesh,bool dual=false);




#endif // DEC_H
