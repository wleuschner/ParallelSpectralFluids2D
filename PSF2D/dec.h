#ifndef DEC_H
#define DEC_H
#include <Eigen/Eigen>
#include "mesh2d.h"
#include "decmesh2d.h"

Eigen::SparseMatrix<float> hodge2(DECMesh2D& mesh,float area,bool dual=false);
Eigen::SparseMatrix<float> hodge1(DECMesh2D& mesh,float area,bool dual=false);
Eigen::SparseMatrix<float> hodge0(DECMesh2D& mesh,float area,bool dual=false);

Eigen::SparseMatrix<float> derivative0(DECMesh2D& mesh,bool dual=false);
Eigen::SparseMatrix<float> derivative1(DECMesh2D& mesh,bool dual=false);




#endif // DEC_H
