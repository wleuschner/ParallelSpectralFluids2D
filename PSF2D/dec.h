#ifndef DEC_H
#define DEC_H
#include <Eigen/Eigen>
#include "mesh2d.h"
#include "decmesh2d.h"

Eigen::SparseMatrix<double> hodge2(DECMesh2D& mesh,double area,bool dual=false);
Eigen::SparseMatrix<double> hodge1(DECMesh2D& mesh,double area,bool dual=false);
Eigen::SparseMatrix<double> hodge0(DECMesh2D& mesh,double area,bool dual=false);

Eigen::SparseMatrix<double> derivative0(DECMesh2D& mesh,bool dual=false);
Eigen::SparseMatrix<double> derivative1(DECMesh2D& mesh,bool dual=false);




#endif // DEC_H
