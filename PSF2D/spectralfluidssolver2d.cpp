#include "spectralfluidssolver2d.h"
#include "Spectra/MatOp/SparseGenRealShiftSolve.h"
#include "Spectra/GenEigsRealShiftSolver.h"
#include "dec.h"
#include <iostream>
#include <glm/gtx/rotate_vector.hpp>

SpectralFluidsSolver2D::SpectralFluidsSolver2D() : AbstractSolver()
{

}

void SpectralFluidsSolver2D::integrate()
{
    std::cout<<"Advect"<<std::endl;

    double e1 = 0.0f;
    double e2 = 0.0f;
    for(unsigned int k=0;k<basisCoeff.rows();k++)
    {
        e1 += basisCoeff(k,0)*basisCoeff(k,0);
    }
    Eigen::VectorXd vel(basisCoeff.rows());

    for(unsigned int k=0;k<basisCoeff.rows();k++)
    {
        vel(k) = (basisCoeff.transpose()*advection[k]*basisCoeff);
    }
    basisCoeff += timeStep*vel;
    for(unsigned int k=0;k<basisCoeff.rows();k++)
    {
        e2 += basisCoeff(k,0)*basisCoeff(k,0);
    }
    basisCoeff *= std::sqrt(e1/e2);
    for(unsigned int k=0;k<basisCoeff.rows();k++)
    {
        basisCoeff(k) *= std::exp(viscosity*-eigenValues(k)*(timeStep));
    }

    velocityField = velBasisField*basisCoeff;
    vorticityField = curl*velocityField;
    std::cout.precision(std::numeric_limits<long double>::digits10 + 1);
    std::cout<<e1<<" "<<e2<<std::endl;
    maxRotation = vorticityField.cwiseAbs().maxCoeff();
    minRotation = vorticityField.minCoeff();
}

void SpectralFluidsSolver2D::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = 1.0*derivative0(decMesh)*hodge2(decMesh,mesh->getResolution()*mesh->getResolution(),true)*derivative1(decMesh,true)*hodge1(decMesh,mesh->getResolution(),false);
    Eigen::SparseMatrix<double> bound = derivative1(decMesh);
    curl = derivative1(decMesh,true)*hodge1(decMesh,mesh->getResolution(),false);
    for(int k=0;k<bound.outerSize();k++)
    {
        unsigned int nFaces=0;
        for(Eigen::SparseMatrix<double>::InnerIterator it(bound,k);it;++it)
        {
            nFaces++;
        }
        if(nFaces!=2)
        {
            mat.prune([k](int i,int j,double v){return !(i==k||j==k);});
        }
    }

    Spectra::SparseGenRealShiftSolve<double> op(mat);
    double nearZ = 1.0f/(mesh->getResolution()*mesh->getResolution());
    Spectra::GenEigsRealShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseGenRealShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions+1,std::numeric_limits<double>::epsilon()+nearZ);
    solver.init();

    int nconv = solver.compute(1000,1e-10f,Spectra::WHICH_SM);
    eigenValues = solver.eigenvalues().real();
    velBasisField = solver.eigenvectors().real();
    vortBasisField = curl*velBasisField;
/*
    vorticityField = Eigen::VectorXd::Zero(decMesh.getNumPoints());
    for(PointIterator it=decMesh.getPointIteratorBegin();it!=decMesh.getPointIteratorEnd();it++)
    {
        unsigned int i=decMesh.getPointIndex(*it);
        glm::dvec2 point=mesh->vertex[*it].pos;
        if((point.x==256&&
            (point.y==128||point.y==384)))
        {
            vorticityField(i) = 500;
        }
    }
    setInitialVorticityField(vorticityField);*/


    velocityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();it++)
    {
        unsigned int i=decMesh.getEdgeIndex(*it);
        glm::dvec2 edge=mesh->vertex[std::get<1>(*it)].pos-mesh->vertex[std::get<0>(*it)].pos;
        if(glm::dot(glm::dvec2(1.0,0.0),edge)>std::numeric_limits<double>::epsilon()&&
           (mesh->vertex[std::get<0>(*it)].pos.x>=128||
           mesh->vertex[std::get<1>(*it)].pos.x<=256))
        {
            velocityField(i) = -(edge.x>0?1:-1)*32.0f;
        }
    }
    setInitialVelocityField(velocityField);

    buildAdvection();
}

void SpectralFluidsSolver2D::buildAdvection()
{
    std::vector<Eigen::MatrixXd> wedges;
    wedges.resize(static_cast<unsigned int>(eigenValues.rows()));
    advection.resize(static_cast<unsigned int>(eigenValues.rows()));
    for(unsigned int i=0;i<velBasisField.cols();i++)
    {
        wedges[i] = Eigen::MatrixXd(decMesh.getNumPoints(),eigenValues.rows());
        wedges[i].setZero();
        advection[i] = Eigen::MatrixXd(decMesh.getNumPoints(),eigenValues.rows());
        advection[i].setZero();
    }
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        unsigned int v1 = std::get<0>(*it);
        unsigned int v2 = std::get<1>(*it);
        unsigned int v3 = std::get<2>(*it);
        unsigned int v4 = std::get<3>(*it);

        unsigned int iv1 = decMesh.getPointIndex(v1);
        unsigned int iv2 = decMesh.getPointIndex(v2);
        unsigned int iv3 = decMesh.getPointIndex(v3);
        unsigned int iv4 = decMesh.getPointIndex(v4);

        unsigned int ie1 = decMesh.getEdgeIndex(v1,v2);
        unsigned int ie2 = decMesh.getEdgeIndex(v2,v3);
        unsigned int ie3 = decMesh.getEdgeIndex(v3,v4);
        unsigned int ie4 = decMesh.getEdgeIndex(v4,v1);

        double sig1=decMesh.getEdgeSignum(v1,v2);
        double sig2=decMesh.getEdgeSignum(v2,v3);
        double sig3=decMesh.getEdgeSignum(v3,v4);
        double sig4=decMesh.getEdgeSignum(v4,v1);

        std::tuple<unsigned int,unsigned int> edge1 = decMesh.getEdge(v1,v2);
        std::tuple<unsigned int,unsigned int> edge2 = decMesh.getEdge(v2,v3);
        std::tuple<unsigned int,unsigned int> edge3 = decMesh.getEdge(v3,v4);
        std::tuple<unsigned int,unsigned int> edge4 = decMesh.getEdge(v4,v1);

        glm::dvec2 e1 = (sig1*(mesh->vertex[std::get<1>(edge1)].pos-mesh->vertex[std::get<0>(edge1)].pos));
        glm::dvec2 e2 = (sig2*(mesh->vertex[std::get<1>(edge2)].pos-mesh->vertex[std::get<0>(edge2)].pos));
        glm::dvec2 e3 = (sig3*(mesh->vertex[std::get<1>(edge3)].pos-mesh->vertex[std::get<0>(edge3)].pos));
        glm::dvec2 e4 = (sig4*(mesh->vertex[std::get<1>(edge4)].pos-mesh->vertex[std::get<0>(edge4)].pos));

        glm::dvec2 n1 = glm::rotate(glm::normalize(e1),glm::radians(-90.0));
        glm::dvec2 n2 = glm::rotate(glm::normalize(e2),glm::radians(-90.0));
        glm::dvec2 n3 = glm::rotate(glm::normalize(e3),glm::radians(-90.0));
        glm::dvec2 n4 = glm::rotate(glm::normalize(e4),glm::radians(-90.0));


        for(unsigned int i=0;i<nEigenFunctions;i++)
        {
            double vel1a = sig1*eigenFunctions[i](ie1);
            double vel2a = sig2*eigenFunctions[i](ie2);
            double vel3a = sig3*eigenFunctions[i](ie3);
            double vel4a = sig4*eigenFunctions[i](ie4);

            for(unsigned int j=0;j<nEigenFunctions;j++)
            {
                double vel1b = sig1*eigenFunctions[j](ie1);
                double vel2b = sig2*eigenFunctions[j](ie2);
                double vel3b = sig3*eigenFunctions[j](ie3);
                double vel4b = sig4*eigenFunctions[j](ie4);

                wedges[i](iv2,j) += (0.25f*mesh->getResolution()*mesh->getResolution())*(vel1a*vel2b-vel1b*vel2a)*(n1.x*n2.y-n1.y*n2.x);
                wedges[i](iv3,j) += (0.25f*mesh->getResolution()*mesh->getResolution())*(vel2a*vel3b-vel2b*vel3a)*(n2.x*n3.y-n2.y*n3.x);
                wedges[i](iv4,j) += (0.25f*mesh->getResolution()*mesh->getResolution())*(vel3a*vel4b-vel3b*vel4a)*(n3.x*n4.y-n3.y*n4.x);
                wedges[i](iv1,j) += (0.25f*mesh->getResolution()*mesh->getResolution())*(vel4a*vel1b-vel4b*vel1a)*(n4.x*n1.y-n4.y*n1.x);
            }
        }
    }
/*
    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            advection[i].col(j) = eigenValues(i)*wedges[i].col(j);
        }
    }*/
    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        advection[i] = (eigenValues(i)*vortBasisField.transpose()*wedges[i]);
    }
}
