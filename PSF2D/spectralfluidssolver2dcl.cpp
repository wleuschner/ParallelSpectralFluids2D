#include "spectralfluidssolver2dcl.h"
#include "Spectra/MatOp/SparseSymShiftSolve.h"
#include "Spectra/SymEigsShiftSolver.h"
#include "dec.h"
#include <iostream>
#include <glm/gtx/rotate_vector.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/linalg/inner_prod.hpp>

SpectralFluidsSolver2DCL::SpectralFluidsSolver2DCL() : AbstractSolver()
{

}

void SpectralFluidsSolver2DCL::integrate()
{
    viennacl::scalar<double> e1 = 0.0f;
    viennacl::scalar<double> e2 = 0.0f;

    e1 = viennacl::linalg::inner_prod(vclBasisCoeff,vclBasisCoeff);

    viennacl::vector<double> vel(nEigenFunctions);
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        viennacl::scalar<double> a = viennacl::linalg::inner_prod(vclBasisCoeff,viennacl::linalg::prod(vclAdvection[k],vclBasisCoeff));
        vel(k) = a;
    }
    vclBasisCoeff += timeStep*vel;

    e2 = viennacl::linalg::inner_prod(vclBasisCoeff,vclBasisCoeff);

    vclBasisCoeff *= std::sqrt(e1/e2);

    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        vclBasisCoeff(k) *= std::exp(viscosity*eigenValues(k)*(timeStep));
    }

    vclVelocityField = viennacl::linalg::prod(vclVelBasisField,vclBasisCoeff);
    vclVorticityField = viennacl::linalg::prod(vclCurl,vclVelocityField);
    maxRotation = vorticityField.cwiseAbs().maxCoeff();
    minRotation = vorticityField.minCoeff();
}

void SpectralFluidsSolver2DCL::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = -1.0f*(derivative0(decMesh)*hodge2(decMesh,mesh->getResolution()*mesh->getResolution(),true)*derivative1(decMesh,true)*hodge1(decMesh,mesh->getResolution(),false));
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

    bool decompositionDone=false;
    double nearZ = 1.0f/(mesh->getResolution()*mesh->getResolution());
    double omega = 0.0f;
    while(!decompositionDone)
    {
        try
        {
            Spectra::SparseSymShiftSolve<double> op(mat);
            Spectra::SymEigsShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseSymShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions,omega);
            solver.init();

            int nconv = solver.compute(1000,1e-10f,Spectra::WHICH_LM);
            if(solver.info()==Spectra::SUCCESSFUL)
            {
                eigenValues = solver.eigenvalues().real();
                velBasisField = solver.eigenvectors().real();
                vortBasisField = (curl*velBasisField);
                decompositionDone = true;
            }
        }
        catch(std::runtime_error e)
        {
            omega+=nearZ;
            std::cout<<"Increase omega"<<std::endl;
        }
    }

    vorticityField = Eigen::VectorXd::Zero(decMesh.getNumPoints());
    for(PointIterator it=decMesh.getPointIteratorBegin();it!=decMesh.getPointIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            glm::dvec2 point=mesh->vertex[it->id].pos;
            if((point.x==512&&
                (point.y==448||point.y==576)))
            {
                vorticityField(it->id) = 1000*2*3.141;
            }
        }
    }
    setInitialVorticityField(vorticityField);

/*
    velocityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();it++)
    {
        unsigned int i=decMesh.getEdgeIndex(*it);
        glm::dvec2 edge=mesh->vertex[std::get<1>(*it)].pos-mesh->vertex[std::get<0>(*it)].pos;
        if(glm::dot(glm::dvec2(1.0,0.0),edge)>std::numeric_limits<double>::epsilon())
        {
            velocityField(i) = -(edge.x>0?1:-1)*16.0f;
        }
    }
    setInitialVelocityField(velocityField);*/

}

void SpectralFluidsSolver2DCL::buildAdvection()
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
    for(FaceIterator fit=decMesh.getFaceIteratorBegin();fit!=decMesh.getFaceIteratorEnd();fit++)
    {
        if(fit->inside==GridState::INSIDE)
        {
            unsigned int v1 = fit->v1;
            unsigned int v2 = fit->v2;
            unsigned int v3 = fit->v3;
            unsigned int v4 = fit->v4;

            unsigned int iv1 = v1;
            unsigned int iv2 = v2;
            unsigned int iv3 = v3;
            unsigned int iv4 = v4;

            unsigned int ie1 = fit->e1;
            unsigned int ie2 = fit->e2;
            unsigned int ie3 = fit->e3;
            unsigned int ie4 = fit->e4;

            double sig1=decMesh.getEdgeSignum(fit->e1,v1,v2);
            double sig2=decMesh.getEdgeSignum(fit->e2,v2,v3);
            double sig3=decMesh.getEdgeSignum(fit->e3,v3,v4);
            double sig4=decMesh.getEdgeSignum(fit->e4,v4,v1);

            Edge2D edge1 = decMesh.getEdge(ie1);
            Edge2D edge2 = decMesh.getEdge(ie2);
            Edge2D edge3 = decMesh.getEdge(ie3);
            Edge2D edge4 = decMesh.getEdge(ie4);

            glm::dvec2 e1 = (sig1*(mesh->vertex[edge1.v2].pos-mesh->vertex[edge1.v1].pos));
            glm::dvec2 e2 = (sig2*(mesh->vertex[edge2.v2].pos-mesh->vertex[edge2.v1].pos));
            glm::dvec2 e3 = (sig3*(mesh->vertex[edge3.v2].pos-mesh->vertex[edge3.v1].pos));
            glm::dvec2 e4 = (sig4*(mesh->vertex[edge4.v2].pos-mesh->vertex[edge4.v1].pos));

            glm::dvec2 n1 = glm::rotate(glm::normalize(e1),glm::radians(-90.0));
            glm::dvec2 n2 = glm::rotate(glm::normalize(e2),glm::radians(-90.0));
            glm::dvec2 n3 = glm::rotate(glm::normalize(e3),glm::radians(-90.0));
            glm::dvec2 n4 = glm::rotate(glm::normalize(e4),glm::radians(-90.0));


            for(unsigned int i=0;i<nEigenFunctions;i++)
            {
                double vel1a = sig1*velBasisField(ie1,i);
                double vel2a = sig2*velBasisField(ie2,i);
                double vel3a = sig3*velBasisField(ie3,i);
                double vel4a = sig4*velBasisField(ie4,i);

                for(unsigned int j=0;j<nEigenFunctions;j++)
                {
                    double vel1b = sig1*velBasisField(ie1,j);
                    double vel2b = sig2*velBasisField(ie2,j);
                    double vel3b = sig3*velBasisField(ie3,j);
                    double vel4b = sig4*velBasisField(ie4,j);

                    wedges[i](iv2,j) += (0.25)*(vel1a*vel2b-vel1b*vel2a)*(n1.x*n2.y-n1.y*n2.x);
                    wedges[i](iv3,j) += (0.25)*(vel2a*vel3b-vel2b*vel3a)*(n2.x*n3.y-n2.y*n3.x);
                    wedges[i](iv4,j) += (0.25)*(vel3a*vel4b-vel3b*vel4a)*(n3.x*n4.y-n3.y*n4.x);
                    wedges[i](iv1,j) += (0.25)*(vel4a*vel1b-vel4b*vel1a)*(n4.x*n1.y-n4.y*n1.x);
                }
            }
        }
    }

    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            advection[i].col(j) = eigenValues(i)*wedges[i].col(j);
        }
    }
    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        advection[i] = vortBasisField.transpose()*advection[i];
    }
}
