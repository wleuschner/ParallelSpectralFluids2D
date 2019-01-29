#include "spectralfluidssolver2domp.h"
#include "Spectra/MatOp/SparseSymShiftSolve.h"
#include "Spectra/SymEigsShiftSolver.h"
#include "dec.h"
#include <iostream>
#include <glm/gtx/rotate_vector.hpp>

SpectralFluidsSolver2DOMP::SpectralFluidsSolver2DOMP() : AbstractSolver()
{

}

void SpectralFluidsSolver2DOMP::integrate()
{
    double e1 = 0.0;
    double e2 = 0.0;

    e1 = basisCoeff.dot(basisCoeff);

    Eigen::VectorXd vel(nEigenFunctions);
    #pragma omp parallel for
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        vel(k) = (basisCoeff.transpose()*advection[k]*basisCoeff);
    }
    basisCoeff += timeStep*vel;

    e2 = basisCoeff.dot(basisCoeff);

    basisCoeff *= std::sqrt(e1/e2);

    #pragma omp parallel for
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        basisCoeff(k) *= std::exp(-viscosity*eigenValues(k)*(timeStep));
    }

    velocityField = velBasisField*basisCoeff;

    //#pragma omp parallel for
    for(std::vector<glm::dvec2>::iterator it=particles.begin();it<particles.end();it++)
    {
        /*unsigned int yOfsMinus = static_cast<unsigned int>(it->y-0.5*resolution)/resolution;
        unsigned int yOfsPlus = static_cast<unsigned int>(it->y+0.5*resolution)/resolution;
        unsigned int xOfsMinus = static_cast<unsigned int>(it->x-0.5*resolution)/resolution;
        unsigned int xOfsPlus = static_cast<unsigned int>(it->x+0.5*resolution)/resolution;

        Face2D f1 = decMesh.getFace(yOfsMinus*((mesh->getWidth()/resolution)+2)+xOfsMinus);
        Face2D f2 = decMesh.getFace(yOfsPlus*((mesh->getWidth()/resolution)+2)+xOfsMinus);
        Face2D f3 = decMesh.getFace(yOfsMinus*((mesh->getWidth()/resolution)+2)+xOfsPlus);
        Face2D f4 = decMesh.getFace(yOfsPlus*((mesh->getWidth()/resolution)+2)+xOfsPlus);


        glm::dvec2 cxf1 = mesh->vertex[f1.v1].pos+0.5*(mesh->vertex[f1.v2].pos-mesh->vertex[f1.v1].pos);
        glm::dvec2 cxf2 = mesh->vertex[f2.v1].pos+0.5*(mesh->vertex[f2.v2].pos-mesh->vertex[f2.v1].pos);
        glm::dvec2 cxf3 = mesh->vertex[f3.v1].pos+0.5*(mesh->vertex[f3.v2].pos-mesh->vertex[f3.v1].pos);
        glm::dvec2 cxf4 = mesh->vertex[f4.v1].pos+0.5*(mesh->vertex[f4.v2].pos-mesh->vertex[f4.v1].pos);

        double vel1x = decMesh.getEdgeSignum(f1.e1,f1.v1,f1.v2)*velocityField(f1.e1);
        double vel2x = decMesh.getEdgeSignum(f2.e1,f2.v1,f2.v2)*velocityField(f2.e1);
        double vel3x = decMesh.getEdgeSignum(f3.e1,f3.v1,f3.v2)*velocityField(f3.e1);
        double vel4x = decMesh.getEdgeSignum(f4.e1,f4.v1,f4.v2)*velocityField(f4.e1);

        glm::dvec2 cyf1 = mesh->vertex[f1.v4].pos+0.5*(mesh->vertex[f1.v1].pos-mesh->vertex[f1.v4].pos);
        glm::dvec2 cyf2 = mesh->vertex[f2.v4].pos+0.5*(mesh->vertex[f2.v1].pos-mesh->vertex[f2.v4].pos);
        glm::dvec2 cyf3 = mesh->vertex[f3.v4].pos+0.5*(mesh->vertex[f3.v1].pos-mesh->vertex[f3.v4].pos);
        glm::dvec2 cyf4 = mesh->vertex[f4.v4].pos+0.5*(mesh->vertex[f4.v1].pos-mesh->vertex[f4.v4].pos);

        double vel1y = decMesh.getEdgeSignum(f1.e4,f1.v4,f1.v1)*velocityField(f1.e4);
        double vel2y = decMesh.getEdgeSignum(f2.e4,f2.v4,f2.v1)*velocityField(f2.e4);
        double vel3y = decMesh.getEdgeSignum(f3.e4,f3.v4,f3.v1)*velocityField(f3.e4);
        double vel4y = decMesh.getEdgeSignum(f4.e4,f4.v4,f4.v1)*velocityField(f4.e4);

        glm::dvec2 particleNormalizedX = (1.0/resolution)*((*it)-(cxf1));
        glm::dvec2 particleNormalizedY = (1.0/resolution)*((*it)-(cyf1));

        glm::dvec2 yVelXInterp = glm::mix(glm::dvec2(vel1x,vel3x),glm::dvec2(vel2x,vel4x),particleNormalizedX.y);
        double xVel = glm::mix(yVelXInterp.x,yVelXInterp.y,particleNormalizedX.x);

        glm::dvec2 yVelYInterp = glm::mix(glm::dvec2(vel1y,vel3y),glm::dvec2(vel2y,vel4y),particleNormalizedY.y);
        double yVel = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.x);

        glm::dvec2 vel = glm::dvec2(xVel,yVel);
        //std::cout<<"PARTICLE: "<<it->x<<" "<<it->y<<std::endl;
        //std::cout<<"ORIGIN: "<<cxf1.x<<" "<<cxf1.y<<std::endl;
        //std::cout<<"X:"<<particleNormalizedX.x<<" "<<particleNormalizedX.y<<std::endl;
        //std::cout<<"Y:"<<particleNormalizedY.x<<" "<<particleNormalizedY.y<<std::endl;
        (*it) = (*it)+timeStep*vel;*/


        unsigned int yOfs = static_cast<unsigned int>(it->y)/resolution;
        unsigned int xOfs = static_cast<unsigned int>(it->x)/resolution;

        unsigned int yOfsMinus = static_cast<unsigned int>(it->y-resolution/2)/resolution;
        unsigned int yOfsPlus = static_cast<unsigned int>(it->y+resolution/2)/resolution;
        unsigned int xOfsMinus = static_cast<unsigned int>(it->x-resolution/2)/resolution;
        unsigned int xOfsPlus = static_cast<unsigned int>(it->x+resolution/2)/resolution;

        Face2D f1x = decMesh.getFace(yOfsMinus*((mesh->getWidth()/resolution)+2)+xOfs);
        Face2D f2x = decMesh.getFace(yOfsPlus*((mesh->getWidth()/resolution)+2)+xOfs);
        Face2D f1y = decMesh.getFace(yOfs*((mesh->getWidth()/resolution)+2)+xOfsMinus);
        Face2D f2y = decMesh.getFace(yOfs*((mesh->getWidth()/resolution)+2)+xOfsPlus);


        glm::dvec2 cf1x = mesh->vertex[f1x.v1].pos+0.5*(mesh->vertex[f1x.v2].pos-mesh->vertex[f1x.v1].pos);
        glm::dvec2 cf2x = mesh->vertex[f2x.v1].pos+0.5*(mesh->vertex[f2x.v2].pos-mesh->vertex[f2x.v1].pos);
        glm::dvec2 cf1y = mesh->vertex[f1y.v4].pos+0.5*(mesh->vertex[f1y.v1].pos-mesh->vertex[f1y.v4].pos);
        glm::dvec2 cf2y = mesh->vertex[f2y.v4].pos+0.5*(mesh->vertex[f2y.v1].pos-mesh->vertex[f2y.v4].pos);

        double vel1x = decMesh.getEdgeSignum(f1x.e1,f1x.v1,f1x.v2)*velocityField(f1x.e1);
        double vel2x = decMesh.getEdgeSignum(f2x.e1,f2x.v1,f2x.v2)*velocityField(f2x.e1);
        double vel3x = -decMesh.getEdgeSignum(f1x.e3,f1x.v3,f1x.v4)*velocityField(f1x.e3);
        double vel4x = -decMesh.getEdgeSignum(f2x.e3,f2x.v3,f2x.v4)*velocityField(f2x.e3);

        double vel1y = decMesh.getEdgeSignum(f1y.e4,f1y.v4,f1y.v1)*velocityField(f1y.e4);
        double vel3y = decMesh.getEdgeSignum(f2y.e4,f2y.v4,f2y.v1)*velocityField(f2y.e4);
        double vel2y = -decMesh.getEdgeSignum(f1y.e2,f1y.v2,f1y.v3)*velocityField(f1y.e2);
        double vel4y = -decMesh.getEdgeSignum(f2y.e2,f2y.v2,f2y.v3)*velocityField(f2y.e2);

        glm::dvec2 particleNormalizedX;
        glm::dvec2 particleNormalizedY;
        glm::dvec2 vel;

        if(cf1x.y<cf2x.y)
        {
            particleNormalizedX = (1.0/resolution)*((*it)-(cf1x));
            glm::dvec2 yVelXInterp = glm::mix(glm::dvec2(vel1x,vel3x),glm::dvec2(vel2x,vel4x),particleNormalizedX.y);
            vel.x = glm::mix(yVelXInterp.x,yVelXInterp.y,particleNormalizedX.x);
        }
        else
        {
            particleNormalizedX = (1.0/resolution)*((*it)-(cf2x));
            glm::dvec2 yVelXInterp = glm::mix(glm::dvec2(vel2x,vel4x),glm::dvec2(vel1x,vel3x),particleNormalizedX.y);
            vel.x = glm::mix(yVelXInterp.x,yVelXInterp.y,particleNormalizedX.x);
        }

        if(cf1y.x<cf2y.x)
        {
            particleNormalizedY = (1.0/resolution)*((*it)-(cf1y));
            glm::dvec2 yVelYInterp = glm::mix(glm::dvec2(vel1y,vel3y),glm::dvec2(vel2y,vel4y),particleNormalizedY.y);
            vel.y = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.x);
        }
        else
        {
            particleNormalizedY = (1.0/resolution)*((*it)-(cf2y));
            glm::dvec2 yVelYInterp = glm::mix(glm::dvec2(vel2y,vel4y),glm::dvec2(vel1y,vel3y),particleNormalizedY.y);
            vel.y = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.x);
        }


        /*std::cout<<"PARTICLE: "<<it->x<<" "<<it->y<<std::endl;
        std::cout<<"ORIGIN: "<<cf1x.x<<" "<<cf1x.y<<std::endl;
        std::cout<<"X:"<<particleNormalizedX.x<<" "<<particleNormalizedX.y<<std::endl;
        std::cout<<"Y:"<<particleNormalizedY.x<<" "<<particleNormalizedY.y<<std::endl;*/
        (*it) = (*it)+timeStep*vel;


        /*unsigned int yOfs = static_cast<unsigned int>(it->y)/resolution;
        unsigned int xOfs = static_cast<unsigned int>(it->x)/resolution;
        Face2D face = decMesh.getFace(yOfs*((mesh->getWidth()/resolution)+2)+xOfs);

        if(face.inside==GridState::INSIDE)
        {
            unsigned int iv1 = face.v1;
            unsigned int iv2 = face.v2;
            unsigned int iv3 = face.v3;
            unsigned int iv4 = face.v4;

            unsigned int ie1 = face.e1;
            unsigned int ie2 = face.e2;
            unsigned int ie3 = face.e3;
            unsigned int ie4 = face.e4;

            Vertex2D v1 = getMesh()->vertex[iv1];
            Vertex2D v2 = getMesh()->vertex[iv2];
            Vertex2D v3 = getMesh()->vertex[iv3];
            Vertex2D v4 = getMesh()->vertex[iv4];

            glm::dvec2 e1 = v2.pos-v1.pos;
            glm::dvec2 e2 = v3.pos-v2.pos;
            glm::dvec2 e3 = v4.pos-v3.pos;
            glm::dvec2 e4 = v1.pos-v4.pos;

            unsigned int fidx = face.id;

            double s1=decMesh.getEdgeSignum(face.e1,iv1,iv2);
            double s2=decMesh.getEdgeSignum(face.e2,iv2,iv3);
            double s3=decMesh.getEdgeSignum(face.e3,iv3,iv4);
            double s4=decMesh.getEdgeSignum(face.e4,iv4,iv1);

            glm::dvec2 vel = glm::rotate(0.5*(s1*velocityField(ie1)*glm::normalize(e1)+s3*velocityField(ie3)*glm::normalize(e3))+
                                        0.5*(s2*velocityField(ie2)*glm::normalize(e2)+s4*velocityField(ie4)*glm::normalize(e4)),
                                        glm::radians(-90.0));

            (*it) = (*it)+timeStep*vel;
        }*/
    }
    vorticityField = curl*velocityField;
    maxRotation = vorticityField.cwiseAbs().maxCoeff();
    minRotation = vorticityField.minCoeff();
}

void SpectralFluidsSolver2DOMP::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = -1.0*(derivative0(decMesh)*hodge2(decMesh,1.0,true)*derivative1(decMesh,true)*hodge1(decMesh,1.0,false));
    Eigen::SparseMatrix<double> bound = derivative1(decMesh);
    curl = derivative1(decMesh,true)*hodge1(decMesh,1.0,false);
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
    mat.pruned();

    bool decompositionDone=false;
    double omega = -0.1;
    while(!decompositionDone)
    {
        try
        {
            Spectra::SparseSymShiftSolve<double> op(mat);
            Spectra::SymEigsShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseSymShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions,omega);
            solver.init();

            int nconv = solver.compute(1000,1e-1,Spectra::WHICH_SM);
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
            omega+=std::numeric_limits<double>::epsilon();
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
                (point.y==512||point.y==512)))
            {
                vorticityField(it->id) = 1000*2*3.141;
            }
        }
    }
    setInitialVorticityField(vorticityField);


    /*velocityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            glm::dvec2 edge=mesh->vertex[it->v2].pos-mesh->vertex[it->v1].pos;
            if(glm::dot(glm::dvec2(1.0,0.0),edge)>std::numeric_limits<double>::epsilon())
            {
                velocityField(it->id) = (edge.x>0?1:-1)*8.0;
            }
        }
    }
    setInitialVelocityField(velocityField);*/

}

void SpectralFluidsSolver2DOMP::buildAdvection()
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
