#include "abstractsolver.h"
#include <iostream>

AbstractSolver::AbstractSolver()
{
    mesh = NULL;
    nEigenFunctions = 32;
    viscosity = 0.0f;
    timeStep = 1.0f/60.0f;
}

void AbstractSolver::setMesh(Mesh2D* mesh)
{
    if(mesh!=NULL)
    {
        delete this->mesh;
    }
    this->mesh = mesh;
    decMesh = mesh->voxelize();
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setInitialVelocityField(const Eigen::VectorXf& field)
{
    basisCoeff = velBasisField.transpose()*field;
}

void AbstractSolver::setInitialVorticityField(const Eigen::VectorXf& field)
{
    basisCoeff = (vortBasisField.transpose()*field);
}

void AbstractSolver::setNumberEigenFunctions(unsigned int n)
{
    this->nEigenFunctions = n;
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setResolution(unsigned int res)
{
    mesh->setResolution(res);
    decMesh = mesh->voxelize();
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setTimestep(float timestep)
{
    this->timeStep = timestep;
}

void AbstractSolver::setViscosity(float visc)
{
    this->viscosity = visc;
}

Mesh2D* AbstractSolver::getMesh()
{
    return mesh;
}

DECMesh2D& AbstractSolver::getDECMesh()
{
    return decMesh;
}

unsigned int AbstractSolver::getNumEigenFunctions()
{
    return nEigenFunctions;
}

float AbstractSolver::getTimestep()
{
    return timeStep;
}

const Eigen::VectorXf& AbstractSolver::getEigenFunction(unsigned int n)
{
    return eigenFunctions[n];
}

const Eigen::MatrixXf& AbstractSolver::getVelocityBasisField()
{
    return velBasisField;
}

const Eigen::MatrixXf& AbstractSolver::getVorticityBasisField()
{
    return vortBasisField;
}

const Eigen::VectorXf& AbstractSolver::getBasisCoefficients()
{
    return basisCoeff;
}

const Eigen::VectorXf& AbstractSolver::getVelocityField()
{
    return velocityField;
}

const Eigen::VectorXf& AbstractSolver::getVorticityField()
{
    return vorticityField;
}

float AbstractSolver::getMaxVorticity()
{
    return maxRotation;
}

float AbstractSolver::getMinVorticity()
{
    return minRotation;
}

void AbstractSolver::buildEigenFunctions()
{
    eigenFunctions.resize(nEigenFunctions);
    for(unsigned int i=0;i<velBasisField.cols();i++)
    {
        Eigen::VectorXf e = velBasisField.col(i);
        std::cout<<velBasisField.cols()<<" "<<e.rows()<<std::endl;
        eigenFunctions[i] = e;
    }
}
