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
    resolution = mesh->getResolution();
    this->mesh = mesh;
    decMesh = mesh->voxelize();
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setInitialVelocityField(const Eigen::VectorXd& field)
{
    basisCoeff = velBasisField.transpose()*field;
}

void AbstractSolver::setInitialVorticityField(const Eigen::VectorXd& field)
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
    this->resolution = res;
    mesh->setResolution(res);
    decMesh = mesh->voxelize();
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setTimestep(double timestep)
{
    this->timeStep = timestep;
}

void AbstractSolver::setViscosity(double visc)
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

double AbstractSolver::getTimestep()
{
    return timeStep;
}

const Eigen::VectorXd& AbstractSolver::getEigenFunction(unsigned int n)
{
    return eigenFunctions[n];
}

const Eigen::MatrixXd& AbstractSolver::getVelocityBasisField()
{
    return velBasisField;
}

const Eigen::MatrixXd& AbstractSolver::getVorticityBasisField()
{
    return vortBasisField;
}

const Eigen::VectorXd& AbstractSolver::getBasisCoefficients()
{
    return basisCoeff;
}

const Eigen::VectorXd& AbstractSolver::getVelocityField()
{
    return velocityField;
}

const Eigen::VectorXd& AbstractSolver::getVorticityField()
{
    return vorticityField;
}

double AbstractSolver::getMaxVorticity()
{
    return maxRotation;
}

double AbstractSolver::getMinVorticity()
{
    return minRotation;
}

void AbstractSolver::buildEigenFunctions()
{
    eigenFunctions.resize(nEigenFunctions);
    for(unsigned int i=0;i<velBasisField.cols();i++)
    {
        Eigen::VectorXd e = velBasisField.col(i);
        eigenFunctions[i] = e;
    }
}

const std::vector<glm::dvec2>& AbstractSolver::getParticles()
{
    return particles;
}

void AbstractSolver::clearParticles()
{
    particles.clear();
}

void AbstractSolver::addParticle(glm::dvec2 particle)
{
    if(decMesh.isPointInside(particle))
    {
        particles.push_back(particle);
    }
}

unsigned int AbstractSolver::getNumParticles()
{
    return particles.size();
}
