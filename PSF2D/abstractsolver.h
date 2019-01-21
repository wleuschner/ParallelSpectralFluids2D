#ifndef ABSTRACTSOLVER_H
#define ABSTRACTSOLVER_H
#include<Eigen/Eigen>
#include<vector>
#include"mesh2d.h"
#include"decmesh2d.h"

class AbstractSolver
{
public:
    AbstractSolver();

    virtual void integrate()=0;
    void setMesh(Mesh2D* mesh);

    void setInitialVelocityField(const Eigen::VectorXf& field);
    void setInitialVorticityField(const Eigen::VectorXf& field);

    void setNumberEigenFunctions(unsigned int n);
    void setResolution(unsigned int res);
    void setViscosity(float visc);
    void setTimestep(float timestep);

    Mesh2D* getMesh();
    DECMesh2D& getDECMesh();

    unsigned int getNumEigenFunctions();
    float getTimestep();
    const Eigen::VectorXf& getEigenFunction(unsigned int n);
    const Eigen::MatrixXf& getVelocityBasisField();
    const Eigen::MatrixXf& getVorticityBasisField();
    const Eigen::VectorXf& getBasisCoefficients();
    const Eigen::VectorXf& getVelocityField();
    const Eigen::VectorXf& getVorticityField();
    float getMaxVorticity();
    float getMinVorticity();
protected:
    virtual void buildLaplace()=0;
    virtual void buildAdvection()=0;

    void buildEigenFunctions();

    float minRotation;
    float maxRotation;

    unsigned int nEigenFunctions;
    unsigned int resolution;
    float timeStep;
    float viscosity;

    Mesh2D* mesh;
    DECMesh2D decMesh;

    Eigen::SparseMatrix<float> curl;

    std::vector<Eigen::VectorXf> eigenFunctions;

    Eigen::VectorXf vorticityField;
    Eigen::VectorXf velocityField;

    std::vector<Eigen::MatrixXf> advection;

    Eigen::VectorXf eigenValues;
    Eigen::VectorXf basisCoeff;
    Eigen::MatrixXf velBasisField;
    Eigen::MatrixXf vortBasisField;
private:
};

#endif // ABSTRACTSOLVER_H
