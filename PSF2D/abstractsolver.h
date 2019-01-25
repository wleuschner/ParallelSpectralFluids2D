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

    void setInitialVelocityField(const Eigen::VectorXd& field);
    void setInitialVorticityField(const Eigen::VectorXd& field);

    void setNumberEigenFunctions(unsigned int n);
    void setResolution(unsigned int res);
    void setViscosity(double visc);
    void setTimestep(double timestep);

    Mesh2D* getMesh();
    DECMesh2D& getDECMesh();

    unsigned int getNumEigenFunctions();
    double getTimestep();
    const Eigen::VectorXd& getEigenFunction(unsigned int n);
    const Eigen::MatrixXd& getVelocityBasisField();
    const Eigen::MatrixXd& getVorticityBasisField();
    const Eigen::VectorXd& getBasisCoefficients();
    const Eigen::VectorXd& getVelocityField();
    const Eigen::VectorXd& getVorticityField();
    double getMaxVorticity();
    double getMinVorticity();
protected:
    virtual void buildLaplace()=0;
    virtual void buildAdvection()=0;

    void buildEigenFunctions();

    double minRotation;
    double maxRotation;

    unsigned int nEigenFunctions;
    unsigned int resolution;
    double timeStep;
    double viscosity;

    Mesh2D* mesh;
    DECMesh2D decMesh;

    Eigen::SparseMatrix<double> curl;

    std::vector<Eigen::VectorXd> eigenFunctions;

    Eigen::VectorXd vorticityField;
    Eigen::VectorXd velocityField;

    std::vector<Eigen::MatrixXd> advection;

    Eigen::VectorXd eigenValues;
    Eigen::VectorXd basisCoeff;
    Eigen::MatrixXd velBasisField;
    Eigen::MatrixXd vortBasisField;
private:
};

#endif // ABSTRACTSOLVER_H
