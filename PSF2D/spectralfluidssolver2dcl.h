#ifndef SPECTRALFLUIDSSOLVER2DCL_H
#define SPECTRALFLUIDSSOLVER2DCL_H
#include "abstractsolver.h"
#include "decmesh2d.h"
#include "mesh2d.h"
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>

class SpectralFluidsSolver2DCL : public AbstractSolver
{
public:
    SpectralFluidsSolver2DCL();
    void integrate();
protected:
    void buildLaplace();
    void buildAdvection();
private:
    viennacl::compressed_matrix<double> vclCurl;

    std::vector<viennacl::vector<double>> vclEigenFunctions;

    viennacl::vector<double> vclVorticityField;
    viennacl::vector<double> vclVelocityField;

    std::vector<viennacl::matrix<double>> vclAdvection;

    viennacl::vector<double> vclEigenValues;
    viennacl::vector<double> vclBasisCoeff;
    viennacl::matrix<double> vclVelBasisField;
    viennacl::matrix<double> vclVortBasisField;
};
#endif // SPECTRALFLUIDSSOLVER2DCL_H
