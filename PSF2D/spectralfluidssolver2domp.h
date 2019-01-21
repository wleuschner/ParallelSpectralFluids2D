#ifndef SPECTRALFLUIDSSOLVER2DOMP_H
#define SPECTRALFLUIDSSOLVER2DOMP_H
#include "abstractsolver.h"
#include "decmesh2d.h"
#include "mesh2d.h"

class SpectralFluidsSolver2DOMP : public AbstractSolver
{
public:
    SpectralFluidsSolver2DOMP();
    void integrate();
protected:
    void buildLaplace();
    void buildAdvection();
};

#endif // SPECTRALFLUIDSSOLVER2D_H
