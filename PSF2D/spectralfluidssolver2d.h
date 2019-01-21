#ifndef SPECTRALFLUIDSSOLVER2D_H
#define SPECTRALFLUIDSSOLVER2D_H
#include "abstractsolver.h"
#include "decmesh2d.h"
#include "mesh2d.h"

class SpectralFluidsSolver2D : public AbstractSolver
{
public:
    SpectralFluidsSolver2D();
    void integrate();
protected:
    void buildLaplace();
    void buildAdvection();
};

#endif // SPECTRALFLUIDSSOLVER2D_H
