#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <Eigen/Eigen>
#include "imageviswidget.h"
#include "voxelviswidget.h"
#include "abstractsolver.h"
#include "spectralfluidssolver2d.h"
#include "spectralfluidssolver2domp.h"
#include "mesh2d.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    void openFile();
    void changeViscosity(double visc);
    void changeNumberEigenfunctions(int n);
    void changeResolution(int n);
    void changeTimestep(double step);
    void intigrate();

    void setVelocityNormalization(int);
    void showGrid(bool showGrid);
    void showVelocity(bool showVelocity);
    void showVorticity(bool showVorticity);
    void showMesh(bool showMesh);

    void changeToImageView();
    void changeToVoxelView();
private:

    AbstractSolver* solver;
    SpectralFluidsSolver2D* spectralFluidsSolver;
    SpectralFluidsSolver2DOMP* spectralFluidsSolverOMP;

    QTimer updateTimer;

    ImageVisWidget* imageVis;
    VoxelVisWidget* voxelVis;

    QImage originalImage;
    QImage voxelImage;
    QImage flowImage;
    QImage particleImage;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
