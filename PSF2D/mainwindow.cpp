#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <iostream>
#include <QPainter>
#include "dec.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    imageVis=NULL;
    voxelVis=NULL;

    connect(ui->actionOpen,SIGNAL(triggered()),this,SLOT(openFile()));
    connect(ui->action_Original,SIGNAL(triggered()),this,SLOT(changeToImageView()));
    connect(ui->action_Voxel,SIGNAL(triggered()),this,SLOT(changeToVoxelView()));

    connect(ui->spinEigenFunctions,SIGNAL(valueChanged(int)),this,SLOT(changeNumberEigenfunctions(int)));
    connect(ui->spinResolution,SIGNAL(valueChanged(int)),this,SLOT(changeResolution(int)));
    connect(ui->spinViscosity,SIGNAL(valueChanged(double)),this,SLOT(changeViscosity(double)));
    connect(ui->spinTimestep,SIGNAL(valueChanged(double)),this,SLOT(changeTimestep(double)));

    connect(ui->comboNormalization,SIGNAL(currentIndexChanged(int)),this,SLOT(setVelocityNormalization(int)));
    connect(ui->chkShowGrid,SIGNAL(toggled(bool)),this,SLOT(showGrid(bool)));
    connect(ui->chkShowVelocity,SIGNAL(toggled(bool)),this,SLOT(showVelocity(bool)));
    connect(ui->chkShowVorticity,SIGNAL(toggled(bool)),this,SLOT(showVorticity(bool)));

    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(intigrate()));
    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(update()));

    spectralFluidsSolver = new SpectralFluidsSolver2D();
    spectralFluidsSolverOMP = new SpectralFluidsSolver2DOMP();
    solver = spectralFluidsSolverOMP;
    updateTimer.setInterval(1000/60);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete spectralFluidsSolver;
}

void MainWindow::openFile()
{
    updateTimer.stop();
    QString fileName = QFileDialog::getOpenFileName(this,"Mesh File");
    if(!fileName.isEmpty())
    {
        //Convert image to Greyscale
        originalImage = QImage(fileName);
        originalImage = originalImage.convertToFormat(QImage::Format_Grayscale8);

        solver->setMesh(new Mesh2D(originalImage.width(),originalImage.height(),32,originalImage.bits()));

        ui->tabImageView->setImage(originalImage);
        ui->tabSpectrumView->setSolver(solver);
        ui->tabSpectrumView->setImage(originalImage);
        ui->tabFlowView->setImage(originalImage);
        ui->tabFlowView->setSolver(solver);
        //ui->renderWidget->setMinimumSize(originalImage.width(),originalImage.height());
        //ui->renderWidget->setPixmap(QPixmap::fromImage(originalImage));

        updateTimer.start();
    }
}

void MainWindow::changeNumberEigenfunctions(int n)
{
    updateTimer.stop();
    solver->setNumberEigenFunctions(n);
    updateTimer.start();
}

void MainWindow::changeResolution(int n)
{
    updateTimer.stop();
    solver->setResolution(n);
    ui->tabFlowView->resizeGrid();
    ui->tabSpectrumView->resizeGrid();
    updateTimer.start();
}

void MainWindow::changeViscosity(double visc)
{
    updateTimer.stop();
    solver->setViscosity(visc);
    updateTimer.start();
}

void MainWindow::changeTimestep(double step)
{
    updateTimer.stop();
    solver->setTimestep(step);
    updateTimer.start();
}

void MainWindow::intigrate()
{
    updateTimer.stop();
    solver->integrate();
    ui->tabFlowView->update();
    updateTimer.start();
}

void MainWindow::setVelocityNormalization(int idx)
{
    switch(idx)
    {
        case 0:
        {
            ui->tabFlowView->setVelocityNormalization(FlowVisWidget::NO_NORMALIZATION);
            break;
        }
        case 1:
        {
            ui->tabFlowView->setVelocityNormalization(FlowVisWidget::UNIT_NORMALIZATION);
            break;
        }
        case 2:
        {
            ui->tabFlowView->setVelocityNormalization(FlowVisWidget::TIMESTEP_NORMALIZATION);
            break;
        }
    }

}

void MainWindow::showGrid(bool showGrid)
{
    ui->tabFlowView->setGridVisisble(showGrid);
}

void MainWindow::showVelocity(bool showVelocity)
{
    ui->tabFlowView->setVelocityVisible(showVelocity);
}

void MainWindow::showVorticity(bool showVorticity)
{
    ui->tabFlowView->setVorticityVisible(showVorticity);
}

void MainWindow::showMesh(bool showMesh)
{
}

void MainWindow::changeToImageView()
{
    if(imageVis!=NULL)
    {
        setCentralWidget(NULL);
        delete imageVis;
    }
    imageVis = new ImageVisWidget(this);
    setCentralWidget(imageVis);
}

void MainWindow::changeToVoxelView()
{
    if(voxelVis!=NULL)
    {
        setCentralWidget(NULL);
        delete voxelVis;
    }
    voxelVis = new VoxelVisWidget(this);
    setCentralWidget(voxelVis);
}
