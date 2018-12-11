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

    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(update()));
    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(intigrate()));

    mesh = NULL;
    updateTimer.setInterval(1000/30);
}

MainWindow::~MainWindow()
{
    delete ui;
    if(mesh!=NULL)
    {
        delete mesh;
    }
}

void MainWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,"Mesh File");
    if(!fileName.isEmpty())
    {
        updateTimer.stop();
        //Convert image to Greyscale
        originalImage = QImage(fileName);
        originalImage = originalImage.convertToFormat(QImage::Format_Grayscale8);

        //Init Mesh
        if(mesh!=NULL)
        {
            delete mesh;
        }
        mesh = new Mesh2D(originalImage.width(),originalImage.height(),32,32,originalImage.bits());


        ui->tabImageView->setImage(originalImage);
        ui->tabVoxelView->setMesh(mesh);
        ui->tabFlowView->setMesh(mesh);
        //ui->renderWidget->setMinimumSize(originalImage.width(),originalImage.height());
        //ui->renderWidget->setPixmap(QPixmap::fromImage(originalImage));

        updateTimer.start();
    }
}

void MainWindow::intigrate()
{
    mesh->integrate();
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
