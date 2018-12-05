#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <Eigen/Eigen>
#include "imageviswidget.h"
#include "voxelviswidget.h"
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
    void intigrate();

    void changeToImageView();
    void changeToVoxelView();
protected:
private:
    void voxelize(unsigned int w,unsigned int h);
    bool checkVoxel(unsigned char* offs,int w,int h);

    Mesh2D* mesh;

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
