#ifndef FLOWVISWIDGET_H
#define FLOWVISWIDGET_H

#include <QObject>
#include <QImage>
#include <QWidget>
#include <QOpenGLWidget>
#include "abstractsolver.h"

class FlowVisWidget : public QOpenGLWidget
{
    Q_OBJECT
public:
    enum VelocityNormalization
    {
        NO_NORMALIZATION,
        UNIT_NORMALIZATION,
        TIMESTEP_NORMALIZATION
    };

    explicit FlowVisWidget(QWidget *parent = nullptr);
    void setVelocityNormalization(VelocityNormalization state);
    void setVelocityVisible(bool state);
    void setVorticityVisible(bool state);
    void setGridVisisble(bool state);
    void setImage(QImage image);
    void setSolver(AbstractSolver* solver);
    void resizeGrid();
protected:
    void initializeGL();
    //void paintGL();
    void resizeGL(int w, int h);
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
signals:

public slots:
private:
    QImage renderBuffer;
    QImage rotationField;
    unsigned int* rotationFieldPixels;
    QImage originalImage;

    std::vector<glm::dvec2> velocities;
    std::vector<glm::dvec2> faceCenters;

    std::vector<glm::dvec2> gridEdges;

    double meshAlpha;
    bool vorticityVisible;
    bool velocityVisible;
    bool gridVisible;
    unsigned int velocityNormalizationState;

    bool recording;
    unsigned int imageNo;
    AbstractSolver* solver;
};

#endif // FLOWVISWIDGET_H
