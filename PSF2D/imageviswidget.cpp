#include "imageviswidget.h"
#include <QPainter>

ImageVisWidget::ImageVisWidget(QWidget *parent) : QWidget(parent)
{

}

void ImageVisWidget::setImage(QImage image)
{
    originalImage = image;
}

void ImageVisWidget::paintEvent(QPaintEvent* event)
{
    QPainter painter(this);
    painter.drawImage(QPoint(0,0),originalImage);
    painter.end();
}
