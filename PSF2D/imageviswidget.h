#ifndef IMAGEVISWIDGET_H
#define IMAGEVISWIDGET_H

#include <QObject>
#include <QWidget>

class ImageVisWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ImageVisWidget(QWidget *parent = nullptr);
    void setImage(QImage image);
protected:
    void paintEvent(QPaintEvent* event);
signals:

public slots:
private:
    QImage originalImage;
};

#endif // IMAGEVISWIDGET_H
