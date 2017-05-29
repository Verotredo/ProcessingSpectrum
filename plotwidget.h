#ifndef PLOTWIDGET_H
#define PLOTWIDGET_H
#include <QDebug>
#include <QObject>
#include <QWidget>
#include <QtCharts>
#include <processplot.h>
#include <QFileDialog>
#include <matvector.h>
#include <QFile>
using namespace QtCharts;
class PlotWidget : public QWidget
{
    Q_OBJECT
public:
    QVector<QVector<double>> v405, v458, v476, v488,v496,v514,v543,v633;
    explicit PlotWidget(QWidget *parent = 0);
    QVector<QVector<double>>   process(const QString &file, int a,const QString &name);
public slots:
    void addSeries(QVector<QVector<double>> &input,const QString &number, int i);
    void removeSeries(int i);
    void connectMarkers( int i);
    void disconnectMarkers( int i);

    void handleMarkerClicked();
    void f405Chosen();
    void f458Chosen();
    void f476Chosen();
    void f488Chosen();
    void f496Chosen();
    void f514Chosen();
    void f543Chosen();
    void f633Chosen();
    void loaded();
    void wroteResult();

private:

    QVector<QChart *>m_chart;
    int amount;
    QList<QLineSeries *> m_series;
    ProcessPlot *pp;
    MatVector *mv = new MatVector();
    QVector<double> res;
    QVector<QChartView *>m_chartView;
    QGridLayout *m_mainLayout;
    QGridLayout *m_fontLayout;
    QSpinBox *sb, *sbbeg,*sbmid,*sbend;
    QString f405="",f458="",f476="",f488="",f496="",f514="",f543="",f633="";

};

#endif // PLOTWIDGET_H
