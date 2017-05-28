#ifndef PROCESSPLOT_H
#define PROCESSPLOT_H
#include <QDebug>
#include <QObject>
#include <QVarLengthArray>
#include <matvector.h>
#include <QFile>
#include "cpp/src/interpolation.h"
using namespace alglib;
class ProcessPlot
{
    MatVector *mv;
public:
    ProcessPlot();
    double errorS(QVector<double> &i1,QVector<double> &i2);
    QVector<QVector<double> > interplt(QVector<QVector<double> > &input);
    double  lat(QVector<QVector<double> > &input);
    QVector<QVector<double> > maxPoints(QVector<QVector<double> > &input);
    QVector<QVector<double> > normalize(QVector<QVector<double> > &input);
    QVector<double>   statistic(QVector<QVector<double> > &input);
    double ** peaks(double ** m);
    double ** differencial(double ** m);
    QVector<QVector<double> > input(const QString &name, int a);
    void printVector(const QVector< QVector <double> > &aVector);
    void setDataToVector(const QStringList &aStringList,
                         QVector< QVector <double> > &aVector);
    const QString allFileToString(QFile &aFile);
};

#endif // PROCESSPLOT_H
