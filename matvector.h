#ifndef MATVECTOR_H
#define MATVECTOR_H
#include <QVector>
#include <QtMath>
#include <float.h>
#include <QDebug>
class MatVector
{

public:
    MatVector()
    {

    }
    QVector <double> plus(const QVector <double> &i1, const QVector <double> &i2){
        QVector <double> i;
        int i1L=i1.length();
        int i2L=i2.length();
        if(i1L!=i2L) i<<0;
        else {
            for(int j=0;j<i1L;j++){
                i<<i1[j]+i2[j];
            }
        }
        return i;
    }
    QVector <double> minus(const QVector <double> &i1, const QVector <double> &i2){
        QVector <double> i;
        int i1L=i1.length();
        int i2L=i2.length();
        if(i1L!=i2L) i<<0;
        else {
            for(int j=0;j<i1L;j++){
                i<<i1[j]-i2[j];
            }
        }
        return i;
    }
    QVector <double> minus(const QVector <double> &i1, double b){
        QVector <double> i;
        int i1L=i1.length();
        for(int j=0;j<i1L;j++){
            i<<i1[j]-b;
        }
        return i;
    }
    QVector <double> mult(const QVector <double> &i1, double b){
        QVector <double> i;
        int i1L=i1.length();
        for(int j=0;j<i1L;j++){
            i<<i1[j]*b;
        }
        return i;
    }
    double aTb(const QVector <double> &i1, const QVector <double> &i2){
        int i1L=i1.length();
        int i2L=i2.length();
        double d=0;
        if(i1L!=i2L) d=0;
        else {
            for(int j=0;j<i1L;j++){
                d+=i1[j]*i2[j];
            }
        }
        return d;
    }
    QVector <double> divide(const QVector <double> &i1, double b){
        QVector <double> i;
        int i1L=i1.length();
        for(int j=0;j<i1L;j++){
            i<<i1[j]/b;
        }
        return i;
    }
    QVector <double> divide(const QVector <double> &i1, const QVector <double> &i2){
        QVector <double> i;
        int i1L=i1.length();
        for(int j=0;j<i1L;j++){
            i<<i1[j]/i2[j];
        }
        return i;
    }
    QVector <double> powV(const QVector <double> &i1, int b){
        QVector <double> i;
        int i1L=i1.length();
        for(int j=0;j<i1L;j++){
            i<<pow(i1[j],b);
        }
        return i;
    }
    QVector <double> diff(const QVector <double> &i){
        QVector <double> res;
        for(int j=1;j<i.length();j++){
            res<<i[j]-i[j-1];
        }
        return res;
    }
    double norm(const QVector <double> &i){
        double d;
        d=0;
        for(int j=0;j<i.length();j++){
            d=d+i[j]*i[j];
        }
        return qSqrt(d);
    }
    double min(const QVector <double> &i){
        double m=DBL_MAX;
        for(int j=0;j<i.length();j++){
            if(m>i[j]) m=i[j];
        }
        return m;
    }
    double max(const QVector <double> &i){
        double m=DBL_MIN;
        for(int j=0;j<i.length();j++){
            if(m<i[j]) m=i[j];
        }
        return m;
    }
    double sum(const QVector <double> &i){
        double s=0;
        for(int j=0;j<i.length();j++){
            s=s+i[j];
        }
        return s;
    }
};

#endif // MATVECTOR_H
