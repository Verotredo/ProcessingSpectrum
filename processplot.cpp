#include "processplot.h"

ProcessPlot::ProcessPlot()
{
    mv= new MatVector();
}
double ProcessPlot::errorS(QVector<double> &i1,QVector<double> &i2){

    QVector<double> i=mv->minus(i1,i2);
    return mv->norm(i);

}

double ProcessPlot::lat(QVector<QVector<double> > &input){
    int n=input[0].length();
    double l=0, y_max=0;
    for(int i=0;i<n;i++){
        if(input[1][i]>y_max){
            y_max=input[1][i];
            l=input[0][i];
        }
    }
    return l;
}
QVector<QVector<double> > ProcessPlot::maxPoints(QVector<QVector<double> > &input,int begin,int mid, int end){
    QVector<double> x ;
    x<<0<<0<<0<<0;
    QVector<double> y,otn;
    y<<0<<0<<0<<0;
    int n=input[0].length();
    for(int i=1; i<n-1;i++){
        if(input[0][i]<begin){
            if(input[1][i]>y[0]){
                x[0]=input[0][i];
                y[0]=input[1][i];
            }
        }
        else if(input[0][i]<mid){
            if(input[1][i]>y[1]){
                x[1]=input[0][i];
                y[1]=input[1][i];
            }
        }
        else if(input[0][i]<end){
            if(input[1][i]>y[2]){
                x[2]=input[0][i];
                y[2]=input[1][i];
            }
        }
        else
            if(input[1][i]>y[3]){
                x[3]=input[0][i];
                y[3]=input[1][i];
            }

    }
    QVector<QVector<double> > v;
    otn<<y[0]/y[2]<<y[1]/y[2]<<y[3]/y[2];
    v<<x<<y<<otn;
    return v;
}

QVector<QVector<double> > ProcessPlot::normalize(QVector<QVector<double> > &input){
    QVector<QVector<double> >v;
    QVector<double> m=mv->minus(input[1],mv->min(input[1]));
    //qDebug() <<mv->sum(m);
    v<<input[0]<<mv->divide(m,mv->sum(m));
    return v;
}

QVector<double>   ProcessPlot::statistic(QVector<QVector<double> > &input){
    double M=mv->aTb(input[0],input[1]);
    double M2=mv->aTb(mv->powV(mv->minus(input[0],M),2),input[1]);
    double M3=mv->aTb(mv->powV(mv->minus(input[0],M),3),input[1]);
    double M4=mv->aTb(mv->powV(mv->minus(input[0],M),4),input[1]);
    double S=sqrt(M2);
    double A=M3/pow(S,3);
    double E=M4/pow(S,4)-3;
    QVector<double>res;
    res<<M<<A<<E;
    return res;
}

QVector<QVector<double> >   ProcessPlot::interplt(QVector<QVector<double> > &input, int diff){
    double min = mv->min(input[0]);
    double max = mv->max(input[0]);
    QVector<double> x;
    for(int i=min;i<max;i=i+diff){
        x<<i;
        if(i+diff>=max) {
            x<<max;
            break;
        }
    }
    int n=input[0].length();
    QVector<double> y;


    int nx=x.length();
        spline1dinterpolant sp;
        real_1d_array xold,yold,xnew,ynew;
        double *dxold =&input[0][0], *dyold=&input[1][0], *dxnew =&x[0],*dynew;
        xold.setcontent(n,dxold); yold.setcontent(n,dyold); xnew.setcontent(nx,dxnew);
        spline1dconvcubic(xold, yold, xnew, ynew);
        dynew= &ynew[0];
        for(int i=0;i<x.length();i++){
            y<<dynew[i];
        }
    QVector<QVector<double> > res;
    res<<x<<y;
    return res;
}



const QString ProcessPlot::allFileToString(QFile &aFile)
{
    if (!aFile.open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "Error opening file!" << endl;
        return NULL;
    }
    QTextStream in(&aFile);
    return in.readAll();
}

void ProcessPlot::setDataToVector(const QStringList &aStringList,
                     QVector< QVector <double> > &aVector)
{
    size_t x = aStringList.length(); // Count of line
    size_t y = aStringList.at(0).count(" ") + 1; // Count of digits in line
    //qDebug() << x<< y;
    for (size_t i = 0; i < y; ++i) {
        QVector<double> temp_vector;
        for (size_t j = 0; j < x; ++j) {
            QString s=aStringList.at(j).split(" ").at(i);
            //qDebug()<<s;
            temp_vector.push_back(s.toDouble());
        }
        aVector.push_back(temp_vector);
    }
}

void ProcessPlot::printVector(const QVector< QVector <double> > &aVector)
{
    for (int i = 0; i < aVector.size(); ++i) {
        for (int j = 0; j < aVector.at(0).size(); ++j) {
             qDebug() << aVector.at(i).at(j) << "\t";
        }
        qDebug() << endl;
    }
}
QVector<QVector<double> > ProcessPlot::input(const QString &name, int a){
    QVector< QVector <double> > vector,v;
    QFile file(name);
    QString s= allFileToString(file);
    s.replace(","," ");
    s.replace("Channel 1","");
    s.replace("Axis [nm],ROI1 [],ROI2 [],ROI3 [],ROI4 []","");
    s.replace("Axis [nm],ROI1 [],ROI2 [],ROI3 []","");
    s.replace("Axis [nm],ROI1 [],ROI2 []","");
    s.replace("Axis [nm],ROI1 []","");
    QStringList sl = s.split("\n");
    sl.removeAll("");

    setDataToVector(sl, vector);

    //printVector(vector);
    if(a<vector.length()) v<<vector[0]<<vector[a];
    else v=vector;
    return v;

}


