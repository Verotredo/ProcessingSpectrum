#include <QApplication>
#include <processplot.h>
#include <QDebug>
#include <QVarLengthArray>
#include <plotwidget.h>
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QVector<QVector<double>> v405, v458, v476, v488,v496,v514,v543,v633;
    PlotWidget w;
    w.resize(900, 680);
    w.show();
    return a.exec();
}
