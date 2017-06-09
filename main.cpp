#include <QApplication>
#include <processplot.h>
#include <QDebug>
#include <QVarLengthArray>
#include <plotwidget.h>
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    PlotWidget w;
    w.resize(900, 680);
    w.setWindowTitle("Extractor");
    w.show();
    return a.exec();
}
