#include "plotwidget.h"

PlotWidget::PlotWidget(QWidget *parent) :
    QWidget(parent)
{
    pp = new ProcessPlot();
    m_mainLayout = new QGridLayout();
    for(int i=0;i<4;i++){
        m_chart<< new QChart();
        m_chartView<< new QChartView(m_chart[i], this);
        m_chart[i]->legend()->setVisible(true);
        m_chart[i]->legend()->setAlignment(Qt::AlignBottom);
        m_chartView[i]->setRubberBand(QChartView::RectangleRubberBand);
        m_chartView[i]->setRenderHint(QPainter::Antialiasing);
    }
    m_chart[0]->setTitle("Исходные графики");
    m_chart[1]->setTitle("Интерполированные графики");
    m_chart[2]->setTitle("Нормированные графики");
    m_chart[3]->setTitle("Центрированные графики");
    m_mainLayout->addWidget(m_chartView[0], 0, 0, 10, 6);
    m_mainLayout->addWidget(m_chartView[1], 0, 6, 10, 6);
    m_mainLayout->addWidget(m_chartView[2], 10, 0, 7, 6);
    m_mainLayout->addWidget(m_chartView[3], 10, 6, 7, 6);
    setLayout(m_mainLayout);

    b=-70; m=-20; e=30;
    QLabel *label1 = new QLabel("Загрузка данных: ");
    m_mainLayout->addWidget(label1,0,12);

    b405 = new QPushButton("405",this);
    connect( b405, SIGNAL(clicked()),this,SLOT(f405Chosen()));
    m_mainLayout->addWidget(b405,1,12);
    b458 = new QPushButton("458",this);
    connect( b458, SIGNAL(clicked()),this,SLOT(f458Chosen()));
    m_mainLayout->addWidget(b458,2,12);
    b476 = new QPushButton("476",this);
    connect( b476, SIGNAL(clicked()),this,SLOT(f476Chosen()));
    m_mainLayout->addWidget(b476,3,12);
    b488 = new QPushButton("488",this);
    connect( b488, SIGNAL(clicked()),this,SLOT(f488Chosen()));
    m_mainLayout->addWidget(b488,4,12);
    b496 = new QPushButton("496",this);
    connect( b496, SIGNAL(clicked()),this,SLOT(f496Chosen()));
    m_mainLayout->addWidget(b496,5,12);
    b514 = new QPushButton("514",this);
    connect( b514, SIGNAL(clicked()),this,SLOT(f514Chosen()));
    m_mainLayout->addWidget(b514,6,12);
    b543 = new QPushButton("543",this);
    connect( b543, SIGNAL(clicked()),this,SLOT(f543Chosen()));
    m_mainLayout->addWidget(b543,7,12);
    b633 = new QPushButton("633",this);
    connect( b633, SIGNAL(clicked()),this,SLOT(f633Chosen()));
    m_mainLayout->addWidget(b633,8,12);
    QLabel *label2 = new QLabel("Номер столбца - бактерии (>0)");
    m_mainLayout->addWidget(label2,9,12);
    sb1 = new QSpinBox(this);
    m_mainLayout->addWidget(sb1,10,12);

    QLabel *label3 = new QLabel("Размер деления при интерп.");
    m_mainLayout->addWidget(label3,11,12);
    sb2 = new QSpinBox(this);
    sb2->setMinimum(1);
    sb2->setMaximum(5);
    sb2->setValue(1);
    m_mainLayout->addWidget(sb2,12,12);

    QPushButton *plot = new QPushButton("Отобразить графики",this);
    connect( plot, SIGNAL(clicked()),this,SLOT(loaded()));
    m_mainLayout->addWidget(plot,15,12);
    QPushButton *write = new QPushButton("Загрузить в файл",this);
    connect( write, SIGNAL(clicked()),this,SLOT(wroteResult()));
    m_mainLayout->addWidget(write,16,12);
}
QVector<QVector<double>>   PlotWidget::process(const QString &file, int a,const QString &name){
    QVector<QVector<double>> v;
    v=pp->input(file,a);
    if(v.length()>2 &&  amount==0) {
        (new QErrorMessage(this))->showMessage("Неверный номер");
        return v;
    }
    amount++;
    addSeries(v,name,0);
    connectMarkers(0);
    v=pp->interplt(v,sb2->value());
    addSeries(v,name,1);
    connectMarkers(1);
    v=pp->normalize(v);
    addSeries(v,name,2);
    connectMarkers(2);
    return v;
}
void  PlotWidget::loaded()
{
    if(amount>0){
        removeSeries();
    }
    amount=0;
    double l=0;
    QVector<QVector<double>> p405,p458,p476,p488,p496,p514,p543,p633;
    double e_405,e_458,e_476,e_496,e_514,e_543,e_633;
    QVector<double> st405,st458,st476,st488,st496,st514,st543,st633;
    if(sb1->value()>0){
        if(f405!="") {
            v405=process(f405,sb1->value(),"405");
        }
        if(f458!="") {
            v458=process(f458,sb1->value(),"458");
        }
        if(f476!="") {
            v476=process(f476,sb1->value(),"476");
        }
        if(f488!="") {
            v488=process(f488,sb1->value(),"488");
        }
        if(f496!="") {
            v496=process(f496,sb1->value(),"496");
        }
        if(f514!="") {
            v514=process(f514,sb1->value(),"514");
        }
        if(f543!="") {
            v543=process(f543,sb1->value(),"543");
        }
        if(f633!="") {
            v633=process(f633,sb1->value(),"633");
        }
        if(f458!="" && f488!="") {
            l=pp->lat(v458);
            if(f405!=""){
                v405[0]=mv->minus(v405[0],l);
                //qDebug()<<v405[0];
                addSeries(v405,"405",3);
                connectMarkers(3);
                p405=pp->maxPoints(v405,b,m, e);
                e_405 = pp->errorS(v488[1],v405[1]);
                st405 = pp->statistic(v405);

            }
            v458[0]=mv->minus(v458[0],l);
            addSeries(v458,"458",3);
            connectMarkers(3);
            p458=pp->maxPoints(v458,b,m, e);
            e_458 = pp->errorS(v488[1],v458[1]);
            st458 = pp->statistic(v458);
            if(f476!="") {
                v476[0]=mv->minus(v476[0],l);
                addSeries(v476,"476",3);
                connectMarkers(3);
                p476=pp->maxPoints(v476,b,m, e);
                e_476 = pp->errorS(v488[1],v476[1]);
                st476 = pp->statistic(v476);
            }
            v488[0]=mv->minus(v488[0],l);
            addSeries(v488,"488",3);
            connectMarkers(3);
            p488=pp->maxPoints(v488,b,m, e);
            st488 = pp->statistic(v488);
            if(f496!="") {
                v496[0]=mv->minus(v496[0],l);
                addSeries(v496,"496",3);
                connectMarkers(3);
                p496=pp->maxPoints(v496,b,m, e);
                e_496 = pp->errorS(v488[1],v496[1]);
                st496 = pp->statistic(v496);
            }
            if(f514!="") {
                v514[0]=mv->minus(v514[0],l);
                addSeries(v514,"514",3);
                connectMarkers(3);
                p514=pp->maxPoints(v514,b,m, e);
                e_514 = pp->errorS(v488[1],v514[1]);
                st514 = pp->statistic(v514);
            }
            if(f543!="") {
                v543[0]=mv->minus(v543[0],l);
                addSeries(v543,"543",3);
                connectMarkers(3);
                p543=pp->maxPoints(v543,b,m, e);
                e_543 = pp->errorS(v488[1],v543[1]);
                st543 = pp->statistic(v543);
            }
            if(f633!="") {
                v633[0]=mv->minus(v633[0],l);
                addSeries(v633,"633",3);
                connectMarkers(3);
                p633=pp->maxPoints(v633,b,m, e);
                e_633 = pp->errorS(v488[1],v633[1]);
                st633 = pp->statistic(v633);
            }
            if(f405!="" && f476!="" && f496!="" && f514!="" && f543!="" && f633!=""){
                res=p405[2]+p458[2]+p476[2]+p488[2]+p496[2]+p514[2]+p543[2]+p633[2];
                res<<e_405<<e_458<<e_476<<e_496<<e_514<<e_543<<e_633;
                res<<st405[0]<<st458[0]<<st476[0]<<st488[0]<<st496[0]<<st514[0]<<st543[0]<<st633[0];
                res<<st405[1]<<st458[1]<<st476[1]<<st488[1]<<st496[1]<<st514[1]<<st543[1]<<st633[1];
                res<<st405[2]<<st458[2]<<st476[2]<<st488[2]<<st496[2]<<st514[2]<<st543[2]<<st633[2];
            }
        }


    }
}

void PlotWidget::wroteResult()
{   if(res.length()>50){
        QString name = QFileDialog::getOpenFileName(0,"Save in ","","*.txt");

        int ind=f488.lastIndexOf("_488");
        QFileInfo *fi=new QFileInfo(f488);
        QFile result(name);
        QString input=fi->baseName(),number,n_date,cell;
        QTextStream stream(&result);
        if(result.size()==0) stream<<"File,Species,Species+Date, otn_405_1,otn_405_2,otn_405_4, otn_458_1,otn_458_2,otn_458_4, otn_476_1,otn_476_2,otn_476_4, otn_488_1,otn_488_2,otn_488_4,otn_496_1,otn_496_2,otn_496_4, otn_514_1,otn_514_2,otn_514_4, otn_543_1,otn_543_2,otn_543_4, otn_633_1,otn_633_2,otn_633_4, e_405, e_458, e_476, e_496, e_514, e_543, e_633, m_405, m_458, m_476,m_488, m_496, m_514, m_543, m_633, a_405, a_458, a_476,a_488, a_496, a_514, a_543, a_633, E_405, E_458, E_476,E_488, E_496, E_514, E_543, E_633 "<<"\r\n";
        cell=f488.mid(ind-2,2);
        input.remove("_488.txt");
        //input.append("_"+QString::number(sb->value()));
        ind=input.indexOf(" ");
        if(ind >=0) number=input.left(ind);
        else {
            ind=input.indexOf("_");
            number=input.left(ind);
        }
        ind=f488.indexOf(".");
        n_date=f488.mid(ind-3,9);
        stream<<"\r\n"<<number+n_date+cell+"_"+QString::number(sb1->value())<<','<<number<<", "<<number+n_date;
        if(result.open(QIODevice::WriteOnly|  QIODevice::Append)){
            for(int i=0; i<res.length();i++)
                stream<<','<<res[i];
        }
        stream<<"\r\n";
        result.close();
    }
}

void PlotWidget::f405Chosen()
{
    f405 = QFileDialog::getOpenFileName(0,"Open 405","","*.txt");
    if(f405!="")b405->setText("405: +");
}
void PlotWidget::f458Chosen()
{
    f458 = QFileDialog::getOpenFileName(0,"Open 458","","*.txt");
    if(f458!="")b458->setText("458: +");
}
void PlotWidget::f476Chosen()
{
    f476 = QFileDialog::getOpenFileName(0,"Open 476","","*.txt");
    if(f476!="")b476->setText("476: +");
}
void PlotWidget::f488Chosen()
{
    f488 = QFileDialog::getOpenFileName(0,"Open 488","","*.txt");
    if(f488!="")b488->setText("488: +");
}
void PlotWidget::f496Chosen()
{
    f496 = QFileDialog::getOpenFileName(0,"Open 496","","*.txt");
    if(f496!="")b496->setText("496: +");
}
void PlotWidget::f514Chosen()
{
    f514 = QFileDialog::getOpenFileName(0,"Open 514","","*.txt");
    if(f514!="")b514->setText("514: +");
}
void PlotWidget::f543Chosen()
{
    f543 = QFileDialog::getOpenFileName(0,"Open 543","","*.txt");
    if(f543!="")b543->setText("543: +");
}
void PlotWidget::f633Chosen()
{
    f633 = QFileDialog::getOpenFileName(0,"Open 633","","*.txt");
    if(f633!="")b633->setText("633: +");
}

void PlotWidget::addSeries(QVector<QVector<double>> &input, const QString &number, int i)
{
    QLineSeries *series = new QLineSeries();
    m_series.append(series);
    series->setName( number);

    QList<QPointF> data;

    for (int i = 0; i < input[0].length(); i++) {
        data.append(QPointF(input[0][i], input[1][i]));
    }

    series->append(data);
    m_chart[i]->addSeries(series);
    m_chart[i]->createDefaultAxes();

}

void PlotWidget::removeSeries()
{
    for(int k=0;k<4;k++)
    {   //qDebug() <<m_chart[k]->series();
        m_chart[k]->removeAllSeries();
    }
}


void PlotWidget::connectMarkers(int i)
{
    foreach (QLegendMarker* marker, m_chart[i]->legend()->markers()) {
        QObject::disconnect(marker, SIGNAL(clicked()), this, SLOT(handleMarkerClicked()));
        QObject::connect(marker, SIGNAL(clicked()), this, SLOT(handleMarkerClicked()));
    }
}

void PlotWidget::disconnectMarkers(int i)
{
    foreach (QLegendMarker* marker, m_chart[i]->legend()->markers()) {
        QObject::disconnect(marker, SIGNAL(clicked()), this, SLOT(handleMarkerClicked()));
    }
}

void PlotWidget::handleMarkerClicked()
{
    QLegendMarker* marker = qobject_cast<QLegendMarker*> (sender());
    Q_ASSERT(marker);
    switch (marker->type())
    {
        case QLegendMarker::LegendMarkerTypeXY:
        {
        marker->series()->setVisible(!marker->series()->isVisible());
        marker->setVisible(true);
        qreal alpha = 1.0;

        if (!marker->series()->isVisible()) {
            alpha = 0.5;
        }

        QColor color;
        QBrush brush = marker->labelBrush();
        color = brush.color();
        color.setAlphaF(alpha);
        brush.setColor(color);
        marker->setLabelBrush(brush);

        brush = marker->brush();
        color = brush.color();
        color.setAlphaF(alpha);
        brush.setColor(color);
        marker->setBrush(brush);

        QPen pen = marker->pen();
        color = pen.color();
        color.setAlphaF(alpha);
        pen.setColor(color);
        marker->setPen(pen);
        break;
        }
    default:
        {
        qDebug() << "Unknown marker type";
        break;
        }
    }
}
