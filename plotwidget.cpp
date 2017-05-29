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


    QPushButton *b405 = new QPushButton("405",this);
    connect( b405, SIGNAL(clicked()),this,SLOT(f405Chosen()));
    m_mainLayout->addWidget(b405,0,12);
    QPushButton *b458 = new QPushButton("458",this);
    connect( b458, SIGNAL(clicked()),this,SLOT(f458Chosen()));
    m_mainLayout->addWidget(b458,1,12);
    QPushButton *b476 = new QPushButton("476",this);
    connect( b476, SIGNAL(clicked()),this,SLOT(f476Chosen()));
    m_mainLayout->addWidget(b476,2,12);
    QPushButton *b488 = new QPushButton("488",this);
    connect( b488, SIGNAL(clicked()),this,SLOT(f488Chosen()));
    m_mainLayout->addWidget(b488,3,12);
    QPushButton *b496 = new QPushButton("496",this);
    connect( b496, SIGNAL(clicked()),this,SLOT(f496Chosen()));
    m_mainLayout->addWidget(b496,4,12);
    QPushButton *b514 = new QPushButton("514",this);
    connect( b514, SIGNAL(clicked()),this,SLOT(f514Chosen()));
    m_mainLayout->addWidget(b514,5,12);
    QPushButton *b543 = new QPushButton("543",this);
    connect( b543, SIGNAL(clicked()),this,SLOT(f543Chosen()));
    m_mainLayout->addWidget(b543,6,12);
    QPushButton *b633 = new QPushButton("633",this);
    connect( b633, SIGNAL(clicked()),this,SLOT(f633Chosen()));
    m_mainLayout->addWidget(b633,7,12);
    QLabel *label1 = new QLabel("Номер столбца");
    m_mainLayout->addWidget(label1,8,12);
    sb = new QSpinBox(this);
    m_mainLayout->addWidget(sb,9,12);

    QLabel *label2 = new QLabel("Границы");
    m_mainLayout->addWidget(label2,10,12);
    sbbeg = new QSpinBox(this);
    sbbeg->setRange(-100,100);
    sbbeg->setValue(-70);
    m_mainLayout->addWidget(sbbeg,11,12);
    sbmid = new QSpinBox(this);
    sbmid->setRange(-100,100);
    sbmid->setValue(-20);
    m_mainLayout->addWidget(sbmid,12,12);
    sbend = new QSpinBox(this);
    sbend->setRange(-100,100);
    sbend->setValue(30);
    m_mainLayout->addWidget(sbend,13,12);


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
    if(v.length()>2) {
        (new QErrorMessage(this))->showMessage("Неверный номер");
        return v;
    }
    amount++;
    addSeries(v,name,0);
    connectMarkers(0);
    v=pp->interplt(v);
    addSeries(v,name,1);
    connectMarkers(1);
    v=pp->normalize(v);
    addSeries(v,name,2);
    connectMarkers(2);
    return v;
}
void  PlotWidget::loaded()
{
    if(amount==9) amount--;
    for(int j=0;j<3;j++){
        for(int i=0; i<amount; i++){
            removeSeries(j);
        }
    }
    amount=0;
    double l=0;
    QVector<QVector<double>> p405,p458,p476,p488,p496,p514,p543,p633;
    double e_405,e_458,e_476,e_496,e_514,e_543,e_633;
    QVector<double> st405,st458,st476,st488,st496,st514,st543,st633;
    if(sb->value()>0){
        if(f405!="") {
            v405=process(f405,sb->value(),"405");
        }
        if(f458!="") {
            v458=process(f458,sb->value(),"458");
        }
        if(f476!="") {
            v476=process(f476,sb->value(),"476");
        }
        if(f488!="") {
            v488=process(f488,sb->value(),"488");
        }
        if(f496!="") {
            v496=process(f496,sb->value(),"496");
        }
        if(f514!="") {
            v514=process(f514,sb->value(),"514");
        }
        if(f543!="") {
            v543=process(f543,sb->value(),"543");
        }
        if(f633!="") {
            v633=process(f633,sb->value(),"633");
        }
        if(f458!="") {
            l=pp->lat(v458);
            amount++;
            if(f405!=""){
                v405[0]=mv->minus(v405[0],l);
                addSeries(v405,"405",3);
                connectMarkers(3);
                p405=pp->maxPoints(v405,sbbeg->value(),sbmid->value(), sbend->value());
                e_405 = pp->errorS(v488[1],v405[1]);
                st405 = pp->statistic(v405);
            }
            v458[0]=mv->minus(v458[0],l);
            addSeries(v458,"458",3);
            connectMarkers(3);
            p458=pp->maxPoints(v458,sbbeg->value(),sbmid->value(), sbend->value());
            e_458 = pp->errorS(v488[1],v458[1]);
            st458 = pp->statistic(v458);
            if(f476!="") {
                v476[0]=mv->minus(v476[0],l);
                addSeries(v476,"476",3);
                connectMarkers(3);
                p476=pp->maxPoints(v476,sbbeg->value(),sbmid->value(), sbend->value());
                e_476 = pp->errorS(v488[1],v476[1]);
                st476 = pp->statistic(v476);
            }
            if(f488!="") {
                v488[0]=mv->minus(v488[0],l);
                addSeries(v488,"488",3);
                connectMarkers(3);
                p488=pp->maxPoints(v488,sbbeg->value(),sbmid->value(), sbend->value());
                st488 = pp->statistic(v488);
            }
            if(f496!="") {
                v496[0]=mv->minus(v496[0],l);
                addSeries(v496,"496",3);
                connectMarkers(3);
                p496=pp->maxPoints(v496,sbbeg->value(),sbmid->value(), sbend->value());
                e_496 = pp->errorS(v488[1],v496[1]);
                st496 = pp->statistic(v496);
            }
            if(f514!="") {
                v514[0]=mv->minus(v514[0],l);
                addSeries(v514,"514",3);
                connectMarkers(3);
                p514=pp->maxPoints(v514,sbbeg->value(),sbmid->value(), sbend->value());
                e_514 = pp->errorS(v488[1],v514[1]);
                st514 = pp->statistic(v514);
            }
            if(f543!="") {
                v543[0]=mv->minus(v543[0],l);
                addSeries(v543,"543",3);
                connectMarkers(3);
                p543=pp->maxPoints(v543,sbbeg->value(),sbmid->value(), sbend->value());
                e_543 = pp->errorS(v488[1],v543[1]);
                st543 = pp->statistic(v543);
            }
            if(f633!="") {
                v633[0]=mv->minus(v633[0],l);
                addSeries(v633,"633",3);
                connectMarkers(3);
                p633=pp->maxPoints(v633,sbbeg->value(),sbmid->value(), sbend->value());
                e_633 = pp->errorS(v488[1],v633[1]);
                st633 = pp->statistic(v633);
            }
            res=p405[2]+p458[2]+p476[2]+p488[2]+p496[2]+p514[2]+p543[2]+p633[2];
            res<<e_405<<e_458<<e_476<<e_496<<e_514<<e_543<<e_633;
            res<<st405[0]<<st458[0]<<st476[0]<<st488[0]<<st496[0]<<st514[0]<<st543[0]<<st633[0];
            res<<st405[1]<<st458[1]<<st476[1]<<st488[1]<<st496[1]<<st514[1]<<st543[1]<<st633[1];
            res<<st405[2]<<st458[2]<<st476[2]<<st488[2]<<st496[2]<<st514[2]<<st543[2]<<st633[2];
        }


    }
}

void PlotWidget::wroteResult()
{
    QString name = QFileDialog::getOpenFileName(0,"Open ","","*.txt");
    QFile result(name);
    QTextStream stream(&result);
    QString f=f488;
    f.remove("_488.txt");
    f.append("_"+QString::number(sb->value()));
    stream<<f;
    if(result.open(QIODevice::WriteOnly|  QIODevice::Append)){
        for(int i=0; i<res.length();i++)
        stream<<','<<res[i];
    }
    stream<<'\n';
   result.close();
}

void PlotWidget::f405Chosen()
{
    f405 = QFileDialog::getOpenFileName(0,"Open 405","","*.txt");
}
void PlotWidget::f458Chosen()
{
    f458 = QFileDialog::getOpenFileName(0,"Open 458","","*.txt");
}
void PlotWidget::f476Chosen()
{
    f476 = QFileDialog::getOpenFileName(0,"Open 476","","*.txt");
}
void PlotWidget::f488Chosen()
{
    f488 = QFileDialog::getOpenFileName(0,"Open 488","","*.txt");
}
void PlotWidget::f496Chosen()
{
    f496 = QFileDialog::getOpenFileName(0,"Open 496","","*.txt");
}
void PlotWidget::f514Chosen()
{
    f514 = QFileDialog::getOpenFileName(0,"Open 514","","*.txt");
}
void PlotWidget::f543Chosen()
{
    f543 = QFileDialog::getOpenFileName(0,"Open 543","","*.txt");
}
void PlotWidget::f633Chosen()
{
    f633 = QFileDialog::getOpenFileName(0,"Open 633","","*.txt");
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

    if (amount == 1 || amount==9) {
        m_chart[i]->createDefaultAxes();
    }
}

void PlotWidget::removeSeries(int i)
{
    if (m_series.count() > 0) {
        QLineSeries *series = m_series.last();
        m_chart[i]->removeSeries(series);
        m_series.removeLast();
        delete series;
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
