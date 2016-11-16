#include "rotationaldialog.h"
#include <QFormLayout>

RotationalDialog::RotationalDialog(QWidget* parent):
    QDialog(parent)
{
    openvdbParameter = new QLineEdit(tr(""), this);
    openvdbParameter->setToolTip(tr("Proporción de la diagonal del objeto que será usada como unidad volumétrica.\n Es un valor entre [0..1]. Mientras más pequeño es el valor, la división de volumen es más fina."));
    openvdbParameter->setEnabled(false);

    simplifyParameter = new QLineEdit(tr(""), this);
    simplifyParameter->setToolTip(tr("Número de triángulos  deseable para simplificar. Se recomienda usar este parámetro cuando se usa la división volumétrica."));
    simplifyParameter->setEnabled(false);

    fractionHarris = new QLineEdit(tr("0.01"), this);
    fractionHarris->setToolTip(tr("Proporción de puntos a ser escogidos para detectar simetrías. Es un valor entre [0..1]."));

    thresholdRANSAC = new QLineEdit(tr("0.05"), this);

    alpha = new QLineEdit(tr("0.9"), this);
    alpha->setToolTip("Valor usado para la alineación.");

    iterNonrigid = new QLineEdit(tr("0"), this);
    localParNonrigid = new QLineEdit(tr("2"), this);
    angleNonRigid = new QLineEdit(tr("0.4"), this);
    rotationsNonrigid = new QLineEdit(tr("2"), this);

    checkVolumetric = new QCheckBox(tr("Volumen"), this);
    connect(checkVolumetric, SIGNAL(toggled(bool)), openvdbParameter, SLOT(setEnabled(bool)));

    checkSimplification = new QCheckBox(tr("Simplificar"), this);
    connect(checkSimplification, SIGNAL(toggled(bool)), simplifyParameter, SLOT(setEnabled(bool)));
    checkICP = new QCheckBox(tr("Registro rígido"), this);

    preprocesamientoGB = new QGroupBox(tr("Preprocesamiento"), this);
    preprocesamientoGB->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 3px;}");
    QFormLayout* preproLayout = new QFormLayout;
    preproLayout->addRow(checkVolumetric, openvdbParameter);
    preproLayout->addRow(checkSimplification, simplifyParameter);
    preprocesamientoGB->setLayout(preproLayout);

    seleccionGB = new QGroupBox(tr("Selección de puntos"), this);
    seleccionGB->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 3px;}");
    QFormLayout* seleccionLayout = new QFormLayout;
    seleccionLayout->addRow(tr("Fracción:"), fractionHarris);
    seleccionLayout->addRow(tr("RANSAC"), thresholdRANSAC);
    seleccionGB->setLayout(seleccionLayout);

    alineacionGB = new QGroupBox(tr("Alineación"), this);
    alineacionGB->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 3px;}");
    QFormLayout* alineacionLayout = new QFormLayout;
    alineacionLayout->addRow(tr("Alpha"), alpha);
    alineacionLayout->addRow(tr("Num. iteraciones:"), iterNonrigid);
    alineacionLayout->addRow(tr("Factor localidad:"), localParNonrigid);
    alineacionLayout->addRow(tr("Angulo"), angleNonRigid);
    alineacionLayout->addRow(tr("Num. rotaciones:"), rotationsNonrigid);
    alineacionLayout->addRow(checkICP);
    alineacionGB->setLayout(alineacionLayout);

    buttonCancel = new QPushButton(tr("&Cancel"), this);
    buttonOK = new QPushButton(tr("&Ok"), this);

    connect(buttonCancel, SIGNAL(clicked()), this, SLOT(reject()));
    connect(buttonOK, SIGNAL(clicked()), this, SLOT(confirmDialog()));

    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(buttonCancel);
    buttonLayout->addWidget(buttonOK);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(preprocesamientoGB);
    layout->addWidget(seleccionGB);
    layout->addWidget(alineacionGB);
    layout->addLayout(buttonLayout);

    setLayout(layout);
}

void RotationalDialog::confirmDialog(){
    //Parse dialog and fill PropertySet
    if(checkVolumetric->isChecked()){
        p->addProperty(std::string("performOpenVDB"), std::string("1"));
        p->addProperty(std::string("factor-volumetric"), openvdbParameter->text().toStdString());
    }else{
        p->addProperty(std::string("performOpenVDB"), std::string("0"));
    }

    if(checkSimplification->isChecked()){
        p->addProperty(std::string("performSimplification"), std::string("1"));
        p->addProperty(std::string("num-triangles-decimation"), simplifyParameter->text().toStdString());
    }else{
        p->addProperty(std::string("performSimplification"), std::string("0"));
    }

    //Parámetros Harris
    p->addProperty(std::string("type-neighborhood"), std::string("rings"));
    p->addProperty(std::string("parameter-neighborhood"), std::string("2"));
    p->addProperty(std::string("K"), std::string("0.04"));
    p->addProperty(std::string("ring-maxima-detection"), std::string("1"));
    p->addProperty(std::string("interest-points-selection"), std::string("fraction"));
    p->addProperty(std::string("parameter-selection"), fractionHarris->text().toStdString());
    p->addProperty(std::string("filtering-steps"), std::string("0"));
    p->addProperty(std::string("threshold-ransac"), thresholdRANSAC->text().toStdString());

    //Parámetros alineación
    p->addProperty(std::string("alpha"), alpha->text().toStdString());
    p->addProperty(std::string("iter-non-rigid"), iterNonrigid->text().toStdString());
    p->addProperty(std::string("factor-spacing"), localParNonrigid->text().toStdString());
    p->addProperty(std::string("compatible-angle"), angleNonRigid->text().toStdString());
    p->addProperty(std::string("rotation-factor"), rotationsNonrigid->text().toStdString());

    if(checkICP->isChecked())
        p->addProperty(std::string("enable-icp"), std::string("1"));
    else
        p->addProperty(std::string("enable-icp"), std::string("0"));

    accept();
}
