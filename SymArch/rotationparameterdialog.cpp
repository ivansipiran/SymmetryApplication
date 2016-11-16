#include "rotationparameterdialog.h"
#include "ui_rotationparameterdialog.h"

RotationParameterDialog::RotationParameterDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RotationParameterDialog)
{
    ui->setupUi(this);
}

RotationParameterDialog::~RotationParameterDialog()
{
    delete ui;
}
