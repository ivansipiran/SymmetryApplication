#ifndef ROTATIONPARAMETERDIALOG_H
#define ROTATIONPARAMETERDIALOG_H

#include <QDialog>

namespace Ui {
class RotationParameterDialog;
}

class RotationParameterDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RotationParameterDialog(QWidget *parent = 0);
    ~RotationParameterDialog();

private:
    Ui::RotationParameterDialog *ui;
};

#endif // ROTATIONPARAMETERDIALOG_H
