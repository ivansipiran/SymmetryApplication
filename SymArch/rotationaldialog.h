#ifndef ROTATIONALDIALOG_H
#define ROTATIONALDIALOG_H

#include <QDialog>
#include <QLineEdit>
#include <QGroupBox>
#include <QCheckBox>
#include <QPushButton>
#include "Util/PropertySet.h"

class RotationalDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RotationalDialog(QWidget *parent = 0);
    void setPropertySet(Util::PropertySet* prop){p = prop;}

public Q_SLOTS:
    void confirmDialog();

private:
    QLineEdit* openvdbParameter;
    QLineEdit* simplifyParameter;
    QLineEdit* fractionHarris;
    QLineEdit* thresholdRANSAC;
    QLineEdit* alpha;
    QLineEdit* iterNonrigid;
    QLineEdit* localParNonrigid;
    QLineEdit* angleNonRigid;
    QLineEdit* rotationsNonrigid;

    QPushButton* buttonCancel;
    QPushButton* buttonOK;

    QGroupBox* preprocesamientoGB;
    QGroupBox* seleccionGB;
    QGroupBox* alineacionGB;

    QCheckBox* checkVolumetric;
    QCheckBox* checkSimplification;
    QCheckBox* checkICP;

    Util::PropertySet* p;


};

#endif // ROTATIONALDIALOG_H
