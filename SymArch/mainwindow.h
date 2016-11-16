#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QListWidget>
#include "glwidget.h"
#include "objectcontainer.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private Q_SLOTS:
    void about();
    void open();
    void deleteMesh();
    void mergeMeshes();
    void repairRotationalSymmetry();
    void segmentateMesh();
    void saveMesh();

private:

    void createActions();
    void createMenus();
    void createToolbar();
    void createStatusBar();
    void createDockWindows();

    Ui::MainWindow *ui;

    //GUI controls
    QToolBar*   m_toolbarMain;
    QMenu*      m_menuFile;
    QMenu*      m_menuProcessing;
    QMenu*      m_menuHelp;

    QAction*    m_actionOpen;
    QAction*    m_actionExit;
    QAction*    m_actionRotSymmetry;

    QAction*    m_actionAbout;
    QListWidget*    m_ListMeshes;

    bool m_bOpenFile;
    GLWidget*   m_glArea;
};

#endif // MAINWINDOW_H
