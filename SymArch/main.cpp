#include "mainwindow.h"
#include <QApplication>
#include <QDesktopWidget>
#include <QSurfaceFormat>
#include <iostream>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QSurfaceFormat fmt;
    fmt.setDepthBufferSize(24);
    fmt.setSamples(4);
    fmt.setVersion(4,5);
    fmt.setProfile(QSurfaceFormat::CoreProfile);

    std::cout << fmt.version().first << " " << fmt.version().second << std::endl;

    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow mainWindow;
    mainWindow.resize(mainWindow.sizeHint());
    mainWindow.showMaximized();
    return app.exec();
}
