#include <QtGui/QApplication>
#include "chipboard.h"


int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    chipboard chipboard;
    chipboard.show();
    return app.exec();
}
