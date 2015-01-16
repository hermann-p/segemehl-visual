#include <QtGui/QApplication>
#include "linplot.h"
#include "genome.h"
#include "utils.h"

#include <string>
#include <QGraphicsView>
#include <QPrinter>
#include <QPainter>
#include <QSize>


int main(int argc, char** argv)
{
    Genome* g(new Genome());
    std::string fileName;
    if (argc > 1) {
      fileName = argv[1];
    }
    else {
      fileName = "out.sam";
    }
    g->read(fileName);
    
    QApplication app(argc, argv);
    LinearPlot lp();
    
    auto seed = g->getReadAt(1, 1);
    assume(seed != nullptr, "Error, no read found!");
    lp.setSceneRect(0, 0, 1200, 300);
    lp.fromRead(seed, g);
    
    QGraphicsView view(&lp);
    view.setWindowTitle("Chipboard demo v0.0.1");
    view.show();
    QPrinter printer;
    printer.setPaperSize(QSize(view.width(), view.height()), QPrinter::Point);
    printer.setFullPage(true);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName("plot.pdf");
    QPainter pdfPainter;
    pdfPainter.begin(&printer);
    view.render(&pdfPainter);
    return app.exec();
}
