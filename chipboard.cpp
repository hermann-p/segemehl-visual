#include "chipboard.h"
#include "plot.h"

#include <QtGui/QLabel>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QAction>
#include <QGraphicsView>

chipboard::chipboard()
{
    QLabel* label = new QLabel( this );
    
    LinearPlot lp;
    
    QGraphicsView view(&lp);
    label->setText( "Hello World!" );
//    setCentralWidget( view );
    view.show();
    QAction* action = new QAction(this);
    action->setText( "Quit" );
    connect(action, SIGNAL(triggered()), SLOT(close()) );
    menuBar()->addMenu( "File" )->addAction( action );
}

chipboard::~chipboard()
{}

#include "chipboard.moc"
