#include "scatter_point_glyph.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ScatterPointGlyph w;

    QFont font;
    font.setFamily("aria");
    font.setPixelSize(9);

	w.showMaximized();
	return a.exec();
}
