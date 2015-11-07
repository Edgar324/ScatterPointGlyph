#include "scatter_point_glyph.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ScatterPointGlyph w;
	w.show();
	return a.exec();
}
