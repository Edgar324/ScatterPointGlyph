#ifndef SCATTER_POINT_VIEW_H_
#define SCATTER_POINT_VIEW_H_

#include <QVTKWidget.h>

class QTimer;

class ScatterPointView : public QVTKWidget
{
	Q_OBJECT

public:
	ScatterPointView();
	~ScatterPointView();

signals:
	void ViewUpdated();
	void GlyphSelected(int x, int y);

protected:
	void wheelEvent(QWheelEvent* event);
	void mousePressEvent(QMouseEvent* event);
	void mouseMoveEvent(QMouseEvent* event);

private:
	QTimer* timer_;
	bool is_wheel_updated_;

	private slots:
		void OnTimeout();
};

#endif