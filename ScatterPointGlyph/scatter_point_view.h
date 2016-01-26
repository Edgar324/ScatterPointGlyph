#ifndef SCATTER_POINT_VIEW_H_
#define SCATTER_POINT_VIEW_H_

#include <QVTKWidget.h>

class QTimer;
class CNode;

class ScatterPointView : public QVTKWidget
{
	Q_OBJECT

public:
	ScatterPointView();
	~ScatterPointView();

	void ShowTooltip(int point_count, QString axis_name, float average_value, float variance_value);
	void HideTooltip();
	void SetHighlightVarIndex(int var_index);

signals:
	void ViewUpdated();
	void GlyphSelected(int x, int y);
	void LeftButtonUp();
	void RightButtonDown();
	void MouseDrag(int x, int y);
	void HighlightVarChanged(int);

protected:
	void wheelEvent(QWheelEvent* event);
	void mousePressEvent(QMouseEvent* event);
	void mouseMoveEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);

private:
	QTimer* timer_;
	bool is_wheel_updated_;
	QPoint global_mouse_pos_;

	private slots:
		void OnTimeout();
};

#endif