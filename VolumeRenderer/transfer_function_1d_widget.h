#ifndef TRANSFER_FUNCTION_1D_WIDGET_H
#define TRANSFER_FUNCTION_1D_WIDGET_H

#include <QWidget>
#include <QPolygonF>
#include <vector>
#include <list>
using std::vector;
using std::list;

#include "volume_renderer_global.h"

class QRenderHoverPoints;

class VOLUMERENDERER_EXPORT TransferFunction1DWidget : public QWidget
{
	Q_OBJECT

public:
	TransferFunction1DWidget();
	~TransferFunction1DWidget();

	enum DataValueType { TYPE_UNKNOWN = 0, TYPE_UNSIGNED_BYTE, TYPE_SHORT, TYPE_FLOAT };
	enum TfLineType { TYPE_LINE = 0, TYPE_CURVE };

	void SetLineType( TfLineType type );

	void SetData( unsigned int data_size, short* data );
	void SetData( unsigned int data_size, unsigned char* data );
    void SetData( unsigned int data_size, float* data );
	void SetViewWindow( double window_level, double window_width );
	
	void GetTfColorAndAlpha( int tf_size, vector< float >& tf_values );
	void GetTfColor( int tf_size, vector< float >& tf_values );
	void GetTfAlpha( int tf_size, vector< float >& tf_values );

    float min_value() { return min_value_; }
    float max_value() { return max_value_; }

	public slots:
		void OnViewWindowChanged( double window_level, double window_width );
		void OnHoverPointChanged();

signals:
	void ValueChanged();

protected:
	void paintEvent( QPaintEvent* event );
	void resizeEvent( QResizeEvent* event );

private:
	float min_value_, max_value_, view_min_value_, view_max_value_;
	vector< unsigned int > histogram_values_;
	int histogram_bin_size_;
	unsigned int max_histogram_height_;

	vector< float > tf_value_array_;

	void* data_;
	unsigned int data_size_;
	DataValueType data_type_;

	QRenderHoverPoints* render_hover_points_;

	TfLineType line_type_;

	void DrawHistogram();
	void DrawTfColorAndAlphaBar();
	void DrawScalePlate();

	template < typename T >
	void GenerateHistogramValues( unsigned int data_size, T* data );
    template < typename T >
    void GetMinMaxValue(unsigned int data_size, T* data, float& min_value, float& max_value);

	void InitWidget();
	void UpdateTfValueArray();
	void GenerateLineColorAndAlphaValue();
	void GenerateCurveColorAndAlphaValue();
};

#endif // TRANSFER_FUNCTION_1D_WIDGET_H
