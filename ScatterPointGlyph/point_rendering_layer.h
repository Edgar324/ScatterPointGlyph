#ifndef POINT_RENDERING_LAYER_H_
#define POINT_RENDERING_LAYER_H_

#include "basic_glyph_layer.h"
#include <vector>

class PointRenderingLayer : public BasicGlyphLayer
{
public:
	PointRenderingLayer();
	~PointRenderingLayer();

	void SetData(vtkPolyData* data);
	void SetHighlightPointIndex(std::vector< int >& index);

private:
	std::vector< int > point_index_;

	void UpdatePointGlyph();
};

#endif