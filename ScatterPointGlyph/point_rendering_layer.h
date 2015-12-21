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
	void SetClusterIndex(int cluster_count, std::vector< int >& point_index);
	void SetHighlightPointIndex(std::vector< int >& index);
	void SetPointValue(std::vector< float >& values);
	std::vector< int >& GetClusterColor() { return cluster_color_; }

private:
	std::vector< int > point_index_;
	std::vector< int > cluster_color_;

	void UpdatePointGlyph();
};

#endif