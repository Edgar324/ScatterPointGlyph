#ifndef MAP_RENDERING_LAYER_H_
#define MAP_RENDERING_LAYER_H_

#include "basic_glyph_layer.h"
#include <vector>

class vtkPolyData;

class MapRenderingLayer : public BasicGlyphLayer
{
public:
	MapRenderingLayer();
	~MapRenderingLayer();

	void SetFieldData(float* data, int w, int h, float wSpacing, float hSpacing, float startX, float startY);
	void Update();

private:
	std::vector< vtkPolyData* > iso_contours_;
	vtkPolyData* border_data_;
};

#endif