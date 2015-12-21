#ifndef BASIC_GLYPH_LAYER_H_
#define BASIC_GLYPH_LAYER_H_

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>

class BasicGlyphLayer
{
public:
	BasicGlyphLayer(vtkRenderer* parent = 0) {
		actor_ = vtkActor::New();
		mapper_ = vtkPolyDataMapper::New();
		poly_data_ = NULL;
	}

	~BasicGlyphLayer() {};

	vtkActor* actor() { return actor_; }

protected:
	vtkRenderer* renderer_;
	vtkActor* actor_;
	vtkPolyDataMapper* mapper_;
	vtkPolyData* poly_data_;
};

#endif