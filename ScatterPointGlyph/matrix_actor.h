#ifndef MATRIX_ACTOR_H_
#define MATRIX_ACTOR_H_

#include "vtkActor2D.h"

class MatrixActor : public vtkActor2D
{
public:
	vtkTypeMacro(MatrixActor, vtkActor2D);
	void PrintSelf(ostream& os, vtkIndent indent);

	static MatrixActor* New();

	// Support the standard render methods.
	virtual int RenderOverlay(vtkViewport *viewport);
	virtual int RenderOpaqueGeometry(vtkViewport *viewport);
	virtual int RenderTranslucentPolygonalGeometry(vtkViewport *viewport);

	// Description:
	// Does this prop have some translucent polygonal geometry?
	virtual int HasTranslucentPolygonalGeometry();

	// Description:
	// Set/Get the vtkMapper2D which defines the data to be drawn.
	virtual void SetMapper(vtkMapper2D *mapper);
	vtkGetObjectMacro(Mapper, vtkMapper2D);

protected:
	MatrixActor();
	~MatrixActor();

	vtkTextActor* title_actor;
	std::vector< vtkTextActor* > x_labels;
	std::vector< vtkTextActor* > y_labels;

	vtkTextProperty* title_text_property;
	vtkTextProperty* label_text_property;

	std::vector< std::string > x_label_names;
	std::vector< std::string > y_label_names;
	std::vector< std::vector< float > > matrix_values;


private:
	MatrixActor(const MatrixActor&);  // Not implemented.
	void operator=(const MatrixActor&);  // Not implemented.

};

#endif