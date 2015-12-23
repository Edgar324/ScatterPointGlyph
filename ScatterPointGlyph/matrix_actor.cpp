#include "matrix_actor.h"

#include "vtkAlgorithm.h"
#include "vtkAlgorithmOutput.h"
#include "vtkAppendPolyData.h"
#include "vtkAxisActor2D.h"
#include "vtkCellArray.h"
#include "vtkDataObjectCollection.h"
#include "vtkDataSetCollection.h"
#include "vtkFieldData.h"
#include "vtkDoubleArray.h"
#include "vtkGlyph2D.h"
#include "vtkGlyphSource2D.h"
#include "vtkIntArray.h"
#include "vtkLegendBoxActor.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPlanes.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkProperty2D.h"
#include "vtkTextActor.h"
#include "vtkTextMapper.h"
#include "vtkTextProperty.h"
#include "vtkTrivialProducer.h"
#include "vtkViewport.h"

vtkStandardNewMacro(MatrixActor);

MatrixActor::MatrixActor() {
	this->PositionCoordinate->SetCoordinateSystemToNormalizedViewport();
	this->PositionCoordinate->SetValue(.25, .25);
	this->Position2Coordinate->SetValue(.5, .5);
}