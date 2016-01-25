#include "scatter_point_glyph.h"
#include <fstream>
#include <QVTKWidget.h>
#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkDataObject.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkWriter.h>
#include <vtkCamera.h>
#include <vtkMatrix4x4.h>
#include <vtkInteractorObserver.h>
#include <vtkDelaunay2D.h>

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QActionGroup>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QSlider>

#include "point_rendering_layer.h"
#include "scatter_point_dataset.h"
#include "scatter_point_view.h"
#include "hierarchical_tree.h"
#include "multi_label_tree.h"
#include "ncut_tree.h"
#include "parallel_coordinate.h"
#include "transmap.h"
#include "transmap_data.h"
#include "path_explore_widget.h"
#include "path_dataset.h"
#include "tour_path_generator.h"
#include "tree_map_view.h"
#include "variable_selection_dialog.h"

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent), scatter_point_dataset_(NULL), sys_mode_(MULTI_LABEL_MODE), cluster_tree_(NULL),
	current_view_level_(0) {

	ui_.setupUi(this);

	this->InitWidget();
}

ScatterPointGlyph::~ScatterPointGlyph() {

}

void ScatterPointGlyph::InitWidget() {
	this->setDockOptions(QMainWindow::AllowNestedDocks);

	main_view_ = new ScatterPointView;
	main_view_->setFocusPolicy(Qt::StrongFocus);

	level_name_label_ = new QLabel("Level: ");
	level_name_label_->setFixedWidth(100);
	level_name_label_->setAlignment(Qt::AlignRight);
	level_slider_ = new QSlider(Qt::Horizontal);
	level_slider_->setMaximumWidth(600);
	level_slider_->setRange(0, 0);
	level_slider_->setSingleStep(1);
	level_slider_->setValue(0);
	level_index_label_ = new QLabel("0");
	level_index_label_->setFixedWidth(50);
	QHBoxLayout* level_layout = new QHBoxLayout;
	level_layout->addWidget(level_name_label_, Qt::AlignLeft);
	level_layout->addWidget(level_slider_, Qt::AlignLeft);
	level_layout->addWidget(level_index_label_, Qt::AlignLeft);
	connect(level_slider_, SIGNAL(valueChanged(int)), this, SLOT(OnViewLevelChanged()));

	parallel_coordinate_ = new ParallelCoordinate;
	parallel_coordinate_->setMinimumHeight(200);
	parallel_dataset_ = new ParallelDataset;

	parallel_coordinate_panel_ = new QDockWidget(QString("Parallel Coordinate"), this);
	parallel_coordinate_panel_->setWidget(parallel_coordinate_);
	this->addDockWidget(Qt::BottomDockWidgetArea, parallel_coordinate_panel_);
	ui_.menuView->addAction(parallel_coordinate_panel_->toggleViewAction());
	parallel_coordinate_panel_->setVisible(false);

	tree_map_view_ = new TreeMapView;
	tree_map_view_->setMinimumWidth(700);
	tree_map_panel_ = new QDockWidget(QString("Tree Map and Table Lens"), this);
	tree_map_panel_->setWidget(tree_map_view_);
	this->addDockWidget(Qt::RightDockWidgetArea, tree_map_panel_);
	ui_.menuView->addAction(tree_map_panel_->toggleViewAction());
	tree_map_panel_->setVisible(false);

	path_explore_view_ = new PathExploreWidget;
	path_explore_view_->setMinimumWidth(700);
	path_dataset = new PathDataset;

	path_explore_panel_ = new QDockWidget(QString("Path Records"), this);
	path_explore_panel_->setWidget(path_explore_view_);
	this->tabifyDockWidget(tree_map_panel_, path_explore_panel_);
	ui_.menuView->addAction(path_explore_panel_->toggleViewAction());
	path_explore_panel_->setVisible(false);

	QVBoxLayout* main_layout = new QVBoxLayout;
	main_layout->addWidget(main_view_);
	main_layout->addLayout(level_layout);
	ui_.centralWidget->setLayout(main_layout);

	main_renderer_ = vtkRenderer::New();
	main_view_->GetRenderWindow()->AddRenderer(main_renderer_);
	main_renderer_->SetViewport(0.0, 0.0, 1.0, 1.0);
	main_renderer_->SetBackground(1.0, 1.0, 1.0);
	connect(main_view_, SIGNAL(ViewUpdated()), this, SLOT(OnMainViewUpdated()));
	connect(main_view_, SIGNAL(GlyphSelected(int, int)), this, SLOT(OnGlyphSelected(int, int)));
	connect(main_view_, SIGNAL(LeftButtonUp()), this, SLOT(OnMainviewLeftButtonUp()));
	connect(main_view_, SIGNAL(MouseDrag(int, int)), this, SLOT(OnMouseDragmove(int, int)));

	trans_map_ = new TransMap(main_view_);
	trans_map_->SetInteractor(main_view_->GetInteractor());
	trans_map_->SetDefaultRenderer(this->main_renderer_);
	transmap_data_ = new TransMapData;

	original_point_rendering_layer_ = new PointRenderingLayer;
	original_point_rendering_layer_->SetInteractor(main_view_->GetInteractor());
	original_point_rendering_layer_->SetDefaultRenderer(this->main_renderer_);
	original_point_rendering_layer_->SetEnabled(true);

	sys_mode_action_group_ = new QActionGroup(this);
	sys_mode_action_group_->addAction(ui_.action_hierarchical_clustering);
	sys_mode_action_group_->addAction(ui_.actionChameleon_Clustering);
	sys_mode_action_group_->addAction(ui_.actionNCuts);
	sys_mode_action_group_->addAction(ui_.actionMulti_Label);
	sys_mode_action_group_->setExclusive(true);
	ui_.actionMulti_Label->setChecked(true);
	connect(sys_mode_action_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnSysmodeChanged()));

	main_view_interaction_mode_group_ = new QActionGroup(this);
	main_view_interaction_mode_group_->addAction(ui_.actionSingle_Selection);
	main_view_interaction_mode_group_->addAction(ui_.actionSequence_Selection);
	main_view_interaction_mode_group_->addAction(ui_.actionSelect_Minimum_Path);
	main_view_interaction_mode_group_->addAction(ui_.actionBrush_Path_Sequence);
	main_view_interaction_mode_group_->addAction(ui_.actionBrush_Cluster);
	main_view_interaction_mode_group_->setExclusive(true);
	connect(main_view_interaction_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnMainViewInteractionModeChanged()));

	transmap_tip_mode_group_ = new QActionGroup(this);
	transmap_tip_mode_group_->addAction(ui_.actionOff);
	transmap_tip_mode_group_->setExclusive(true);
	connect(transmap_tip_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnShowVarTrendTriggered()));

	// load data actions
	connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));

	// actions for tips on the cluster transition map
	ui_.mainToolBar->insertAction(ui_.actionShow_Minimum_Spanning_Tree, ui_.menuShow_Sequence->menuAction());
	connect(ui_.actionShow_Minimum_Spanning_Tree, SIGNAL(triggered()), this, SLOT(OnShowMstTriggered()));

	// actions for manipulating cluster nodes
	connect(ui_.actionExec, SIGNAL(triggered()), this, SLOT(OnExecClusteringTriggered()));
	connect(ui_.actionSplit, SIGNAL(triggered()), this, SLOT(OnSplitClusterTriggered()));
	connect(ui_.actionMerge, SIGNAL(triggered()), this, SLOT(OnMergeClusterTriggered()));

	// action for saving exploration path
	connect(ui_.actionAdd_Path_Sequence, SIGNAL(triggered()), this, SLOT(OnSavePathSequenceTriggered()));
}

void ScatterPointGlyph::OnActionOpenVtkFileTriggered() {
	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName("./TestData/timestep0_0003600.vtk");
	reader->SetReadAllScalars(1);
	reader->SetReadAllVectors(1);
	reader->Update();
	vtkDataObject* data = reader->GetOutput();
	vtkUnstructuredGrid* scatter_point_data = vtkUnstructuredGrid::SafeDownCast(data);

	vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
	polydata->SetPoints(scatter_point_data->GetPoints());

	vtkIntArray* id_array = vtkIntArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(0));
	vtkFloatArray* speed_array = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(1));
	vtkFloatArray* xparray = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(2));
	vtkFloatArray* yparray = vtkFloatArray::SafeDownCast(scatter_point_data->GetPointData()->GetArray(3));

	float x_scale1 = 0.35, x_scale2 = 0.48;
	float y_scale1 = 0.65, y_scale2 = 0.77;
	double bounds[6];
	scatter_point_data->GetBounds(bounds);

	for (int i = 0; i < scatter_point_data->GetPoints()->GetNumberOfPoints(); ++i) {
		float x = scatter_point_data->GetPoint(i)[0];
		float y = scatter_point_data->GetPoint(i)[1];
		if (x > bounds[0] + (bounds[1] - bounds[0]) * x_scale1 && x < bounds[0] + (bounds[1] - bounds[0]) * x_scale2
			&& y > bounds[2] + (bounds[3] - bounds[2]) * y_scale1 && y < bounds[2] + (bounds[3] - bounds[2]) * y_scale2) {
			std::vector< float > pos;
			std::vector< float > value;
			pos.push_back(x);
			pos.push_back(y);
			value.push_back(speed_array->GetValue(i));
			value.push_back(xparray->GetValue(i));
			value.push_back(yparray->GetValue(i));

			scatter_point_dataset_->original_point_pos.push_back(pos);
			scatter_point_dataset_->original_point_values.push_back(value);
		}
	}

	scatter_point_dataset_->DirectConstruct();
	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionOpenScatterFileTriggered() {
	//QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.sc");
	//if (file_path.length() == 0) return;

	if (scatter_point_dataset_ == NULL) scatter_point_dataset_ = new ScatterPointDataset;
	scatter_point_dataset_->ClearData();

	//std::ifstream input_file(file_path.toLocal8Bit());
	std::ifstream input_file("./TestData/wine.sc");
	char char_str[1000];
	input_file.getline(char_str, 1000);
	QString value_str = QString::fromLocal8Bit(char_str);
	QStringList value_list = value_str.split(' ');

	int record_num, var_num;
	record_num = value_list.at(0).toInt();
	var_num = value_list.at(1).toInt();

	input_file.getline(char_str, 1000);
	value_str = QString::fromLocal8Bit(char_str);
	value_list = value_str.split(' ');
	for (int i = 0; i < value_list.size(); ++i) scatter_point_dataset_->var_names.push_back(value_list.at(i));

	VariableSelectionDialog var_dialog;
	var_dialog.SetDatasetInfo(record_num, var_num, scatter_point_dataset_->var_names);
	if (var_dialog.exec() != QDialog::Accepted) {
		input_file.close();
		return;
	}

	scatter_point_dataset_->original_point_pos.resize(record_num);
	scatter_point_dataset_->original_point_values.resize(record_num);
	for (int i = 0; i < record_num; ++i) {
		scatter_point_dataset_->original_point_values[i].resize(var_num);

		input_file.getline(char_str, 1000);
		value_str = QString::fromLocal8Bit(char_str);
		value_list = value_str.split(',');

		for (int j = 0; j < var_num; ++j)
			scatter_point_dataset_->original_point_values[i][j] = value_list.at(j).toFloat();
	}
	input_file.close();

	if (var_dialog.IsAutomaticDimReduction()) {
		int dim_num = var_dialog.GetDimNumber();
		if (dim_num <= 0 && dim_num > var_num) return;
		scatter_point_dataset_->AutoDimReduction(dim_num);
	} else {
		std::vector< float > dim_weights;
		std::vector< bool > is_dim_selected;
		var_dialog.GetSelectionResult(dim_weights);
		scatter_point_dataset_->var_weights.clear();
		for (int i = 0; i < dim_weights.size(); ++i)
			if (dim_weights[i] >= 0) {
				scatter_point_dataset_->var_weights.push_back(dim_weights[i]);
				is_dim_selected.push_back(true);
			}
			else
				is_dim_selected.push_back(false);
		scatter_point_dataset_->ManualSelectDim(is_dim_selected);
	}

	scatter_point_dataset_->ExecMds();
	scatter_point_dataset_->DirectConstruct();

	this->UpdateMenus();

	this->AddPointData2View();
}

void ScatterPointGlyph::OnActionCloseTriggered() {

}

void ScatterPointGlyph::OnActionExitTriggered() {
	
}

void ScatterPointGlyph::OnSysmodeChanged() {
	SystemMode current_mode;
	if (ui_.action_hierarchical_clustering->isChecked())
		current_mode = HIER_MODE;
	else if (ui_.actionChameleon_Clustering->isChecked())
		current_mode = CHAMELEON_MODE;
	else if (ui_.actionNCuts->isChecked())
		current_mode = NCUTS_MODE;
	else if (ui_.actionMulti_Label)
		current_mode = MULTI_LABEL_MODE;

	if (current_mode != sys_mode_) {
		if (cluster_tree_ != NULL) {
			delete cluster_tree_;
			cluster_tree_ = NULL;
		}

		sys_mode_ = current_mode;
	}
}

void ScatterPointGlyph::OnExecClusteringTriggered() {
	if (scatter_point_dataset_ == NULL) {
		QMessageBox::information(this, tr("Warning"), tr("Please load data first."));
		return;
	}

	float dis_per_pixel = this->GetMainViewDisPerPixel();
	
	switch (sys_mode_)
	{
	case ScatterPointGlyph::HIER_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new HierarchicalTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		HierarchicalTree* hier_tree = dynamic_cast<HierarchicalTree*>(cluster_tree_);
		if (hier_tree != NULL) {
			hier_tree->SetExpectedClusterNum(8);
			hier_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::CHAMELEON_MODE:
	{
		
	}
		break;
	case ScatterPointGlyph::NCUTS_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new NCutTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		NCutTree* ncut_tree = dynamic_cast<NCutTree*>(cluster_tree_);
		if (ncut_tree != NULL) {
			ncut_tree->start();
		}
	}
		break;
	case ScatterPointGlyph::MULTI_LABEL_MODE:
	{
		if (cluster_tree_ == NULL) {
			cluster_tree_ = new MultiLabelTree(scatter_point_dataset_);

			connect(cluster_tree_, SIGNAL(finished()), this, SLOT(OnClusterFinished()));
		}

		MultiLabelTree* un_tree = dynamic_cast< MultiLabelTree* >(cluster_tree_);
		if (un_tree != NULL) {
			float dis_per_pixel = this->GetMainViewDisPerPixel();
			un_tree->SetRadiusThreshold(200.0 * dis_per_pixel / (scatter_point_dataset_->original_pos_ranges[0][1] - scatter_point_dataset_->original_pos_ranges[0][0]));
			un_tree->start();
		}
	}
		break;
	default:
		break;
	}
}

void ScatterPointGlyph::AddPointData2View() {
	original_point_rendering_layer_->SetData(scatter_point_dataset_);

	main_renderer_->ResetCamera();

	this->UpdateParallelCoordinate();

	main_view_->update();
}

void ScatterPointGlyph::OnClusterFinished() {
	int max_level = cluster_tree_->GetMaxLevel();
	level_slider_->setRange(0, max_level);
	level_slider_->setValue(1);
	level_index_label_->setText(QString("%0").arg(1));
	level_name_label_->setText(QString("Level(0~%0): ").arg(max_level));

	current_view_level_ = 1;

	this->UpdateAllViews();
}

void ScatterPointGlyph::OnMainViewUpdated() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* multi_label_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);
		if (multi_label_tree == NULL) return;

		float dis_per_pixel = this->GetMainViewDisPerPixel();
		current_view_level_ = multi_label_tree->GetRadiusLevel(label_pixel_radius_ * dis_per_pixel);
		level_slider_->setValue(current_view_level_);
		level_index_label_->setText(QString("%0").arg(current_view_level_));
	}

	this->UpdateAllViews();
}

void ScatterPointGlyph::OnViewLevelChanged() {
	current_view_level_ = level_slider_->value();
	level_index_label_->setText(QString("%0").arg(level_slider_->value()));

	this->UpdateAllViews();
}

void ScatterPointGlyph::UpdateAllViews() {
	this->UpdateTransmap();
	this->UpdateParallelCoordinate();
	this->UpdateTreemap();
}

void ScatterPointGlyph::UpdateParallelCoordinate() {
	if (transmap_data_->cluster_nodes.size() == 0) {
		parallel_dataset_->ClearData();
		parallel_dataset_->subset_names.push_back(QString("MDS Data"));
		parallel_dataset_->subset_records.resize(1);
		parallel_dataset_->axis_anchors.resize(scatter_point_dataset_->original_point_values[0].size());
		for (int i = 0; i < scatter_point_dataset_->original_point_values[0].size(); ++i) {
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(scatter_point_dataset_->original_value_ranges[i][0]));
			parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(scatter_point_dataset_->original_value_ranges[i][1]));
			parallel_dataset_->axis_names.push_back(scatter_point_dataset_->var_names[i]);
		}

		for (int i = 0; i < scatter_point_dataset_->point_values.size(); ++i) {
			ParallelRecord* record = new ParallelRecord;
			record->values = scatter_point_dataset_->point_values[i];
			parallel_dataset_->subset_records[0].push_back(record);
		}
		parallel_dataset_->CompleteInput();
		parallel_dataset_->UpdateGaussian();
		parallel_coordinate_->SetDataset(parallel_dataset_);

		QString str = QString("%0 points.").arg(scatter_point_dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	} else {
		std::vector< int > selection_index;
		trans_map_->GetSelectedClusterIndex(selection_index);
		if (selection_index.size() == 0) {
			for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i) selection_index.push_back(i);
		}
		this->GenerateParallelDataset(parallel_dataset_, selection_index);
		parallel_coordinate_->SetDataset(parallel_dataset_);
		parallel_coordinate_->update();

		int selected_count = 0;
		for (int i = 0; i < parallel_dataset_->subset_records.size(); ++i)
			selected_count += parallel_dataset_->subset_records[i].size();
		QString str = QString("%0 out of %1 (%2%) points are selected.").arg(selected_count).arg(scatter_point_dataset_->original_point_pos.size()).arg((int)selected_count * 100 / scatter_point_dataset_->original_point_pos.size());
		ui_.statusBar->showMessage(str);
	}
}

void ScatterPointGlyph::GenerateParallelDataset(ParallelDataset* pdata, std::vector< int >& cluster_ids) {
	pdata->is_updating = true;
	pdata->ClearData();

	std::vector< int > cluster_color = original_point_rendering_layer_->GetClusterColor();
	parallel_dataset_->subset_names.resize(cluster_ids.size());
	parallel_dataset_->subset_colors.resize(cluster_ids.size());
	parallel_dataset_->subset_records.resize(cluster_ids.size());
	parallel_dataset_->axis_names.resize(scatter_point_dataset_->original_value_ranges.size());
	parallel_dataset_->axis_anchors.resize(scatter_point_dataset_->original_value_ranges.size());

	for (int i = 0; i < cluster_ids.size(); ++i) {
		parallel_dataset_->subset_names[i] = QString("Cluster %0").arg(i);
		parallel_dataset_->subset_colors[i] = QColor(cluster_color[3 * cluster_ids[i]], cluster_color[3 * cluster_ids[i] + 1], cluster_color[3 * cluster_ids[i] + 2]);
		for (int j = 0; j < scatter_point_dataset_->point_num; ++j)
			if (cluster_index[j] == cluster_ids[i]) {
				ParallelRecord* record = new ParallelRecord;
				record->values = scatter_point_dataset_->point_values[j];
				parallel_dataset_->subset_records[i].push_back(record);
			}
	}
	for (int i = 0; i < scatter_point_dataset_->original_value_ranges.size(); ++i) {
		parallel_dataset_->axis_names[i] = scatter_point_dataset_->var_names[i];
		parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(scatter_point_dataset_->original_value_ranges[i][0]));
		parallel_dataset_->axis_anchors[i].push_back(QString("%0").arg(scatter_point_dataset_->original_value_ranges[i][1]));
	}
	parallel_dataset_->CompleteInput();
	parallel_dataset_->UpdateGaussian();

	pdata->is_updating = false;
}



void ScatterPointGlyph::UpdateTransmap() {
	if (transmap_data_ != NULL) {
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			transmap_data_->cluster_nodes[i]->is_highlighted = false;
		transmap_data_->cluster_nodes.clear();
	}

	float dis_per_pixel = this->GetMainViewDisPerPixel();
	// dis_per_pixel * 100.0
	cluster_tree_->GetClusterResult(current_view_level_, transmap_data_->cluster_nodes);
	cluster_tree_->GetClusterResult(current_view_level_, cluster_num, cluster_index);
	
	std::vector< QColor > colors;
	for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
		colors.push_back(transmap_data_->cluster_nodes[i]->color);

	original_point_rendering_layer_->SetClusterIndex(cluster_num, cluster_index, colors);
	original_point_rendering_layer_->SetHighlightCluster(-1);

	transmap_data_->ClearData();
	transmap_data_->cluster_num = cluster_num;
	transmap_data_->var_num = scatter_point_dataset_->var_weights.size();
	transmap_data_->dataset = scatter_point_dataset_;
	transmap_data_->ProcessData();

	trans_map_->SetNodeRadius(dis_per_pixel * 50);
	trans_map_->SetData(scatter_point_dataset_, transmap_data_);

	main_view_->update();
}


void ScatterPointGlyph::UpdatePathMap() {

}

void ScatterPointGlyph::UpdateTreemap() {
	std::vector< int > selected_ids;
	trans_map_->GetSelectedClusterIds(selected_ids);

	//UncertaintyTree* un_tree = dynamic_cast<UncertaintyTree*>(cluster_tree_);
	cluster_tree_->SortTree(selected_ids);

	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	std::vector< bool > is_selected;
	is_selected.resize(transmap_data_->cluster_nodes.size(), false);

	std::vector< int > var_order = parallel_coordinate_->GetAxisOrder();
	trans_map_->SetAxisOrder(var_order);

	std::vector< CNode* > selected_nodes;
	if (selected_ids.size() != 0) {
		for (int i = 0; i < selection_index.size(); ++i) {
			selected_nodes.push_back(transmap_data_->cluster_nodes[selection_index[i]]);
			is_selected[selection_index[i]] = true;
		}
		for (int i = 0; i < transmap_data_->cluster_nodes.size(); ++i)
			if (!is_selected[i]) selected_nodes.push_back(transmap_data_->cluster_nodes[i]);
		tree_map_view_->SetData(cluster_tree_->root(), scatter_point_dataset_->var_num, selected_nodes, selected_ids.size(), var_order, scatter_point_dataset_->var_names);
	}
	else {
		selected_nodes = transmap_data_->cluster_nodes;
		tree_map_view_->SetData(cluster_tree_->root(), scatter_point_dataset_->var_num, selected_nodes, selected_nodes.size(), var_order, scatter_point_dataset_->var_names);
	}

	tree_map_view_->scene()->update();
}

float ScatterPointGlyph::GetMainViewDisPerPixel() {
	double point_one[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 1, 0, 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeWorldToDisplay(this->main_renderer_, 2, 0, 0, point_two);
	return (1.0 / abs(point_one[0] - point_two[0]));
}

void ScatterPointGlyph::GetSceneRange(float& left, float& right, float& bottom, float& top) {
	double viewport[4];
	this->main_renderer_->GetViewport(viewport);
	double point_one[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[0] * this->main_view_->width(), viewport[1] * this->main_view_->height(), 0, point_one);
	double point_two[4];
	vtkInteractorObserver::ComputeDisplayToWorld(this->main_renderer_, viewport[2] * this->main_view_->width(), viewport[3] * this->main_view_->height(), 0, point_two);
	left = point_one[0];
	right = point_two[0];
	bottom = point_one[1];
	top = point_two[1];
}

void ScatterPointGlyph::OnGlyphSelected(int x, int y) {
	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	UpdateParallelCoordinate();
	UpdateTreemap();
}

void ScatterPointGlyph::OnMainviewLeftButtonUp() {
	if (trans_map_ != NULL) {
		trans_map_->OnMouseReleased();
		if (ui_.actionBrush_Cluster->isChecked() || ui_.actionBrush_Path_Sequence->isChecked()) {
			std::vector< int > selection_index;
			trans_map_->GetSelectedClusterIndex(selection_index);
			if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
				original_point_rendering_layer_->SetHighlightClusters(selection_index);
			}

			UpdateParallelCoordinate();
			UpdateTreemap();

			this->main_view_->update();
		}
	}
}

void ScatterPointGlyph::OnMainViewRightButtonDown() {
}

void ScatterPointGlyph::OnMouseDragmove(int x, int y) {
	if (trans_map_ != NULL) trans_map_->OnMouseMove(x, y);
}

void ScatterPointGlyph::OnSplitClusterTriggered() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* un_tree = dynamic_cast< MultiLabelTree* >(cluster_tree_);

		int temp_cluster_index = trans_map_->GetSelectedClusterIndex();
		if (temp_cluster_index != -1) {
			un_tree->SplitCluster(transmap_data_->cluster_nodes[temp_cluster_index]->id());
			UpdateTransmap();
			UpdateTreemap();
			UpdateParallelCoordinate();
		}
	}
}

void ScatterPointGlyph::OnMergeClusterTriggered() {
	if (sys_mode_ == MULTI_LABEL_MODE && cluster_tree_ != NULL) {
		MultiLabelTree* un_tree = dynamic_cast<MultiLabelTree*>(cluster_tree_);

		std::vector< int > merged_clusters;
		trans_map_->GetSelectedClusterIndex(merged_clusters);

		std::vector< int > cluster_ids;
		for (int i = 0; i < merged_clusters.size(); ++i)
			cluster_ids.push_back(transmap_data_->cluster_nodes[merged_clusters[i]]->id());
		un_tree->MergeClusters(cluster_ids);

		UpdateTransmap();
		UpdateTreemap();
		UpdateParallelCoordinate();
	}
}

void ScatterPointGlyph::OnTreemapNodeSelected(int node_id) {
	trans_map_->OnNodeSelected(node_id);

	std::vector< int > selection_index;
	trans_map_->GetSelectedClusterIndex(selection_index);
	if (original_point_rendering_layer_ != NULL && trans_map_ != NULL) {
		original_point_rendering_layer_->SetHighlightClusters(selection_index);
	}

	this->GenerateParallelDataset(parallel_dataset_, selection_index);

	parallel_coordinate_->update();

	this->main_view_->update();
}

void ScatterPointGlyph::OnMainViewInteractionModeChanged() {
	if (ui_.actionSingle_Selection->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_SINGLE_CLUSTER);
	} else if (ui_.actionSequence_Selection->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_MULTI_CLUSTERS);
	} else if (ui_.actionSelect_Minimum_Path->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_MINIMUM_PATH);
	} else if (ui_.actionBrush_Path_Sequence->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_BRUSHED_PATH_SEQUENCE);
	} else if (ui_.actionBrush_Cluster->isChecked()) {
		trans_map_->SetInteractionState(TransMap::SELECT_BRUSHED_CLUSTERS);
	} else {
		trans_map_->SetInteractionState(TransMap::NORMAL);
	}
}

void ScatterPointGlyph::OnSavePathSequenceTriggered() {
	std::list< CNode* > node_seq = trans_map_->GetNodeSequence();

	PathRecord* record = new PathRecord;
	record->item_values.resize(node_seq.size());
	record->item_color.resize(node_seq.size());
	std::list< CNode* >::iterator node_iter = node_seq.begin();

	int count = 0;
	while (node_iter != node_seq.end()) {
		record->item_values[count] = (*node_iter)->average_values;
		record->item_color[count] = QColor(128, 128, 128);
		count++;
		node_iter++;
	}

	record->change_values.resize(node_seq.size() - 1);
	for (int i = 0; i < node_seq.size() - 1; ++i) {
		record->change_values[i].resize(record->item_values[0].size());
		for (int j = 0; j < record->item_values[0].size(); ++j)
			record->change_values[i][j] = record->item_values[i + 1][j] - record->item_values[i][j];
	}
	for (int k = 0; k < record->item_values[0].size(); ++k)
		record->var_names.push_back("Test");

	path_dataset->path_records.push_back(record);
	path_explore_view_->SetData(path_dataset);
}

void ScatterPointGlyph::OnShowMstTriggered() {
	if (ui_.actionShow_Minimum_Spanning_Tree->isChecked()) {
		//ui_.actionShow_Sequence->setChecked(false);

		trans_map_->ShowMinimumSpanningTree(true);
	} else {
		trans_map_->ShowMinimumSpanningTree(false);
	}

	main_view_->update();
}

void ScatterPointGlyph::OnShowVarTrendTriggered() {
	if (!ui_.actionOff->isChecked()) {
		ui_.actionShow_Minimum_Spanning_Tree->setChecked(false);

		int var_index = -1;
		QList< QAction* > actions = transmap_tip_mode_group_->actions();
		for (int i = 0; i < actions.size(); ++i)
			if (actions.at(i)->isChecked()) {
				var_index = i;
				break;
			}

		trans_map_->ShowVarTrend(var_index - 1);
	} else {
		trans_map_->ShowVarTrend(-1);
	}

	main_view_->update();
}

void ScatterPointGlyph::UpdateMenus() {
	QList< QAction* > actions = transmap_tip_mode_group_->actions();
	for (int i = 1; i < actions.size(); ++i) {
		transmap_tip_mode_group_->removeAction(actions.at(i));
		ui_.menuShow_Sequence->removeAction(actions.at(i));
		delete actions.at(i);
	}

	ui_.actionOff->setChecked(true);
	for (int i = 0; i < scatter_point_dataset_->var_names.size(); ++i) {
		QAction* action = ui_.menuShow_Sequence->addAction(scatter_point_dataset_->var_names[i]);
		action->setCheckable(true);
		action->setChecked(false);
		transmap_tip_mode_group_->addAction(action);
	}
	//transmap_tip_mode_group_->setExclusive(true);
}