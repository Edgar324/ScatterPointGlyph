#include "scatter_point_glyph.h"

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QActionGroup>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QSlider>
#include <QtWidgets/QTableView>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QProgressDialog>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "scatter_point_dataset.h"
#include "point_data_reader.h"
#include "geo_point_data_reader.h"

#include "hierarchical_tree.h"
#include "multi_label_tree.h"
#include "ncut_tree.h"
#include "view_dependent_tree.h"
#include "cluster_projection_tree.h"

#include "glyph_rendering_widget.h"
#include "parallel_coordinate.h"
#include "tree_map.h"
#include "var_selection_widget.h"

#include "glyph_dataset_builder.h"
#include "parallel_dataset_builder.h"
#include "scatter_volume_dataset.h"
#include "volume_data_reader.h"
#include "volume_renderer.h"
#include "table_lens.h"
#include "glyph_object.h"
#include "quality_metric.h"
#include "utility.h"

//#define USE_QUALITY_METRIC
//#define SAVE_PROJECTION
//#define USE_SAVED_PROJECTION

ScatterPointGlyph::ScatterPointGlyph(QWidget *parent)
	: QMainWindow(parent) {

	ui_.setupUi(this);

	this->InitWidget();
}

ScatterPointGlyph::~ScatterPointGlyph() {
    
}

void ScatterPointGlyph::InitWidget() {
	this->setDockOptions(QMainWindow::AllowNestedDocks);

    glyph_widget_ = new GlyphRenderingWidget;
    glyph_dataset_ = new GlyphDataset;
    
    projection_control_widget_ = new QWidget;
    ui_projection_control_.setupUi(projection_control_widget_);
    ui_projection_control_.axisSelectionWidget->setVisible(false);

    QVBoxLayout* central_layout = new QVBoxLayout;
	central_layout->addWidget(glyph_widget_);
    central_layout->addWidget(projection_control_widget_);
	ui_.centralWidget->setLayout(central_layout);

    variable_selection_widget_ = new VariableSelectionWidget;
    variable_selection_widget_->setMinimumWidth(700);
    var_selection_panel_ = new QDockWidget(QString("Variable Selection"), this);
	var_selection_panel_->setWidget(variable_selection_widget_);
	this->addDockWidget(Qt::RightDockWidgetArea, var_selection_panel_);
	var_selection_panel_->setVisible(true);

	tree_map_ = new TreeMap;
	tree_map_->setMinimumWidth(700);
	tree_map_panel_ = new QDockWidget(QString("Tree Map"), this);
	tree_map_panel_->setWidget(tree_map_);
	this->tabifyDockWidget(var_selection_panel_, tree_map_panel_);
	tree_map_panel_->setVisible(false);

    detailed_data_tableview_ = new QTableView;
    detailed_data_tableview_->setMinimumWidth(700);
    detailed_data_model_ = new QStandardItemModel;
    data_table_panel_ = new QDockWidget(QString("Data"), this);
    data_table_panel_->setWidget(detailed_data_tableview_);
    this->tabifyDockWidget(var_selection_panel_, data_table_panel_);
    data_table_panel_->setVisible(false);

    volume_renderer_ = new VolumeRenderer;
    volume_renderer_->setMinimumWidth(700);
    volume_render_panel_ = new QDockWidget(QString("Volume"), this);
    volume_render_panel_->setWidget(volume_renderer_);
    this->tabifyDockWidget(var_selection_panel_, volume_render_panel_);
    volume_render_panel_->setVisible(false);

    table_lens_ = new TableLens;
	table_lens_->setMinimumWidth(700);
	table_lens_panel_ = new QDockWidget(QString("Table Lens"), this);
	table_lens_panel_->setWidget(table_lens_);
	this->tabifyDockWidget(var_selection_panel_, table_lens_panel_);
	table_lens_panel_->setVisible(false);

    parallel_coordinate_ = new ParallelCoordinate;
    parallel_coordinate_->setFixedHeight(250);
	parallel_dataset_ = new ParallelDataset;
	parallel_coordinate_panel_ = new QDockWidget(QString("Parallel Coordinate"), this);
	parallel_coordinate_panel_->setWidget(parallel_coordinate_);
	this->addDockWidget(Qt::BottomDockWidgetArea, parallel_coordinate_panel_);
	parallel_coordinate_panel_->setVisible(false);

    connect(glyph_widget_, SIGNAL(ViewportUpdated()), this, SLOT(OnGlyphViewportUpdated()));
    connect(glyph_widget_, SIGNAL(SelectionChanged()), this, SLOT(OnGlyphSelectionChanged()));
    connect(glyph_widget_, SIGNAL(BrushSelectionChanged()), this, SLOT(OnBrushSelectionChanged()));

    connect(ui_.actionOpen_File, SIGNAL(triggered()), this, SLOT(OnActionOpenScatterFileTriggered()));
	connect(ui_.action_close, SIGNAL(triggered()), this, SLOT(OnActionCloseTriggered()));
	connect(ui_.actionExit, SIGNAL(triggered()), this, SLOT(OnActionExitTriggered()));
    connect(ui_projection_control_.applyProjectionButton, SIGNAL(clicked()), this, SLOT(OnApplyProjectionTriggered()));
    connect(ui_projection_control_.projectionComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(OnProjectionMethodChanged()));

    connect(variable_selection_widget_, SIGNAL(SelectionChanged()), this, SLOT(OnVariableSelectionChanged()));

    connect(ui_.actionSplit, SIGNAL(triggered()), this, SLOT(OnSplitClusterTriggered()));
    connect(ui_.actionRecursive_Splitting, SIGNAL(triggered()), this, SLOT(OnRecursiveSplittingTriggered()));
	connect(ui_.actionMerge, SIGNAL(triggered()), this, SLOT(OnMergeClusterTriggered()));

    color_mapping_group_ = new QActionGroup(this);
	color_mapping_group_->addAction(ui_.actionColor_Mapping_Off);
	color_mapping_group_->setExclusive(true);
	connect(color_mapping_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnVariableMappingTriggered()));
    ui_.mainToolBar->insertAction(ui_.actionShow_Minimum_Spanning_Tree, ui_.menuColor_Mapping->menuAction());

    clustering_mode_action_group_ = new QActionGroup(this);
	clustering_mode_action_group_->addAction(ui_.action_hierarchical_clustering);
	clustering_mode_action_group_->addAction(ui_.actionChameleon_Clustering);
	clustering_mode_action_group_->addAction(ui_.actionNCuts);
	clustering_mode_action_group_->addAction(ui_.actionMulti_Label);
    clustering_mode_action_group_->addAction(ui_.actionView_Dependent_Clustering);
	clustering_mode_action_group_->setExclusive(true);
	ui_.actionMulti_Label->setChecked(true);
	connect(clustering_mode_action_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnClusteringModeChanged()));
    connect(ui_.actionClusterOptions, SIGNAL(triggered()), this, SLOT(OnClusteringOptionTriggered()));

    glyph_view_interaction_mode_group_ = new QActionGroup(this);
	glyph_view_interaction_mode_group_->addAction(ui_.actionSingle_Selection);
	glyph_view_interaction_mode_group_->addAction(ui_.actionSequence_Selection);
	glyph_view_interaction_mode_group_->addAction(ui_.actionSelect_Minimum_Path);
	glyph_view_interaction_mode_group_->addAction(ui_.actionBrush_Path_Sequence);
	glyph_view_interaction_mode_group_->addAction(ui_.actionBrush_Cluster);
	glyph_view_interaction_mode_group_->setExclusive(true);
	connect(glyph_view_interaction_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnGlyphViewInteractionModeChanged()));

    // action for view visibility
    connect(ui_.actionShow_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowScatterPlotTriggered()));
	connect(ui_.actionShow_Glyph, SIGNAL(triggered()), this, SLOT(OnActionShowGlyphTriggered()));
	connect(ui_.actionShow_Tree_Map, SIGNAL(triggered()), this, SLOT(OnActionShowTreemapTriggered()));
	connect(ui_.actionShow_Table_Lens, SIGNAL(triggered()), this, SLOT(OnActionShowTableLensTriggerd()));
	connect(ui_.actionShow_PCP, SIGNAL(triggered()), this, SLOT(OnActionShowParallelCoordinateTriggered()));
    connect(ui_.actionShow_Density_Map, SIGNAL(triggered()), this, SLOT(OnActionShowDensityMapTriggered()));
    connect(ui_.actionShow_Data_Table, SIGNAL(triggered()), this, SLOT(OnActionShowDataTableTriggered()));
    connect(ui_.actionShow_Map, SIGNAL(triggered()), this, SLOT(OnActionShowMapTriggered()));
    connect(ui_.actionShow_Volume, SIGNAL(triggered()), this, SLOT(OnActionShowVolumeRendererTriggered()));

    connect(ui_.actionEvaluate_Quality, SIGNAL(triggered()), this, SLOT(OnActionEvaluateQualityTriggered()));

	//
	//connect(parallel_coordinate_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnPcpHighlightVarChanged(int)));
	//connect(tree_map_view_, SIGNAL(HighlightVarChanged(int)), this, SLOT(OnTreemapHighlightVarChagned(int)));
 //   connect(tree_map_view_, SIGNAL(NodeSelected(int)), this, SLOT(OnTreemapNodeSelected(int)));

	//progress_dialog_ = new QProgressDialog;
	//progress_dialog_->setWindowModality(Qt::WindowModal);
	//progress_dialog_->setLabelText("Processing");
	//progress_dialog_->setCancelButton(0);
	//progress_dialog_->setMinimumDuration(0);
	//progress_dialog_->setRange(0, 0);

	//move_focus_timer_.setSingleShot(true);
	//connect(&move_focus_timer_, SIGNAL(timeout()), this, SLOT(OnFocusTimerRunOut()));

	//transmap_tip_mode_group_ = new QActionGroup(this);
	//transmap_tip_mode_group_->addAction(ui_.actionOff);
	//transmap_tip_mode_group_->setExclusive(true);
	//connect(transmap_tip_mode_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnShowVarTrendTriggered()));

	example_action_group_ = new QActionGroup(this);
	example_action_group_->addAction(ui_.actionIris);
	example_action_group_->addAction(ui_.actionWine);
	example_action_group_->addAction(ui_.actionAuto_MPG);
	example_action_group_->addAction(ui_.actionWdbc);
	example_action_group_->addAction(ui_.actionMeteo_Case);
	example_action_group_->setExclusive(true);
	connect(example_action_group_, SIGNAL(triggered(QAction*)), this, SLOT(OnActionOpenExampleDataTriggered()));

	//// actions for tips on the cluster transition map
 
	//ui_.mainToolBar->insertAction(ui_.actionShow_Minimum_Spanning_Tree, ui_.menuShow_Sequence->menuAction());
	//connect(ui_.actionShow_Minimum_Spanning_Tree, SIGNAL(triggered()), this, SLOT(OnShowMstTriggered()));

	//// actions for manipulating cluster nodes
	//connect(ui_.actionBegin_Clustering, SIGNAL(triggered()), this, SLOT(OnBeginClusteringTriggered()));

	//// action for saving exploration path
	//connect(ui_.actionAdd_Path_Sequence, SIGNAL(triggered()), this, SLOT(OnSavePathSequenceTriggered()));

 //   connect(ui_.actionShow_2D_Result_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowResult2DSPTriggered()));
 //   connect(ui_.actionShow_Result_3D_Scatter_Plot, SIGNAL(triggered()), this, SLOT(OnActionShowResult3DSpTriggered()));
}

void ScatterPointGlyph::OnActionOpenScatterFileTriggered() {
    if (scatter_point_dataset_ != NULL) {
        QMessageBox::warning(this, tr("Warning"), tr("Close the current dataset before opening a new dataset!"));
        return;
    }

	/*QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.sc");
	if (file_path.length() == 0) return;

    PointDataReader point_reader;
    scatter_point_dataset_ = point_reader.LoadFile(file_path.toLocal8Bit().data());*/

    /*QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.gsc");
	if (file_path.length() == 0) return;

	GeoPointDataReader point_reader;
    scatter_point_dataset_ = point_reader.LoadFile(file_path.toLocal8Bit().data());*/


    /*VolumeDataReader reader;
    scatter_point_dataset_ = reader.LoadFile("");*/

    this->OnActionOpenVtkFileTriggered();

    /*PointDataReader point_reader;
    scatter_point_dataset_ = point_reader.LoadFile("./TestData/auto-mpg.sc");*/

    this->UpdateMenus();
}

void ScatterPointGlyph::OnActionOpenGscFileTriggered() {
    if (scatter_point_dataset_ != NULL) {
        QMessageBox::warning(this, tr("Warning"), tr("Close the current dataset before opening a new dataset!"));
        return;
    }

	QString file_path = QFileDialog::getOpenFileName(this, tr("Open Scatter Point Data"), ".", "*.gsc");
	if (file_path.length() == 0) return;

	GeoPointDataReader point_reader;
    scatter_point_dataset_ = point_reader.LoadFile(file_path.toLocal8Bit().data());

    this->InitExploration();
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

	float x_scale1 = -0.1, x_scale2 = 1.1;
	float y_scale1 = -0.1, y_scale2 = 1.1;
	double bounds[6];
	scatter_point_data->GetBounds(bounds);

    scatter_point_dataset_->var_names.resize(3, QString("var"));
    scatter_point_dataset_->var_weights.resize(3, 1.0 / 3);

    scatter_point_dataset_->original_point_pos.resize(2);
    for (int i = 0; i < 2; ++i)
        scatter_point_dataset_->original_point_pos[i].resize(scatter_point_data->GetPoints()->GetNumberOfPoints());
    scatter_point_dataset_->original_point_values.resize(3);
    for (int i = 0; i < 3; ++i)
        scatter_point_dataset_->original_point_values[i].resize(scatter_point_data->GetPoints()->GetNumberOfPoints());

	for (int i = 0; i < scatter_point_data->GetPoints()->GetNumberOfPoints(); ++i) {
		float x = scatter_point_data->GetPoint(i)[0];
		float y = scatter_point_data->GetPoint(i)[1];
        scatter_point_dataset_->original_point_pos[0][i] = x;
        scatter_point_dataset_->original_point_pos[1][i] = y;
        scatter_point_dataset_->original_point_values[0][i] = speed_array->GetValue(i);
        scatter_point_dataset_->original_point_values[1][i] = xparray->GetValue(i);
        scatter_point_dataset_->original_point_values[2][i] = yparray->GetValue(i);

		/*if (x > bounds[0] + (bounds[1] - bounds[0]) * x_scale1 && x < bounds[0] + (bounds[1] - bounds[0]) * x_scale2
			&& y > bounds[2] + (bounds[3] - bounds[2]) * y_scale1 && y < bounds[2] + (bounds[3] - bounds[2]) * y_scale2) {
			vector<float> pos;
			vector<float> value;
			pos.push_back(x);
			pos.push_back(y);
			value.push_back(speed_array->GetValue(i));
			value.push_back(xparray->GetValue(i));
			value.push_back(yparray->GetValue(i));

			scatter_point_dataset_->original_point_pos.push_back(pos);
			scatter_point_dataset_->original_point_values.push_back(value);
		}*/
	}

	scatter_point_dataset_->DirectConstruct();
}

void ScatterPointGlyph::OnActionOpenExampleDataTriggered() {
	//this->OnActionCloseTriggered();

	QString file_path = "./TestData/";
	if (ui_.actionIris->isChecked())
		file_path += "iris.sc";
	else if (ui_.actionAuto_MPG->isChecked())
		file_path += "auto-mpg.sc";
	else if (ui_.actionWine->isChecked())
		file_path += "wine.sc";
	else if (ui_.actionWdbc->isChecked())
		file_path += "wdbc.sc";
	else if (ui_.actionMeteo_Case->isChecked())
        file_path += "plot.gsc";
		//file_path += "agent_step100_status.gsc";

	if (!ui_.actionMeteo_Case->isChecked()) {
		PointDataReader point_reader;
        scatter_point_dataset_ = point_reader.LoadFile(file_path.toLocal8Bit().data());

	    QString mds_path = file_path.left(file_path.length() - 2);
	    mds_path = mds_path + QString("mds");
	    ifstream input(mds_path.toLocal8Bit().data());
	    if (input.good()) {
		    for (int i = 0; i < scatter_point_dataset_->original_point_pos[0].size(); ++i)
			    for (int j = 0; j < scatter_point_dataset_->original_point_pos.size(); ++j)
				    input >> scatter_point_dataset_->original_point_pos[j][i];
	    }
	    input.close();
	    scatter_point_dataset_->DirectConstruct();

	    this->UpdateMenus();

        this->InitExploration();

	//	this->AddPointData2View();

	//	this->label_size_factor_ = 3.0;
	} else {
		GeoPointDataReader point_reader;
        scatter_point_dataset_ = point_reader.LoadFile(file_path.toLocal8Bit().data());

		this->UpdateMenus();
        this->InitExploration();
	}

 //   var_selection_widget_->SetData(scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);
 //   selected_var_index_.resize(scatter_point_dataset_->var_num);
 //   for (int i = 0; i < scatter_point_dataset_->var_num; ++i) selected_var_index_[i] = i;
}

void ScatterPointGlyph::OnActionCloseTriggered() {
	this->ClearAllViews();

    if (cluster_tree_ != NULL ) delete cluster_tree_;
	cluster_tree_ = NULL;

	if (scatter_point_dataset_ != NULL) delete scatter_point_dataset_;
	scatter_point_dataset_ = NULL;
}

void ScatterPointGlyph::OnActionExitTriggered() {
	this->OnActionCloseTriggered();
	exit(0);
}

void ScatterPointGlyph::OnApplyProjectionTriggered() {
    if (scatter_point_dataset_->type() == ScatterPointDataset::POINT_DATA && scatter_point_dataset_->point_num < 1000) {
        QString projection_method = ui_projection_control_.projectionComboBox->currentText();
        if (projection_method.compare(QString("TSNE")) == 0) {
            scatter_point_dataset_->ApplyTsne();
        }
        else if (projection_method.compare(QString("MDS")) == 0) {
            scatter_point_dataset_->ApplyMds();
        }
        else if (projection_method.compare(QString("NORMAL")) == 0) {
            int axis_one = ui_projection_control_.xAxisComboBox->currentIndex();
            int axis_two = ui_projection_control_.yAxisComboBox->currentIndex();
            scatter_point_dataset_->ApplyNormal(axis_one, axis_two);
        }
    }
    this->InitExploration();
}

void ScatterPointGlyph::OnProjectionMethodChanged() {
    QString projection_method = ui_projection_control_.projectionComboBox->currentText();
    if (projection_method.compare(QString("NORMAL")) == 0) {
        ui_projection_control_.axisSelectionWidget->setVisible(true);
    } else {
        ui_projection_control_.axisSelectionWidget->setVisible(false);
    }
}

void ScatterPointGlyph::InitExploration() {
    if (scatter_point_dataset_ == NULL) {
        cout << "Point dataset is not loaded before exploration!" << endl;
        return;
    }

    if (cluster_tree_ != NULL) {
        cout << "Cluster tree is not cleared before exploration!" << endl;
        cluster_tree_->Clear();
    }

    selected_var_index_.resize(scatter_point_dataset_->var_num);
    for (int i = 0; i < selected_var_index_.size(); ++i)
        selected_var_index_[i] = i;

    switch (clustering_mode_) {
	case ScatterPointGlyph::HIER_MODE:
        cluster_tree_ = new HierarchicalTree(scatter_point_dataset_);
        cluster_tree_->AutoConstructTree(multi_label_threshold_);
		break;
	case ScatterPointGlyph::CHAMELEON_MODE:
		break;
	case ScatterPointGlyph::NCUTS_MODE:
	    cluster_tree_ = new NCutTree(scatter_point_dataset_);
        cluster_tree_->AutoConstructTree(ncuts_threshold_);
		break;
	case ScatterPointGlyph::MULTI_LABEL_MODE:
	    cluster_tree_ = new MultiLabelTree(scatter_point_dataset_);
        //cluster_tree_->SplitNodeOnce(cluster_tree_->root()->id());
        cluster_tree_->AutoConstructTree(multi_label_threshold_);
		break;
    case ScatterPointGlyph::VIEW_DEPENDENT_MODE:
        cluster_tree_ = new ViewDependentTree(scatter_point_dataset_);
        break;
    case ScatterPointGlyph::CLUSTER_PROJECTION_MODE:
        cluster_tree_ = new ClusterProjectionTree(scatter_point_dataset_);
        break;
	default:
		break;
	}

    /*if (cluster_tree_ == NULL) {
        cout << "Initializing clustering tree failed!" << endl;
        return;
    }*/

    this->InitAllViews();
    this->UpdateAllViews();
}

void ScatterPointGlyph::InitAllViews() {
    // init variable selection
    variable_selection_widget_->SetData(scatter_point_dataset_->var_names, scatter_point_dataset_->var_colors);

    // init glyph widget
    glyph_widget_->SetPointData(scatter_point_dataset_);
    glyph_widget_->SetData(glyph_dataset_);

    tree_map_->SetData(cluster_tree_);

    parallel_coordinate_->SetData(parallel_dataset_);
}

void ScatterPointGlyph::UpdateAllViews() {
	if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;
    
    this->UpdateGlyphWidget();
    this->UpdatePcp();
    this->UpdateTableLens();
    this->UpdateTreemap();
    this->UpdateDataTable();
}

void ScatterPointGlyph::ClearAllViews() {
	glyph_dataset_->Clear();
    parallel_dataset_->Clear();
}

void ScatterPointGlyph::UpdateGlyphWidget() {
    float left, right, bottom, top;
    glyph_widget_->GetViewPort(left, right, bottom, top);
    GlyphDatasetBuilder::Build(scatter_point_dataset_, selected_var_index_, 
        cluster_tree_, left, right, bottom, top, this->glyph_size_ * glyph_widget_->GetDistancePerPixel(), glyph_dataset_);
    glyph_dataset_->Modified();
}

void ScatterPointGlyph::UpdatePcp() {
    ParallelDatasetBuilder::Build(scatter_point_dataset_, selected_var_index_, selected_point_ids_, parallel_dataset_);
    parallel_dataset_->Modified();
}

void ScatterPointGlyph::UpdateTreemap() {
    tree_map_->UpdateView();
}

void ScatterPointGlyph::UpdateTableLens() {
    table_lens_->SetData(cluster_tree_, selected_cluster_ids_, selected_var_index_);
}

void ScatterPointGlyph::UpdateVolumeRender() {
    // update volume rendering
    if (this->scatter_point_dataset_->type() == ScatterPointDataset::VOLUME_DATA) {
        ScatterVolumeDataset* scatter_volume = (ScatterVolumeDataset*)scatter_point_dataset_;
        voxel_values_.resize(scatter_volume->dims()[0] * scatter_volume->dims()[1] * scatter_volume->dims()[2]);
        memset(voxel_values_.data(), 0, voxel_values_.size() * sizeof(float));

        for (int i = 0; i < selected_point_ids_.size(); ++i) {
            for (int j = 0; j < selected_point_ids_[i].size(); ++j)
                voxel_values_[selected_point_ids_[i][j]] = (i +  1) * 255;
        }
      
        float spacing[] = {0.1, 0.1, 0.2};
        volume_renderer_->SetData(scatter_volume->dims(), spacing, GL_FLOAT, voxel_values_.data());
        volume_renderer_->show();
    }
}

void ScatterPointGlyph::UpdateDataTable() {
    /*detailed_data_model_->clear();
    QStringList headers;
    for (int i = 0; i < scatter_point_dataset_->var_num; ++i)
        headers << scatter_point_dataset_->var_names[i];
    detailed_data_model_->setHorizontalHeaderLabels(headers);

    for (int i = 0; i < selected_point_ids_.size(); ++i) {
        for (int j = 0; j < selected_point_ids_[i].size(); ++j){
            int row_count = detailed_data_model_->rowCount();
            detailed_data_model_->insertRow(row_count);
            for (int k = 0; k < scatter_point_dataset_->var_num; ++k) {
                detailed_data_model_->setData(detailed_data_model_->index(row_count, k), scatter_point_dataset_->original_point_values[k][selected_point_ids_[i][j]], Qt::DisplayRole);
            }
        }
	}
    detailed_data_tableview_->setModel(detailed_data_model_);
    detailed_data_tableview_->setSortingEnabled(true);
    detailed_data_tableview_->update();*/
}

void ScatterPointGlyph::OnClusteringModeChanged() {
	ClusteringMode current_mode = clustering_mode_;

	if (cluster_tree_ != NULL) {
		QMessageBox msgbox;
		msgbox.setText("The tree structure will be destroyed if you change the clustering mode!");
		msgbox.setInformativeText("Do you want to continue?");
		msgbox.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
		int ret = msgbox.exec();
		switch (ret)
		{
		case QMessageBox::Yes:
			if (ui_.action_hierarchical_clustering->isChecked())
				current_mode = HIER_MODE;
			else if (ui_.actionChameleon_Clustering->isChecked())
				current_mode = CHAMELEON_MODE;
			else if (ui_.actionNCuts->isChecked())
				current_mode = NCUTS_MODE;
			else if (ui_.actionMulti_Label->isChecked())
				current_mode = MULTI_LABEL_MODE;
            else if (ui_.actionView_Dependent_Clustering->isChecked())
                current_mode = VIEW_DEPENDENT_MODE;
			break;
		default:
			switch (current_mode)
			{
			case ScatterPointGlyph::HIER_MODE:
				ui_.action_hierarchical_clustering->setChecked(true);
				break;
			case ScatterPointGlyph::CHAMELEON_MODE:
				ui_.action_hierarchical_clustering->setChecked(true);
				break;
			case ScatterPointGlyph::NCUTS_MODE:
				ui_.actionNCuts->setChecked(true);
				break;
			case ScatterPointGlyph::MULTI_LABEL_MODE:
				ui_.actionMulti_Label->setChecked(true);
				break;
            case ScatterPointGlyph::VIEW_DEPENDENT_MODE:
                ui_.actionView_Dependent_Clustering->setChecked(true);
                break;
			default:
				break;
			}
			return;			
		}
	} else {
		if (ui_.action_hierarchical_clustering->isChecked())
			current_mode = HIER_MODE;
		else if (ui_.actionChameleon_Clustering->isChecked())
			current_mode = CHAMELEON_MODE;
		else if (ui_.actionNCuts->isChecked())
			current_mode = NCUTS_MODE;
		else if (ui_.actionMulti_Label->isChecked())
			current_mode = MULTI_LABEL_MODE;
        else if (ui_.actionView_Dependent_Clustering->isChecked())
            current_mode = VIEW_DEPENDENT_MODE;
	}

	if (current_mode != clustering_mode_) {
		if (cluster_tree_ != NULL) {
			this->ClearAllViews();
			delete cluster_tree_;
			cluster_tree_ = NULL;
		}
		clustering_mode_ = current_mode;
	}

    //if (scatter_point_dataset_ != NULL) this->InitExploration();
}

void ScatterPointGlyph::OnClusteringOptionTriggered() {
	QDialog option_dialog;
	Ui::OptionDialog option_ui;
	option_ui.setupUi(&option_dialog);
	option_ui.ml_threshold_spinbox->setValue(multi_label_threshold_);
	option_ui.ncuts_threshold_spinbox->setValue(ncuts_threshold_);
	option_ui.hier_number_spinbox->setValue(hier_number_threshold_);
	option_ui.factor_spinbox->setValue(glyph_size_);

	if (option_dialog.exec() == QDialog::Accepted) {
		multi_label_threshold_ = option_ui.ml_threshold_spinbox->value();
		ncuts_threshold_ = option_ui.ncuts_threshold_spinbox->value();
		hier_number_threshold_ = option_ui.hier_number_spinbox->value();
		glyph_size_ = option_ui.factor_spinbox->value();
	}
}

void ScatterPointGlyph::OnGlyphViewportUpdated() {
    if (scatter_point_dataset_ == NULL || cluster_tree_ == NULL) return;
    this->UpdateGlyphWidget();
    this->UpdateTreemap();
}

void ScatterPointGlyph::OnGlyphSelectionChanged() {
    this->glyph_widget_->GetSelectedClusters(selected_cluster_ids_);

    selected_point_ids_.clear();
    for (int i = 0; i < selected_cluster_ids_.size(); ++i) {
        CNode* node = cluster_tree_->GetNode(selected_cluster_ids_[i]);
        vector<int> point_ids;
        cluster_tree_->GetNodePoints(node, point_ids);
        if (point_ids.size() != 0) selected_point_ids_.push_back(point_ids);
    }

    if (scatter_point_dataset_->type() == ScatterPointDataset::VOLUME_DATA)
        this->UpdateVolumeRender();

    if (scatter_point_dataset_->point_num < 1e5)
        this->glyph_widget_->SetSelectedPointIds(selected_point_ids_);
    this->UpdatePcp();
    this->UpdateTableLens();
    this->UpdateTreemap();
    this->UpdateDataTable();
}

void ScatterPointGlyph::OnBrushSelectionChanged() {
    if (scatter_point_dataset_ == NULL) return;

    vector<float> path;
    this->glyph_widget_->GetBrushingPath(path);

    selected_point_ids_.clear();
    selected_point_ids_.resize(1);
    float x, y;
    for (int i = 0; i < scatter_point_dataset_->point_num; ++i) {
        x = scatter_point_dataset_->original_point_pos[0][i];
        y = scatter_point_dataset_->original_point_pos[1][i];
        if (Utility::CheckInside(path, x, y)) selected_point_ids_[0].push_back(i);
    }

    if (scatter_point_dataset_->type() == ScatterPointDataset::VOLUME_DATA)
        this->UpdateVolumeRender();

    this->glyph_widget_->SetSelectedPointIds(selected_point_ids_);
    this->UpdatePcp();
    this->UpdateDataTable();
}

void ScatterPointGlyph::OnVariableSelectionChanged()
{
    vector<bool> is_selected;
    variable_selection_widget_->GetSelection(is_selected);
    selected_var_index_.clear();
    for (int i = 0; i < is_selected.size(); ++i)
        if (is_selected[i]) selected_var_index_.push_back(i);
    this->UpdateAllViews();
}

void ScatterPointGlyph::OnSplitClusterTriggered() {
    if (selected_cluster_ids_.size() != 1) {
        QMessageBox::warning(this, tr("Warning"), tr("Splitting is only supported for one selected cluster!"));
        return;
    }

    if (cluster_tree_ == NULL) {
        QMessageBox::warning(this, tr("Warning"), tr("Cluster tree is not constructed!"));
        return;
    }

    cluster_tree_->SplitNodeOnce(selected_cluster_ids_[0]);
    selected_cluster_ids_.clear();
    selected_point_ids_.clear();

    this->UpdateAllViews();
}

void ScatterPointGlyph::OnRecursiveSplittingTriggered() {
	if (selected_cluster_ids_.size() != 1) {
        QMessageBox::warning(this, tr("Warning"), tr("Splitting is only supported for one selected cluster!"));
        return;
    }

    if (cluster_tree_ == NULL) {
        QMessageBox::warning(this, tr("Warning"), tr("Cluster tree is not constructed!"));
        return;
    }

    cluster_tree_->SplitNodeRecursively(selected_cluster_ids_[0], multi_label_threshold_);
    selected_cluster_ids_.clear();
    selected_point_ids_.clear();

    this->UpdateAllViews();
}


void ScatterPointGlyph::OnMergeClusterTriggered() {
    if (selected_cluster_ids_.size() <= 1) {
        QMessageBox::warning(this, tr("Warning"), tr("Merging is only supported for multiple selected clusters!"));
        return;
    }

    if (cluster_tree_ == NULL) {
        QMessageBox::warning(this, tr("Warning"), tr("Cluster tree is not constructed!"));
        return;
    }

    cluster_tree_->MergeNodes(selected_cluster_ids_);
    selected_cluster_ids_.clear();
    selected_point_ids_.clear();

    this->UpdateAllViews();
}

void ScatterPointGlyph::OnVariableMappingTriggered() {
    int var_index = -1;
    QList< QAction* > actions = color_mapping_group_->actions();
    for (int i = 0; i < actions.size(); ++i)
	    if (actions.at(i)->isChecked()) {
		    var_index = i - 1;
		    break;
	    }
        
    if (var_index == -1) {
        this->glyph_widget_->SetColorMappingOff();
    } else {
        vector<float> values;
        values.resize(scatter_point_dataset_->point_num);
        for (int i = 0; i < scatter_point_dataset_->point_num; ++i) 
            values[i] = scatter_point_dataset_->original_point_values[var_index][i];
        this->glyph_widget_->SetColorMapping(scatter_point_dataset_->var_names[var_index], values, scatter_point_dataset_->original_value_ranges[var_index][0], scatter_point_dataset_->original_value_ranges[var_index][1]);
    }
}

void ScatterPointGlyph::UpdateMenus() {
	QList< QAction* > actions_mapping = color_mapping_group_->actions(); 
	for (int i = 1; i < actions_mapping.size(); ++i) {
		color_mapping_group_->removeAction(actions_mapping.at(i));
		ui_.menuColor_Mapping->removeAction(actions_mapping.at(i));
		delete actions_mapping.at(i);
	}

	ui_.actionColor_Mapping_Off->setChecked(true);
	for (int i = 0; i < scatter_point_dataset_->var_names.size(); ++i) {
		QAction* action = ui_.menuColor_Mapping->addAction(scatter_point_dataset_->var_names[i]);
		action->setCheckable(true);
		action->setChecked(false);
		color_mapping_group_->addAction(action);
	}

    // update axis names
    ui_projection_control_.xAxisComboBox->clear();
    ui_projection_control_.yAxisComboBox->clear();
    for (int i = 0; i < scatter_point_dataset_->var_names.size(); ++i) {
        ui_projection_control_.xAxisComboBox->addItem(scatter_point_dataset_->var_names[i]);
        ui_projection_control_.yAxisComboBox->addItem(scatter_point_dataset_->var_names[i]);
    }
    ui_projection_control_.xAxisComboBox->setCurrentIndex(0);
    ui_projection_control_.yAxisComboBox->setCurrentIndex(1);
}


void ScatterPointGlyph::OnGlyphViewInteractionModeChanged() {
	if (ui_.actionSingle_Selection->isChecked()) {

	} else if (ui_.actionSequence_Selection->isChecked()) {

	} else if (ui_.actionSelect_Minimum_Path->isChecked()) {
		
	} else if (ui_.actionBrush_Path_Sequence->isChecked()) {
		
	} else {
		
	}

    if (ui_.actionBrush_Cluster->isChecked()) {
		this->glyph_widget_->SetWidgetState(GlyphRenderingWidget::BRUSH_SELECTION);
    } else {
        this->glyph_widget_->SetWidgetState(GlyphRenderingWidget::NORMAL);
    }
}


void ScatterPointGlyph::OnActionShowScatterPlotTriggered() {
    glyph_widget_->SetPointMapVisibility(ui_.actionShow_Scatter_Plot->isChecked());
}

void ScatterPointGlyph::OnActionShowGlyphTriggered() {
	glyph_widget_->SetGlyphVisibility(ui_.actionShow_Glyph->isChecked());
}

void ScatterPointGlyph::OnActionShowDensityMapTriggered() {
    if (ui_.actionShow_Density_Map->isChecked()) {
        vector<vector<float>> point_pos;
        vector<int> cluster_index;
        point_pos.resize(2);

        vector<GlyphObject*> object_vec;
        glyph_dataset_->GetAllGlyphObjects(object_vec);
        for (int i = 0; i < object_vec.size(); ++i) {
            int cluster_id = object_vec[i]->cluster_id();
            vector<int> point_ids;
            cluster_tree_->GetNodePoints(cluster_tree_->GetNode(cluster_id), point_ids);
            for (int j = 0; j < point_ids.size(); ++j) {
                point_pos[0].push_back(scatter_point_dataset_->original_point_pos[0][point_ids[j]]);
                point_pos[1].push_back(scatter_point_dataset_->original_point_pos[1][point_ids[j]]);
                cluster_index.push_back(i);
            }
        }
        if (point_pos[0].size() > 1000) return;
        glyph_widget_->SetPointDensityInfo(point_pos, cluster_index);
    }
    glyph_widget_->SetDensityMapVisibility(ui_.actionShow_Density_Map->isChecked());
}

void ScatterPointGlyph::OnActionShowTreemapTriggered() {
	tree_map_panel_->setVisible(ui_.actionShow_Tree_Map->isChecked());
}

void ScatterPointGlyph::OnActionShowTableLensTriggerd() {
	table_lens_panel_->setVisible(ui_.actionShow_Table_Lens->isChecked());
}

void ScatterPointGlyph::OnActionShowParallelCoordinateTriggered() {
	parallel_coordinate_panel_->setVisible(ui_.actionShow_PCP->isChecked());
}

void ScatterPointGlyph::OnActionShowDataTableTriggered() {
    data_table_panel_->setVisible(ui_.actionShow_Data_Table->isChecked());
}

void ScatterPointGlyph::OnActionShowMapTriggered() {
    
}

void ScatterPointGlyph::OnActionShowVolumeRendererTriggered() {
    volume_render_panel_->setVisible(ui_.actionShow_Volume->isChecked());
}

void ScatterPointGlyph::OnActionEvaluateQualityTriggered() {
    if (cluster_tree_ == NULL) return;

    QualityMetric quality_metric;
    quality_metric.GenerateQualityMeasures(cluster_tree_);
    quality_metric.SaveMeasures("quality.txt");
}

//void ScatterPointGlyph::OnTransmapHighlightVarChanged(int var_index)
//{
//	parallel_coordinate_->SetHighlightAxis(var_index);
//	tree_map_view_->SetHighlightVarIndex(var_index);
//}
//
//void ScatterPointGlyph::OnPcpHighlightVarChanged(int var_index)
//{
//	trans_map_->HighlightVar(var_index);
//	tree_map_view_->SetHighlightVarIndex(var_index);
//}
//
//void ScatterPointGlyph::OnTreemapHighlightVarChagned(int var_index)
//{
//	trans_map_->HighlightVar(var_index);
//	parallel_coordinate_->SetHighlightAxis(var_index);
//}
//
//void ScatterPointGlyph::OnGlyphSizeChanged()
//{
//	glyph_pixel_radius_ = map_control_ui_.glyph_size_spinbox->value();
//
//	this->UpdateTransmap();
//}