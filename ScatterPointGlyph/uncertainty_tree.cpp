#include "uncertainty_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"
#include <queue>
#include "tour_path_generator.h"

UncertaintyTree::UncertaintyTree(ScatterPointDataset* data)
	: TreeCommon(data),
	max_radius_threshold_(0.2),
	un_threshold_(0.2),
	sample_num_(300),
	is_segment_uncertainty_applied_(false){

	processor_ = new MultiLabelProcessor;

	common_parent_node_ = NULL;
}

UncertaintyTree::~UncertaintyTree() {

}

void UncertaintyTree::SetRadiusThreshold(float max_radius) {
	this->max_radius_threshold_ = max_radius;
}

void UncertaintyTree::SetSampleSize(int num) {
	sample_num_ = num;
}

void UncertaintyTree::SetUncertaintyThreshold(float un_threshold) {
	this->un_threshold_ = un_threshold;
}

void UncertaintyTree::GetClusterResult(float dis_per_pixel, std::vector< std::vector< int > >& cluster_index) {
	cluster_index.clear();

	std::vector< CNode* > level_node;
	this->Traverse(root_->level() - 1, level_node);

	cluster_index.resize(level_node.size());
	for (int i = 0; i < level_node.size(); ++i)
		this->Traverse(level_node[i], cluster_index[i]);
}

void UncertaintyTree::GetClusterResult(float dis_per_pixel, int& cluster_num, std::vector< int >& cluster_index) {
	cluster_index.resize(dataset_->original_point_pos.size());
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) cluster_index[i] = -1;

	int level = 0;
	float temp_radius = max_radius_threshold_;
	while (temp_radius / factor_ > dis_per_pixel) {
		temp_radius /= factor_;
		level++;
	}
	if (level == 0) level = 1;

	std::vector< CNode* > level_node;
	this->Traverse(level, level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[dataset_->sample_index[point_vec[j]]] = i;
	}

	std::cout << "Cluster Count: " << cluster_num << std::endl;
}

void UncertaintyTree::GetClusterResult(float radius, std::vector< CNode* >& level_nodes) {
	int level = 0;
	float temp_radius = max_radius_threshold_;
	while (temp_radius / factor_ > radius) {
		temp_radius /= factor_;
		level++;
	}
	if (level == 0) level = 1;

	this->Traverse(level, level_nodes);
}

void UncertaintyTree::GetActiveClusterResult(std::vector< CNode* >& level_nodes) {
	this->ActiveTraverse(level_nodes);
}

void UncertaintyTree::GetActiveClusterResult(int& cluster_num, std::vector< int >& cluster_index) {
	cluster_index.resize(dataset_->original_point_pos.size());
	for (int i = 0; i < dataset_->original_point_pos.size(); ++i) cluster_index[i] = -1;

	std::vector< CNode* > level_node;
	this->ActiveTraverse(level_node);

	cluster_num = level_node.size();

	for (int i = 0; i < level_node.size(); ++i) {
		std::vector< int > point_vec;
		this->Traverse(level_node[i], point_vec);
		for (int j = 0; j < point_vec.size(); ++j) cluster_index[dataset_->sample_index[point_vec[j]]] = i;
	}
}

void UncertaintyTree::SplitCluster(int cluster_index) {
	std::map< int, CNode* >::iterator node_iter = id_node_map_.find(cluster_index);
	if (node_iter == id_node_map_.end()) {
		std::cout << "Cluster not found!" << std::endl;
		return;
	}

	if (node_iter->second->type() != CNode::BRANCH) {
		std::cout << "Leaf cluster is not allowed to be split!" << std::endl;
		return;
	}

	CBranch* branch = (CBranch*)(node_iter->second);
	bool is_child_leaf = true;
	for (int i = 0; i < branch->linked_nodes.size(); ++i) 
		if (branch->linked_nodes[i]->type() != CNode::LEAF) {
			is_child_leaf = false;
			break;
		}
	if (is_child_leaf) this->GenerateCluster(branch);

	CBranch* temp_parent = branch->parent;
	if (temp_parent == NULL) return;

	int node_index = -1;
	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
		if (temp_parent->linked_nodes[j] == branch) {
			node_index = j;
			break;
		}
	if (node_index != -1) {
		for (int j = node_index; j < temp_parent->linked_nodes.size() - 1; ++j)
			temp_parent->linked_nodes[j] = temp_parent->linked_nodes[j + 1];
		int original_size = temp_parent->linked_nodes.size();
		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() - 1 + branch->linked_nodes.size());
		for (int i = original_size - 1; i > node_index; --i)
			temp_parent->linked_nodes[i + branch->linked_nodes.size() - 1] = temp_parent->linked_nodes[i];
		for (int i = 0; i < branch->linked_nodes.size(); ++i) {
			temp_parent->linked_nodes[i + node_index] = branch->linked_nodes[i];
			branch->linked_nodes[i]->parent = temp_parent;
		}

		UpdateChildLevel(temp_parent);
	}
}

void UncertaintyTree::MergeClusters(std::vector< int >& cluster_index) {
	std::vector< CNode* > cluster_nodes;
	for (int i = 0; i < cluster_index.size(); ++i) {
		std::map< int, CNode* >::iterator node_iter = id_node_map_.find(cluster_index[i]);
		if (node_iter == id_node_map_.end()) {
			std::cout << "Cluster not found!" << std::endl;
			return;
		}
		cluster_nodes.push_back(node_iter->second);
	}

	bool is_all_branch = true;
	for (int i = 0; i < cluster_nodes.size(); ++i) {
		if (cluster_nodes[i]->type() != CNode::BRANCH) {
			is_all_branch = false;
			break;
		}
	}

	if (is_all_branch) {
		// find the minimum level node
		int max_level_index = -10000;
		CNode* max_level_node = NULL;
		for (int i = 0; i < cluster_nodes.size(); ++i) 
			if (cluster_nodes[i]->level() > max_level_index) {
				max_level_index = cluster_nodes[i]->level();
				max_level_node = cluster_nodes[i];
			}
		if (max_level_node != NULL) {
			CBranch* parent_node = max_level_node->parent;
			CBranch* new_branch = new CBranch;
			new_branch->set_level(max_level_node->level());
			new_branch->radius = parent_node->radius / factor_;
			new_branch->parent = parent_node;
			parent_node->linked_nodes.push_back(new_branch);
			id_node_map_.insert(std::map< int, CNode* >::value_type(new_branch->id, new_branch));

			this->FindCommonParent(root_, cluster_index);

			for (int i = 0; i < cluster_nodes.size(); ++i) {
				RemoveChildNode(cluster_nodes[i], true);
				cluster_nodes[i]->parent = new_branch;
				new_branch->linked_nodes.push_back(cluster_nodes[i]);
			}

			ProgressNodeAndParent(new_branch);
			UpdateChildLevel(new_branch->parent);
		}
		return;
	}
}

void UncertaintyTree::RemoveChildNode(CNode* node, bool is_empty_deleted) {
	CBranch* temp_parent = node->parent;
	if (temp_parent == NULL) return;

	int node_index = -1;
	for (int j = 0; j < temp_parent->linked_nodes.size(); ++j)
		if (temp_parent->linked_nodes[j] == node) {
			node_index = j;
			break;
		}
	if (node_index != -1) {
		for (int j = node_index; j < temp_parent->linked_nodes.size() - 1; ++j)
			temp_parent->linked_nodes[j] = temp_parent->linked_nodes[j + 1];
		temp_parent->linked_nodes.resize(temp_parent->linked_nodes.size() - 1);
		
		if (temp_parent->linked_nodes.size() == 0 && is_empty_deleted)
			this->RemoveChildNode(temp_parent, true);
		else
			ProgressNodeAndParent(temp_parent);
	}
}

void UncertaintyTree::UpdateChildLevel(CBranch* node) {
	std::queue< CNode* > node_queue;
	for (int i = 0; i < node->linked_nodes.size(); ++i)
		node_queue.push(node->linked_nodes[i]);

	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();
		if (temp_node->type() == CNode::BRANCH) {
			CBranch* temp_branch = dynamic_cast<CBranch*>(temp_node);
			for (int i = 0; i < temp_branch->linked_nodes.size(); ++i)
				node_queue.push(temp_branch->linked_nodes[i]);
		}
		//temp_node->radius = temp_node->parent->radius / factor_;
		temp_node->set_level(temp_node->parent->level() + 1);
	}

	AssignColor(node, node->hstart, node->hend);
}

int UncertaintyTree::FindCommonParent(CNode* node, std::vector< int >& node_ids) {
	int value = 0;
	for (int i = 0; i < node_ids.size(); ++i)
		if (node->id == node_ids[i]) {
			value = 1;
		}

	if (node->type() == CNode::BRANCH) {
		CBranch* branch = dynamic_cast<CBranch*>(node);
		for (int i = 0; i < branch->linked_nodes.size(); ++i)
			value += FindCommonParent(branch->linked_nodes[i], node_ids);
	}
	if (value == node_ids.size()) {
		common_parent_node_ = node;
		value = 0;
	}
	return value;
}

void UncertaintyTree::ProgressNodeAndParent(CNode* node) {
	if (node == NULL || node == common_parent_node_) return;

	// update the statistics
	std::vector< int > point_index;
	this->Traverse(node, point_index);

	std::vector< float > average, variance, center_pos;
	average.resize(dataset_->var_num, 0);
	variance.resize(dataset_->var_num, 0);
	center_pos.resize(2, 0);
	for (int i = 0; i < point_index.size(); ++i) {
		center_pos[0] += dataset_->point_pos[point_index[i]][0];
		center_pos[1] += dataset_->point_pos[point_index[i]][1];

		for (int j = 0; j < dataset_->var_num; ++j)
			average[j] += dataset_->point_values[point_index[i]][j];
	}
	for (int i = 0; i < dataset_->var_num; ++i)
		average[i] /= point_index.size();
	center_pos[0] /= point_index.size();
	center_pos[1] /= point_index.size();

	for (int i = 0; i < point_index.size(); ++i)
		for (int j = 0; j < dataset_->var_num; ++j)
			variance[j] += pow(dataset_->point_values[point_index[i]][j] - average[j], 2);
	for (int i = 0; i < dataset_->var_num; ++i)
		variance[i] = sqrt(variance[i] / point_index.size());
	node->average_values = average;
	node->value_variance = variance;
	node->point_count = point_index.size();
	node->center_pos = center_pos;

	ProgressNodeAndParent(node->parent);
}

void UncertaintyTree::AddUserDefinedCluster(int origin_cluster, std::vector< int >& point_index) {
	if (origin_cluster == -1) {
		// do it on the top level

	} else {
		// create a new cluster based on the original cluster

	}
}

void UncertaintyTree::run() {
	if (root_->linked_nodes.size() != 0) {
		for (int i = 0; i < root_->linked_nodes.size(); ++i) delete root_->linked_nodes[i];
		root_->linked_nodes.clear();
	}
	
	this->ConstructDirectly();

	/*if (dataset_->point_num < sample_num_)
		this->ConstructDirectly();
		else
		this->ConstructOnRandomSample(sample_num_);*/

	id_node_map_.insert(std::map< int, CNode* >::value_type(root_->id, root_));

	root_->set_level(0);
	root_->is_expanded = true;
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		root_->linked_nodes[i]->set_level(1);
		root_->linked_nodes[i]->parent = root_;
		root_->linked_nodes[i]->radius = this->min_edge_length_;
		id_node_map_.insert(std::map< int, CNode* >::value_type(root_->linked_nodes[i]->id, root_->linked_nodes[i]));
	}
	root_->radius = max_radius_threshold_;

	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch == NULL || branch->radius / factor_ < this->min_edge_length_) continue;
			bool is_un_fit = false;
			for (int i = 0; i < dataset_->var_num; ++i)
				if (branch == root_ || branch->value_variance[i] > un_threshold_) {
					is_un_fit = true;
					break;
				}
			if (is_un_fit) GenerateCluster(branch);

			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				node_queue.push(branch->linked_nodes[i]);
		}
	}

	AssignColor(root_, 0, 1.0);
}

void UncertaintyTree::GenerateSegmentUncertainty(std::vector< CNode* >& nodes, std::vector< std::vector< bool > >& connecting_status, std::vector< std::vector< float > >& edge_weight) {
	std::vector< std::vector< int > > var_labels;
	var_labels.resize(dataset_->weights.size());

	for (int v = 0; v < dataset_->weights.size(); ++v) {
		std::vector< std::vector< float > > pos;
		std::vector< std::vector< float > > value;
		value.resize(nodes.size());
		for (int i = 0; i < nodes.size(); ++i) {
			pos.push_back(nodes[i]->center_pos);
			value[i].push_back(nodes[i]->average_values[v]);
		}

		std::vector< float > weights;
		weights.push_back(1.0);

		processor_->SetLabelEstimationRadius(max_radius_threshold_);
		processor_->SetData(pos, value, weights, connecting_status);
		processor_->GenerateCluster();
		var_labels[v] = processor_->GetResultLabel();
	}

	edge_weight.resize(nodes.size());
	for (int i = 0; i < nodes.size(); ++i) {
		edge_weight[i].resize(nodes.size());
		edge_weight[i].assign(nodes.size(), 0);
	}

	for (int i = 0; i < nodes.size() - 1; ++i)
		for (int j = i + 1; j < nodes.size(); ++j)
			if (connecting_status[i][j]) {
				float average = 0;
				for (int k = 0; k < dataset_->weights.size(); ++k)
					average += (var_labels[k][i] == var_labels[k][j] ? 0 : 1);
				average /= dataset_->weights.size();
				edge_weight[i][j] = average;
				edge_weight[j][i] = edge_weight[i][j];
			}
}

void UncertaintyTree::GenerateCluster(CBranch* node) {
	if (node == NULL || node->linked_nodes.size() == 0) {
		std::cout << "Need data initialization first." << std::endl;
		return;
	}

	std::vector< std::vector< float > > pos;
	std::vector< std::vector< float > > value;
	for (int i = 0; i < node->linked_nodes.size(); ++i) {
		pos.push_back(node->linked_nodes[i]->center_pos);
		value.push_back(node->linked_nodes[i]->average_values);
	}

	std::vector< std::vector< bool > > connecting_status;
	this->VtkTriangulation(node->linked_nodes, connecting_status);
	for (int i = 0; i < node->linked_nodes.size() - 1; ++i)
		for (int j = i + 1; j < node->linked_nodes.size(); ++j)
			if (connecting_status[i][j]) {
				float dis = sqrt(pow(node->linked_nodes[i]->center_pos[0] - node->linked_nodes[j]->center_pos[0], 2) + pow(node->linked_nodes[i]->center_pos[1] - node->linked_nodes[j]->center_pos[1], 2));
				if (dis > 0.5 * max_radius_threshold_) {
					connecting_status[i][j] = false;
					connecting_status[j][i] = false;
				}
			}

	std::vector< std::vector< float > > edge_weights;
	edge_weights.resize(node->linked_nodes.size());
	for (int i = 0; i < node->linked_nodes.size(); ++i)
		edge_weights[i].resize(node->linked_nodes.size(), 0.5);

	processor_->SetLabelEstimationRadius(node->radius / factor_);
	processor_->SetData(pos, value, dataset_->weights, connecting_status);
	if (is_segment_uncertainty_applied_)
		this->GenerateSegmentUncertainty(node->linked_nodes, connecting_status, edge_weights);
	processor_->SetEdgeWeights(edge_weights);
	processor_->GenerateCluster();
	std::vector< int > node_cluster_index = processor_->GetResultLabel();

	std::vector< CBranch* > label_nodes;
	label_nodes.resize(node->linked_nodes.size(), NULL);
	for (int i = 0; i < node_cluster_index.size(); ++i) {
		int label = node_cluster_index[i];
		if (label < 0) continue;
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
			label_nodes[label]->set_level(node->level() + 1);
			label_nodes[label]->radius = node->radius / factor_;
			label_nodes[label]->parent = node;

			id_node_map_.insert(std::map< int, CNode* >::value_type(label_nodes[label]->id, label_nodes[label]));
		}
		label_nodes[label]->linked_nodes.push_back(node->linked_nodes[i]);
		node->linked_nodes[i]->set_level(node->level() + 2);
		node->linked_nodes[i]->parent = label_nodes[label];
	}

	if (node->level() + 1 > max_level_) max_level_ = node->level() + 1;

	common_parent_node_ = node;

	node->linked_nodes.clear();

	std::vector< CNode* > linked_nodes;
	for (int i = 0; i < node_cluster_index.size(); ++i)
		if (label_nodes[i] != NULL) {
			node->linked_nodes.push_back(label_nodes[i]);
			linked_nodes.push_back(label_nodes[i]);
		}
	for (int i = 0; i < linked_nodes.size(); ++i)
		ProgressNodeAndParent(linked_nodes[i]);
	std::vector< int > tour_list;
	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
	node->linked_nodes.clear();
	for (int i = 0; i < tour_list.size(); ++i)
		node->linked_nodes.push_back(linked_nodes[tour_list[i]]);

	/*for (int i = 0; i < node_cluster_index.size(); ++i) 
		if (label_nodes[i] != NULL) {
			node->linked_nodes.push_back(label_nodes[i]);
			ProgressNodeAndParent(label_nodes[i]);
		}*/
}