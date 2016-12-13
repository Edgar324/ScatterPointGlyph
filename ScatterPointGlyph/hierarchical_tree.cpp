#include "hierarchical_tree.h"
#include <queue>
#include <iostream>
#include <fstream>
#include "multivariate_dataset.h"
#include "utility.h"
#include "quality_metric.h"

HierarchicalTree::HierarchicalTree(MultivariateDataset* data) 
	: TreeCommon(data), expected_cluster_num_(6), type_(CENTER_DISTANCE) {
}

HierarchicalTree::~HierarchicalTree() {

}

void HierarchicalTree::SetExpectedClusterNum(int num) {
	this->expected_cluster_num_ = num;
}

void HierarchicalTree::SetDistanceType(DistanceType type) {
	this->type_ = type;
}

void HierarchicalTree::SplitNode(CBranch* node) {

}

void HierarchicalTree::AutoConstructTree(float std_dev_threshold) {
 //   std::vector<CNode*> leaf_nodes = root_->children;
 //   root_->children.clear();
 //   for (int i = 0; i < leaf_nodes.size(); ++i) {
 //       CBranch* branch = new CBranch;
 //       branch->children.push_back(leaf_nodes[i]);
 //       ProgressNodeData(branch);
 //       root_->children.push_back(branch);

 //       id_node_map_.insert(std::map< int, CNode*>::value_type(branch->id(), branch));
 //   }

 //   // TODO: add connecting status
 //   std::vector<std::vector<bool>> connecting_status;
 //   float min_edge_length;
 //   Utility::Triangulation(root_->children, connecting_status, min_edge_length);
 //   //std::vector<std::vector<bool>> connecting_status = this->node_connecting_status_;

 //   vector<int> cluster_num = vector<int>{1, 2, 4, 8, 16, 32, 64, 128};
 //   vector<vector<float>> quality_results;
 //   quality_results.resize(cluster_num.size());

	//while (root_->children.size() > expected_cluster_num_) {
 //       cout << "Cluster Number: " << root_->children.size() << endl;

 //       bool is_exist = false;
 //       int existing_index = -1;
 //       for (int i = cluster_num.size() - 1; i >= 0; --i) {
 //           if (root_->children.size() > cluster_num[i]) break;
 //           if (root_->children.size() == cluster_num[i]) {
 //               is_exist = true;
 //               existing_index = i;
 //               break;
 //           }
 //       }
 //       if (is_exist) {
 //           this->ResetLevel(root_, 0);
 //           QualityMetric metric;
 //           vector<float> temp_quality;
 //           metric.GenerateLvelMeasure(this, 1, temp_quality);
 //           quality_results[existing_index] = temp_quality;
 //       }

	//	int min_node_index[2];
	//	float min_node_dis = 1e20;
	//	min_node_index[0] = -1;
	//	min_node_index[1] = -1;

	//	for (int i = 0; i < root_->children.size() - 1; ++i)
	//		for (int j = i + 1; j < root_->children.size(); ++j) 
 //               if (connecting_status[i][j]) {
	//			    float temp_dis = 0;

	//			    float value_dis = 0;
	//			    for (int k = 0; k < mv_dataset_->var_num; ++k)
	//				    value_dis += abs(root_->children[i]->mean_values[k] - root_->children[j]->mean_values[k]) * mv_dataset_->var_weights[k];

	//			    float pos_dis = sqrt(pow(root_->children[i]->mean_pos[0] - root_->children[j]->mean_pos[0], 2) + pow(root_->children[i]->mean_pos[1] - root_->children[j]->mean_pos[1], 2));

	//			    temp_dis = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;

	//			    if (temp_dis < min_node_dis) {
	//				    min_node_index[0] = i;
	//				    min_node_index[1] = j;
	//				    min_node_dis = temp_dis;
	//			    }
	//		    }

	//	if (min_node_index[0] != -1 && min_node_index[1] != -1) {
	//		CBranch* branch = new CBranch;
	//		branch->children.push_back(root_->children[min_node_index[0]]);
	//		branch->children.push_back(root_->children[min_node_index[1]]);
	//		ProgressNodeData(branch);

 //           id_node_map_.insert(std::map< int, CNode*>::value_type(branch->id(), branch));

	//		root_->children[min_node_index[0]] = branch;
	//		for (int i = min_node_index[1]; i < root_->children.size() - 1; ++i)
	//			root_->children[i] = root_->children[i + 1];
	//		root_->children.resize(root_->children.size() - 1);

 //           // update connecting status
 //           for (int i = 0; i < connecting_status.size(); ++i) {
 //               connecting_status[i][min_node_index[0]] = connecting_status[i][min_node_index[0]] || connecting_status[i][min_node_index[1]];
 //               connecting_status[min_node_index[0]][i] = connecting_status[i][min_node_index[0]];
 //           }
 //           for (int i = 0; i < connecting_status.size(); ++i) {
 //               for (int j = min_node_index[1]; j < connecting_status[i].size() - 1; ++j)
 //                   connecting_status[i][j] = connecting_status[i][j + 1];
 //               connecting_status[i].resize(connecting_status[i].size() - 1);
 //           }
 //           for (int i = min_node_index[1]; i < connecting_status.size() - 1; ++i)
 //               connecting_status[i] = connecting_status[i + 1];
 //           connecting_status.resize(connecting_status.size() - 1);
	//	}
	//}

    /*std::ofstream output("quality.txt");

	if (output.good()) {
		output << quality_results.size() << std::endl;
		for (int i = 1; i < quality_results.size(); ++i)
			output << quality_results[i][0] << " " << quality_results[i][1] << " " << quality_results[i][2] << std::endl;
		output.close();
	}*/
}
