#include "hierarchical_tree.h"
#include <queue>
#include "scatter_point_dataset.h"

HierarchicalTree::HierarchicalTree(ScatterPointDataset* data) 
	: TreeCommon(data), expected_cluster_num_(7), type_(CENTER_DISTANCE) {
}

HierarchicalTree::~HierarchicalTree() {

}

void HierarchicalTree::SetExpectedClusterNum(int num) {
	this->expected_cluster_num_ = num;
}

void HierarchicalTree::SetDistanceType(DistanceType type) {
	this->type_ = type;
}

void HierarchicalTree::GenerateClusters() {
    std::vector< CNode* > leaf_nodes = root_->linked_nodes;
    root_->linked_nodes.clear();
    for (int i = 0; i < leaf_nodes.size(); ++i) {
        CBranch* branch = new CBranch;
        branch->linked_nodes.push_back(leaf_nodes[i]);
        ProgressNodeData(branch);
        root_->linked_nodes.push_back(branch);
    }

    std::vector< std::vector< bool > > connecting_status = this->node_connecting_status_;

	while (root_->linked_nodes.size() > expected_cluster_num_) {
		int min_node_index[2];
		float min_node_dis = 1e20;
		min_node_index[0] = -1;
		min_node_index[1] = -1;

		for (int i = 0; i < root_->linked_nodes.size() - 1; ++i)
			for (int j = i + 1; j < root_->linked_nodes.size(); ++j) 
                if (connecting_status[i][j]) {
				    float temp_dis = 0;

				    float value_dis = 0;
				    for (int k = 0; k < dataset_->var_num; ++k)
					    value_dis += abs(root_->linked_nodes[i]->average_values[k] - root_->linked_nodes[j]->average_values[k]) * dataset_->var_weights[k];

				    float pos_dis = sqrt(pow(root_->linked_nodes[i]->center_pos[0] - root_->linked_nodes[j]->center_pos[0], 2) + pow(root_->linked_nodes[i]->center_pos[1] - root_->linked_nodes[j]->center_pos[1], 2));

				    temp_dis = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;

				    if (temp_dis < min_node_dis) {
					    min_node_index[0] = i;
					    min_node_index[1] = j;
					    min_node_dis = temp_dis;
				    }
			    }

		if (min_node_index[0] != -1 && min_node_index[1] != -1) {
			CBranch* branch = new CBranch;
			branch->linked_nodes.push_back(root_->linked_nodes[min_node_index[0]]);
			branch->linked_nodes.push_back(root_->linked_nodes[min_node_index[1]]);
			ProgressNodeData(branch);

			root_->linked_nodes[min_node_index[0]] = branch;
			for (int i = min_node_index[1]; i < root_->linked_nodes.size() - 1; ++i)
				root_->linked_nodes[i] = root_->linked_nodes[i + 1];
			root_->linked_nodes.resize(root_->linked_nodes.size() - 1);

            // update connecting status
            for (int i = 0; i < connecting_status.size(); ++i) {
                connecting_status[i][min_node_index[0]] = connecting_status[i][min_node_index[0]] || connecting_status[i][min_node_index[1]];
                connecting_status[min_node_index[0]][i] = connecting_status[i][min_node_index[0]];
            }
            for (int i = 0; i < connecting_status.size(); ++i) {
                for (int j = min_node_index[1]; j < connecting_status[i].size() - 1; ++j)
                    connecting_status[i][j] = connecting_status[i][j + 1];
                connecting_status[i].resize(connecting_status[i].size() - 1);
            }
            for (int i = min_node_index[1]; i < connecting_status.size() - 1; ++i)
                connecting_status[i] = connecting_status[i + 1];
            connecting_status.resize(connecting_status.size() - 1);
		}
	}
}

void HierarchicalTree::SplitNode(CBranch* node) {

}