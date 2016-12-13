#include "multi_label_tree.h"
#include <queue>
#include "multivariate_dataset.h"
#include "multi_label_processor.h"
#include "utility.h"

MultiLabelTree::MultiLabelTree(MultivariateDataset* data)
	: TreeCommon(data) {
	processor_ = new MultiLabelProcessor;
    /*if (this->type() == TreeCommon::MULTI_LABEL_TREE)
        this->SplitNodeOnce(root_->id());*/
}

MultiLabelTree::~MultiLabelTree() {

}

void MultiLabelTree::SplitNode(CBranch* root_node) {
    const vector<CNode*>& children = root_node->children();
    vector<CNode*> children_copy(children.begin(), children.end());

    vector<vector<double>> pos;
	vector<vector<double>> value;
    pos.resize(children.size());
    value.resize(children.size());

	for (int i = 0; i < children.size(); ++i) {
        const vector<double>& mean_pos = children[i]->mean_pos();
        const vector<double>& mean_values = children[i]->mean_values();
        pos[i].assign(mean_pos.begin(), mean_pos.end());
        value[i].assign(mean_values.begin(), mean_values.end());
	}

    double left = 1e10, right = -1e10, bottom = 1e10, top = -1e10;
    for (int i = 0; i < pos.size(); i++) {
        if (pos[i][0] < left) left = pos[i][0];
        if (pos[i][0] > right) right = pos[i][0];
        if (pos[i][1] < bottom) bottom = pos[i][1];
        if (pos[i][1] > top) top = pos[i][1];
    }

    // Evaluate possible label radius
    double child_radius = sqrt((right - left) * (top - bottom) / EXPECTED_CLUSTER_NUM);
    // Step 1: Construct graph
    vector<vector<bool>> connecting;
    double min_edge_length;
    Utility::Triangulation(pos, connecting, min_edge_length);

    // Step 2: Clustering
    vector<vector<int>> clusters;

    // Step 3: Generate clusters
    vector<CNode*> new_children(clusters.size());
    for (int i = 0; i < clusters.size(); i++) {
        vector<CNode*> temp_nodes(clusters[i].size());
        for (int j = 0; j < clusters[i].size(); j++)
            temp_nodes[i] = children_copy[clusters[i][j]];
        new_children[i] = new CBranch(mv_dataset_, temp_nodes);
        new_children[i]->Update();
    }
    root_node->ClearChildren();
    root_node->AppendChildren(new_children);

    /// TODO: update children sequence
	/*vector<int> tour_list;
	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
	node->children.clear();
	for (int i = 0; i < tour_list.size(); ++i)
		node->children.push_back(linked_nodes[tour_list[i]]);*/
}

void MultiLabelTree::AutoConstructTree(float std_dev_threshold) {
    queue<CNode*> node_queue;
	node_queue.push(root_);

	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch == NULL || branch->level() >= MAX_ALLOWED_LEVEL) continue;

			float accu_std_dev = 0.0;
			for (int i = 0; i < mv_dataset_->var_num(); ++i) {
				accu_std_dev += branch->std_deviations()[i] * mv_dataset_->var_weights()[i];
			}
            if ((branch == root_ || accu_std_dev > std_dev_threshold) && branch->point_count() > 10) {
                SplitNode(branch);

                if (branch->level() > this->max_level_) this->max_level_ = branch->level();

                for (int i = 0; i < branch->children().size(); ++i) {
                    id_node_map_.insert(map<int, CNode*>::value_type(branch->children()[i]->id(), branch->children()[i]));
                }
            }

			for (int i = 0; i < branch->children().size(); ++i)
				node_queue.push(branch->children()[i]);
		}
	}
}

void MultiLabelTree::ConstructTree(float left, float right, float bottom, float top, float glyph_radius) {
    
}

void MultiLabelTree::Clear() {

}

void MultiLabelTree::BuildLeafs() {

}

//void MultiLabelTree::SplitOnSlic(CBranch* node) {
//    vector<vector<float>> pos;
//	vector<vector<float>> value;
//	for (int i = 0; i < node->children.size(); ++i) {
//		pos.push_back(node->children[i]->mean_pos);
//		value.push_back(node->children[i]->mean_values);
//	}
//
//    vector<vector<vector<float>>> pixel_data;
//    vector<vector<int>> node_count;
//    int imw = (int)sqrt(node->children.size()), imh = 200;
//    imh = (node->top - node->bottom) / (node->right - node->left) * imw;
//    while (imw * imh < 16 * SLIC_PIXEL_THRESHOLD) {
//        imh *= 2;
//        imw *= 2;
//    }
//
//    pixel_data.resize(imh);
//    node_count.resize(imh);
//    for (int i = 0; i < imh; ++i) {
//        pixel_data[i].resize(imw);
//        node_count[i].resize(imw, 0);
//        for (int j = 0; j < imw; ++j)
//            pixel_data[i][j].resize(mv_dataset_->var_num, 0);
//    }
//
//    float node_width = node->right - node->left;
//    float node_height = node->top - node->bottom;
//    for (int i = 0; i < pos.size(); i++) {
//        int xi = (int)((pos[i][0] - node->left) / node_width * (imw - 1));
//        int yi = (int)((pos[i][1] - node->bottom) / node_height * (imh - 1));
//        node_count[yi][xi]++;
//        for (int j = 0; j < mv_dataset_->var_num; j++)
//            pixel_data[yi][xi][j] += value[i][j];
//    }
//
//        
//    for (int i = 0; i < imh; ++i)
//        for (int j = 0; j < imw; ++j)
//            if (node_count[i][j] != 0) {
//                for (int k = 0; k < mv_dataset_->var_num; ++k)
//                    pixel_data[i][j][k] /= node_count[i][j];
//            }
//
//    VectorFieldData* vfd = new VectorFieldData(pixel_data);
//
//    double step = sqrt(imw * imh / (double)SLIC_PIXEL_THRESHOLD);
//    int nc = 1.0;
//
//    Slic slic;
//    slic.GenerateSuperPixels(vfd, step, nc);
//    vector<vector<int>>& spixel_index = slic.GetClusters();
//    vector<vector<float>>& spixel_centers = slic.GetCenters();
//    vector<int> spixel_count;
//    spixel_count.resize(spixel_centers.size(), 0);
//    for (int i = 0; i < imh; ++i)
//        for (int j = 0; j < imw; ++j) 
//            spixel_count[spixel_index[j][i]] += node_count[i][j];
//
//    vector<vector<float>> new_pos;
//    vector<vector<float>> new_value;
//    vector<int> spixel_to_node_index;
//
//    float scale_rate = node->right - node->left;
//    for (int i = 0; i < spixel_count.size(); ++i)
//        if (spixel_count[i] != 0) {
//            new_pos.push_back(vector<float>{spixel_centers[i][0] / (imw - 1) * scale_rate, spixel_centers[i][1] / (imh - 1) * scale_rate});
//            vector<float> temp_value;
//            for (int j = 2; j < spixel_centers[i].size(); ++j)
//                temp_value.push_back(spixel_centers[i][j]);
//            new_value.push_back(temp_value);
//            spixel_to_node_index.push_back(new_pos.size() - 1);
//        }
//        else {
//            spixel_to_node_index.push_back(-1);
//        }
//
//
//    vector<int> clusters;
//    //float child_radius = max((node->right - node->left), (node->top - node->bottom)) / 2;
//    //float child_radius = ((node->right - node->left) + (node->top - node->bottom)) / 4;
//    //float child_radius = node->radius / factor_;
//    float child_radius = sqrt((node->right - node->left) * (node->top - node->bottom) / EXPECTED_CLUSTER_NUM);
//    SplitPoints(new_pos, new_value, child_radius, clusters);
//
//    vector<CBranch*> label_nodes;
//	label_nodes.resize(clusters.size(), NULL);
//    for (int i = 0; i < pos.size(); i++) {
//        int xi = (int)((pos[i][0] - node->left) / node_width * (imw - 1));
//        int yi = (int)((pos[i][1] - node->bottom) / node_height * (imh - 1));
//        int label = clusters[spixel_to_node_index[spixel_index[xi][yi]]];
//
//        if (label_nodes[label] == NULL) {
//			label_nodes[label] = new CBranch;
//			label_nodes[label]->set_level(node->level() + 1);
//			label_nodes[label]->radius = child_radius;
//			label_nodes[label]->parent = node;
//		}
//
//		label_nodes[label]->linked_nodes.push_back(node->children[i]);
//		node->children[i]->set_level(node->level() + 2);
//		node->children[i]->parent = label_nodes[label];
//    }
//
//	vector<CNode*> linked_nodes;
//	for (int i = 0; i < clusters.size(); ++i)
//		if (label_nodes[i] != NULL) {
//			linked_nodes.push_back(label_nodes[i]);
//		}
//	for (int i = 0; i < linked_nodes.size(); ++i)
//		ProgressNodeData(linked_nodes[i]);
//
//	vector<int> tour_list;
//	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
//	node->children.clear();
//	for (int i = 0; i < tour_list.size(); ++i)
//		node->children.push_back(linked_nodes[tour_list[i]]);
//}
//
//void MultiLabelTree::DirectSplit(CBranch* node) {
//    vector<vector<float>> pos;
//	vector<vector<float>> value;
//	for (int i = 0; i < node->children.size(); ++i) {
//		pos.push_back(node->children[i]->mean_pos);
//		value.push_back(node->children[i]->mean_values);
//	}
//
//    vector<int> clusters;
//    //float child_radius = min((node->right - node->left), (node->top - node->bottom)) / 2;
//    //float child_radius = ((node->right - node->left) + (node->top - node->bottom)) / 4;
//    //float child_radius = node->radius / factor_;
//    float child_radius = sqrt((node->right - node->left) * (node->top - node->bottom) / EXPECTED_CLUSTER_NUM);
//    SplitPoints(pos, value, child_radius, clusters);
//
//    vector<CBranch*> label_nodes;
//	label_nodes.resize(node->children.size(), NULL);
//	for (int i = 0; i < clusters.size(); ++i) {
//		int label = clusters[i];
//		if (label < 0) continue;
//		if (label_nodes[label] == NULL) {
//			label_nodes[label] = new CBranch;
//			label_nodes[label]->set_level(node->level() + 1);
//			label_nodes[label]->radius = child_radius;
//			label_nodes[label]->parent = node;
//		}
//		label_nodes[label]->linked_nodes.push_back(node->children[i]);
//		node->children[i]->set_level(node->level() + 2);
//		node->children[i]->parent = label_nodes[label];
//	}
//
//	vector<CNode*> linked_nodes;
//	for (int i = 0; i < clusters.size(); ++i)
//		if (label_nodes[i] != NULL) {
//			linked_nodes.push_back(label_nodes[i]);
//		}
//	for (int i = 0; i < linked_nodes.size(); ++i)
//		ProgressNodeData(linked_nodes[i]);
//
//	vector<int> tour_list;
//	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
//	node->children.clear();
//	for (int i = 0; i < tour_list.size(); ++i)
//		node->children.push_back(linked_nodes[tour_list[i]]);
//}
//
//void MultiLabelTree::SplitPoints(vector<vector<float>>& pos, vector<vector<float>>& value, float radius, vector<int>& clusters) {
//    vector<vector<bool>> connecting_status;
//	float temp_min_length;
//    //this->GetConnectionStatus(node->linked_nodes, connecting_status, temp_min_length);
//	Utility::VtkTriangulation(pos, connecting_status, temp_min_length);
//	for (int i = 0; i < pos.size() - 1; ++i)
//		for (int j = i + 1; j < pos.size(); ++j)
//			if (connecting_status[i][j]) {
//				float dis = sqrt(pow(pos[i][0] - pos[j][0], 2) + pow(pos[i][1] - pos[j][1], 2));
//				if (dis > 0.5 * max_radius_threshold_) {
//					connecting_status[i][j] = false;
//					connecting_status[j][i] = false;
//				}
//			}
//
//    bool is_splitted = false;
//    float temp_radius = radius;
//
//	// as long as the node is not split, split again with a smaller radius for the label estimation
//    while (!is_splitted) {
//        processor_->SetLabelEstimationRadius(temp_radius);
//	    processor_->SetData(pos, value, mv_dataset_->var_weights, connecting_status);
//	    //this->GenerateSegmentUncertainty(node->linked_nodes, connecting_status, edge_weights);
//	    //processor_->SetEdgeWeights(edge_weights);
//	    processor_->GenerateCluster();
//        clusters = processor_->GetResultLabel();
//
//        int first_label = -1, second_label = -1;
//        for (int i = 0; i < clusters.size(); ++i) 
//            if (clusters[i] >= 0) {
//                first_label = clusters[i];
//                for (int j = i + 1; j < clusters.size(); ++j) 
//                    if (clusters[j] >= 0 && clusters[j] != first_label) {
//                        second_label = clusters[j];
//                        break;
//                    }
//                break;
//            }
//        if (first_label != -1 && second_label != -1 && first_label != second_label) {
//            is_splitted = true;
//        }
//        else {
//            temp_radius /= factor_;
//        }
//    }
//}