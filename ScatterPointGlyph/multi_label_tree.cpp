#include "multi_label_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"
#include <queue>
#include "tour_path_generator.h"
#include "utility.h"
#include "slic_vf.h"

MultiLabelTree::MultiLabelTree(ScatterPointDataset* data)
	: TreeCommon(data) {
	processor_ = new MultiLabelProcessor;
    /*if (this->type() == TreeCommon::MULTI_LABEL_TREE)
        this->SplitNodeOnce(root_->id());*/
}

MultiLabelTree::~MultiLabelTree() {

}

void MultiLabelTree::AutoConstructTree(float std_dev_threshold) {
    if (root_ != NULL) this->SplitNodeRecursively(root_->id(), std_dev_threshold);
}

void MultiLabelTree::SplitNode(CBranch* node) {
    // if size of linked_node > SLIC_PIXEL_THRESHOLD, using slic to generate super pixels and clustering based on super pixels
    if (node->linked_nodes.size() > SLIC_PIXEL_THRESHOLD) {
        SplitOnSlic(node);
    } else {
        DirectSplit(node);
    }
}

void MultiLabelTree::SplitOnSlic(CBranch* node) {
    vector<vector<float>> pos;
	vector<vector<float>> value;
	for (int i = 0; i < node->linked_nodes.size(); ++i) {
		pos.push_back(node->linked_nodes[i]->mean_pos);
		value.push_back(node->linked_nodes[i]->mean_values);
	}

    vector<vector<vector<float>>> pixel_data;
    vector<vector<int>> node_count;
    int imw = (int)sqrt(node->linked_nodes.size()), imh = 200;
    imh = (node->top - node->bottom) / (node->right - node->left) * imw;
    while (imw * imh < 16 * SLIC_PIXEL_THRESHOLD) {
        imh *= 2;
        imw *= 2;
    }

    pixel_data.resize(imh);
    node_count.resize(imh);
    for (int i = 0; i < imh; ++i) {
        pixel_data[i].resize(imw);
        node_count[i].resize(imw, 0);
        for (int j = 0; j < imw; ++j)
            pixel_data[i][j].resize(point_dataset_->var_num, 0);
    }

    float node_width = node->right - node->left;
    float node_height = node->top - node->bottom;
    for (int i = 0; i < pos.size(); i++) {
        int xi = (int)((pos[i][0] - node->left) / node_width * (imw - 1));
        int yi = (int)((pos[i][1] - node->bottom) / node_height * (imh - 1));
        node_count[yi][xi]++;
        for (int j = 0; j < point_dataset_->var_num; j++)
            pixel_data[yi][xi][j] += value[i][j];
    }

        
    for (int i = 0; i < imh; ++i)
        for (int j = 0; j < imw; ++j)
            if (node_count[i][j] != 0) {
                for (int k = 0; k < point_dataset_->var_num; ++k)
                    pixel_data[i][j][k] /= node_count[i][j];
            }

    VectorFieldData* vfd = new VectorFieldData(pixel_data);

    double step = sqrt(imw * imh / (double)SLIC_PIXEL_THRESHOLD);
    int nc = 1.0;

    Slic slic;
    slic.GenerateSuperPixels(vfd, step, nc);
    vector<vector<int>>& spixel_index = slic.GetClusters();
    vector<vector<float>>& spixel_centers = slic.GetCenters();
    vector<int> spixel_count;
    spixel_count.resize(spixel_centers.size(), 0);
    for (int i = 0; i < imh; ++i)
        for (int j = 0; j < imw; ++j) 
            spixel_count[spixel_index[j][i]] += node_count[i][j];

    vector<vector<float>> new_pos;
    vector<vector<float>> new_value;
    vector<int> spixel_to_node_index;

    float scale_rate = node->right - node->left;
    for (int i = 0; i < spixel_count.size(); ++i)
        if (spixel_count[i] != 0) {
            new_pos.push_back(vector<float>{spixel_centers[i][0] / (imw - 1) * scale_rate, spixel_centers[i][1] / (imh - 1) * scale_rate});
            vector<float> temp_value;
            for (int j = 2; j < spixel_centers[i].size(); ++j)
                temp_value.push_back(spixel_centers[i][j]);
            new_value.push_back(temp_value);
            spixel_to_node_index.push_back(new_pos.size() - 1);
        }
        else {
            spixel_to_node_index.push_back(-1);
        }


    vector<int> clusters;
    //float child_radius = max((node->right - node->left), (node->top - node->bottom)) / 2;
    //float child_radius = ((node->right - node->left) + (node->top - node->bottom)) / 4;
    //float child_radius = node->radius / factor_;
    float child_radius = sqrt((node->right - node->left) * (node->top - node->bottom) / EXPECTED_CLUSTER_NUM);
    SplitPoints(new_pos, new_value, child_radius, clusters);

    vector<CBranch*> label_nodes;
	label_nodes.resize(clusters.size(), NULL);
    for (int i = 0; i < pos.size(); i++) {
        int xi = (int)((pos[i][0] - node->left) / node_width * (imw - 1));
        int yi = (int)((pos[i][1] - node->bottom) / node_height * (imh - 1));
        int label = clusters[spixel_to_node_index[spixel_index[xi][yi]]];

        if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
			label_nodes[label]->set_level(node->level() + 1);
			label_nodes[label]->radius = child_radius;
			label_nodes[label]->parent = node;
		}

		label_nodes[label]->linked_nodes.push_back(node->linked_nodes[i]);
		node->linked_nodes[i]->set_level(node->level() + 2);
		node->linked_nodes[i]->parent = label_nodes[label];
    }

	vector<CNode*> linked_nodes;
	for (int i = 0; i < clusters.size(); ++i)
		if (label_nodes[i] != NULL) {
			linked_nodes.push_back(label_nodes[i]);
		}
	for (int i = 0; i < linked_nodes.size(); ++i)
		ProgressNodeData(linked_nodes[i]);

	vector<int> tour_list;
	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
	node->linked_nodes.clear();
	for (int i = 0; i < tour_list.size(); ++i)
		node->linked_nodes.push_back(linked_nodes[tour_list[i]]);
}

void MultiLabelTree::DirectSplit(CBranch* node) {
    vector<vector<float>> pos;
	vector<vector<float>> value;
	for (int i = 0; i < node->linked_nodes.size(); ++i) {
		pos.push_back(node->linked_nodes[i]->mean_pos);
		value.push_back(node->linked_nodes[i]->mean_values);
	}

    vector<int> clusters;
    //float child_radius = min((node->right - node->left), (node->top - node->bottom)) / 2;
    //float child_radius = ((node->right - node->left) + (node->top - node->bottom)) / 4;
    //float child_radius = node->radius / factor_;
    float child_radius = sqrt((node->right - node->left) * (node->top - node->bottom) / EXPECTED_CLUSTER_NUM);
    SplitPoints(pos, value, child_radius, clusters);

    vector<CBranch*> label_nodes;
	label_nodes.resize(node->linked_nodes.size(), NULL);
	for (int i = 0; i < clusters.size(); ++i) {
		int label = clusters[i];
		if (label < 0) continue;
		if (label_nodes[label] == NULL) {
			label_nodes[label] = new CBranch;
			label_nodes[label]->set_level(node->level() + 1);
			label_nodes[label]->radius = child_radius;
			label_nodes[label]->parent = node;
		}
		label_nodes[label]->linked_nodes.push_back(node->linked_nodes[i]);
		node->linked_nodes[i]->set_level(node->level() + 2);
		node->linked_nodes[i]->parent = label_nodes[label];
	}

	vector<CNode*> linked_nodes;
	for (int i = 0; i < clusters.size(); ++i)
		if (label_nodes[i] != NULL) {
			linked_nodes.push_back(label_nodes[i]);
		}
	for (int i = 0; i < linked_nodes.size(); ++i)
		ProgressNodeData(linked_nodes[i]);

	vector<int> tour_list;
	TourPathGenerator::GenerateRoundPath(linked_nodes, tour_list);
	node->linked_nodes.clear();
	for (int i = 0; i < tour_list.size(); ++i)
		node->linked_nodes.push_back(linked_nodes[tour_list[i]]);
}

void MultiLabelTree::SplitPoints(vector<vector<float>>& pos, vector<vector<float>>& value, float radius, vector<int>& clusters) {
    vector<vector<bool>> connecting_status;
	float temp_min_length;
    //this->GetConnectionStatus(node->linked_nodes, connecting_status, temp_min_length);
	Utility::VtkTriangulation(pos, connecting_status, temp_min_length);
	for (int i = 0; i < pos.size() - 1; ++i)
		for (int j = i + 1; j < pos.size(); ++j)
			if (connecting_status[i][j]) {
				float dis = sqrt(pow(pos[i][0] - pos[j][0], 2) + pow(pos[i][1] - pos[j][1], 2));
				if (dis > 0.5 * max_radius_threshold_) {
					connecting_status[i][j] = false;
					connecting_status[j][i] = false;
				}
			}

    bool is_splitted = false;
    float temp_radius = radius;

	// as long as the node is not split, split again with a smaller radius for the label estimation
    while (!is_splitted) {
        processor_->SetLabelEstimationRadius(temp_radius);
	    processor_->SetData(pos, value, point_dataset_->var_weights, connecting_status);
	    //this->GenerateSegmentUncertainty(node->linked_nodes, connecting_status, edge_weights);
	    //processor_->SetEdgeWeights(edge_weights);
	    processor_->GenerateCluster();
        clusters = processor_->GetResultLabel();

        int first_label = -1, second_label = -1;
        for (int i = 0; i < clusters.size(); ++i) 
            if (clusters[i] >= 0) {
                first_label = clusters[i];
                for (int j = i + 1; j < clusters.size(); ++j) 
                    if (clusters[j] >= 0 && clusters[j] != first_label) {
                        second_label = clusters[j];
                        break;
                    }
                break;
            }
        if (first_label != -1 && second_label != -1 && first_label != second_label) {
            is_splitted = true;
        }
        else {
            temp_radius /= factor_;
        }
    }
}