#include "multi_label_tree.h"
#include "scatter_point_dataset.h"
#include "multi_label_processor.h"
#include <queue>
#include "tour_path_generator.h"
#include "utility.h"
#include "slic_vf.h"

MultiLabelTree::MultiLabelTree(ScatterPointDataset* data)
	: TreeCommon(data),
	max_radius_threshold_(0.5),
	un_threshold_(0.1) {

	processor_ = new MultiLabelProcessor;

	common_parent_node_ = NULL;
}

MultiLabelTree::~MultiLabelTree() {

}

void MultiLabelTree::SetRadiusThreshold(float max_radius) {
	this->max_radius_threshold_ = max_radius;
}

void MultiLabelTree::SetUncertaintyThreshold(float un_threshold) {
	this->un_threshold_ = un_threshold;
}

int MultiLabelTree::GetRadiusLevel(float radius) {
	int level = 0;
	float temp_radius = max_radius_threshold_;
	while (temp_radius / (factor_ * 1.25) > radius) {
		temp_radius /= (factor_ * 1.25);
		level++;
	}
	if (level > max_level_) level = max_level_;

	return level;
}

void MultiLabelTree::BeginClustering() {
	id_node_map_.insert(map< int, CNode*>::value_type(root_->id(), root_));
	for (int i = 0; i < root_->linked_nodes.size(); ++i) {
		id_node_map_.insert(map< int, CNode*>::value_type(root_->linked_nodes[i]->id(), root_->linked_nodes[i]));
	}
	root_->radius = max_radius_threshold_;
}

void MultiLabelTree::GenerateClusters() {
	root_->radius = max_radius_threshold_;

	queue<CNode*> node_queue;
	node_queue.push(root_);
	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			if (branch == NULL || branch->level() > 5 || branch->radius / factor_ < this->min_edge_length_) continue;
			float accu_un = 0.0;
			for (int i = 0; i < dataset_->var_num; ++i) {
				accu_un += branch->variable_variances[i] * dataset_->var_weights[i];
			}

			/*bool is_un_fit = false;
			for (int i = 0; i < dataset_->var_num; ++i)
				if (branch == root_ || branch->value_variance[i] > un_threshold_) {
					is_un_fit = true;
					break;
				}
			if (is_un_fit) SplitNode(branch);*/
			if (branch == root_ || accu_un > un_threshold_) SplitNode(branch);

			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				node_queue.push(branch->linked_nodes[i]);
		}
	}
}

void MultiLabelTree::GenerateSegmentUncertainty(vector<CNode*>& nodes, vector<vector<bool>>& connecting_status, vector<vector<float>>& edge_weight) {
	vector<vector<int>> var_labels;
	var_labels.resize(dataset_->var_weights.size());

	for (int v = 0; v < dataset_->var_weights.size(); ++v) {
		vector<vector<float>> pos;
		vector<vector<float>> value;
		value.resize(nodes.size());
		for (int i = 0; i < nodes.size(); ++i) {
			pos.push_back(nodes[i]->center_pos);
			value[i].push_back(nodes[i]->average_values[v]);
		}

		vector<float> weights;
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
				for (int k = 0; k < dataset_->var_weights.size(); ++k)
					average += (var_labels[k][i] == var_labels[k][j] ? 0 : 1);
				average /= dataset_->var_weights.size();
				edge_weight[i][j] = average;
				edge_weight[j][i] = edge_weight[i][j];
			}
}

void MultiLabelTree::SplitNode(CBranch* node) {
	if (node == NULL || node->linked_nodes.size() <= 1) {
		cout << "Need data initialization first." << endl;
		return;
	}

    // check whether all the linked nodes are leaf nodes
    bool is_all_leaf = true;
    for (int i = 0; i < node->linked_nodes.size(); ++i) 
        if (node->linked_nodes[i]->type() != CNode::LEAF) {
            is_all_leaf = false;
            break;
        }
    if (!is_all_leaf) return;

	vector<vector<float>> pos;
	vector<vector<float>> value;
	for (int i = 0; i < node->linked_nodes.size(); ++i) {
		pos.push_back(node->linked_nodes[i]->center_pos);
		value.push_back(node->linked_nodes[i]->average_values);
	}

    // if size of linked_node > num_threshold, using slic to generate super pixels and clustering based on super pixels
    int num_threshold = 500;
    if (pos.size() > num_threshold) {
        vector<vector<vector<float>>> pixel_data;
        vector<vector<int>> node_count;
        int imw = 100, imh = 200;
        float minx = 1e10, maxx = -1e10, miny = 1e10, maxy = -1e10;
        for (int i = 0; i < pos.size(); ++i) {
            if (pos[i][0] > maxx) maxx = pos[i][0];
            if (pos[i][0] < minx) minx = pos[i][0];
            if (pos[i][1] > maxy) maxy = pos[i][1];
            if (pos[i][1] < miny) miny = pos[i][1];
        }
        imh = (maxy - miny) / (maxx - minx) * imw;
        while (imw * imh < 16 * num_threshold) {
            imh *= 2;
            imw *= 2;
        }

        pixel_data.resize(imh);
        node_count.resize(imh);
        for (int i = 0; i < imh; ++i) {
            pixel_data[i].resize(imw);
            node_count[i].resize(imw, 0);
            for (int j = 0; j < imw; ++j)
                pixel_data[i][j].resize(dataset_->var_num, 0);
        }

        for (int i = 0; i < pos.size(); i++) {
            int xi = (int)((pos[i][0] - minx) / (maxx - minx) * (imw - 1));
            int yi = (int)((pos[i][1] - miny) / (maxy - miny) * (imh - 1));
            node_count[yi][xi]++;
            for (int j = 0; j < dataset_->var_num; j++)
                pixel_data[yi][xi][j] += value[i][j];
        }

        
        for (int i = 0; i < imh; ++i)
            for (int j = 0; j < imw; ++j)
                if (node_count[i][j] != 0) {
                    for (int k = 0; k < dataset_->var_num; ++k)
                        pixel_data[i][j][k] /= node_count[i][j];
                }

        VectorFieldData* vfd = new VectorFieldData(pixel_data);

        double step = sqrt(imw * imh / (double)num_threshold);
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

        for (int i = 0; i < spixel_count.size(); ++i)
            if (spixel_count[i] != 0) {
                new_pos.push_back(vector<float>{spixel_centers[i][0] / (imw - 1), spixel_centers[i][1] / (imh - 1)});
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
        SplitPoints(new_pos, new_value, node->radius, clusters);

        vector<CBranch*> label_nodes;
	    label_nodes.resize(clusters.size(), NULL);
        for (int i = 0; i < pos.size(); i++) {
            int xi = (int)((pos[i][0] - minx) / (maxx - minx) * (imw - 1));
            int yi = (int)((pos[i][1] - miny) / (maxy - miny) * (imh - 1));
            int label = clusters[spixel_to_node_index[spixel_index[xi][yi]]];

            if (label_nodes[label] == NULL) {
			    label_nodes[label] = new CBranch;
			    label_nodes[label]->set_level(node->level() + 1);
			    label_nodes[label]->radius = node->radius / factor_;
			    label_nodes[label]->parent = node;

			    id_node_map_.insert(map< int, CNode*>::value_type(label_nodes[label]->id(), label_nodes[label]));
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
    } else {
        vector<int> clusters;
        SplitPoints(pos, value, node->radius, clusters);

        vector<CBranch*> label_nodes;
	    label_nodes.resize(node->linked_nodes.size(), NULL);
	    for (int i = 0; i < clusters.size(); ++i) {
		    int label = clusters[i];
		    if (label < 0) continue;
		    if (label_nodes[label] == NULL) {
			    label_nodes[label] = new CBranch;
			    label_nodes[label]->set_level(node->level() + 1);
			    label_nodes[label]->radius = node->radius / factor_;
			    label_nodes[label]->parent = node;

			    id_node_map_.insert(map< int, CNode*>::value_type(label_nodes[label]->id(), label_nodes[label]));
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
    float temp_radius = radius / factor_;

	// as long as the node is not split, split again with a smaller radius for the label estimation
    while (!is_splitted) {
        processor_->SetLabelEstimationRadius(temp_radius);
	    processor_->SetData(pos, value, dataset_->var_weights, connecting_status);
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