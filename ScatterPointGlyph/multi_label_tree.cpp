#include "multi_label_tree.h"
#include <queue>
#include <iostream>
using namespace std;
#include "multivariate_dataset.h"
#include "multi_label_processor.h"
#include "utility.h"
#include "slic_vf.h"

MultiLabelTree::MultiLabelTree(MultivariateDataset* data)
	: TreeCommon(data) {
	processor_ = new MultiLabelProcessor;
}

MultiLabelTree::~MultiLabelTree() {

}

void MultiLabelTree::GetNodes(int level, vector<CNode*>& level_nodes) {
    level_nodes.clear();

    queue<CNode*> node_queue;
    node_queue.push(root_);

    while (!node_queue.empty()) {
        CNode* temp_node = node_queue.front();
        node_queue.pop();

        if (temp_node->level() == level && temp_node->type() != CNode::LEAF) {
            level_nodes.push_back(temp_node);
        } else if (temp_node->type() == CNode::BRANCH) {
            CBranch* branch_node = (CBranch*)temp_node;
            bool is_children_all_leaf = true;
            for (int i = 0; i < branch_node->children().size(); ++i)
                if (branch_node->children()[i]->type() == CNode::BRANCH) {
                    is_children_all_leaf = false;
                    break;
                }
            if (is_children_all_leaf) {
                level_nodes.push_back(temp_node);
            } else {
                for (int i = 0; i < branch_node->children().size(); ++i)
                node_queue.push(branch_node->children()[i]);
            }
        }
    }
}

void MultiLabelTree::GetNodes(float left, float right, float bottom, float top, vector<CNode*>& nodes) {
    cout << "View dependent GetNodes is not supported in multilabel tree" << endl;
    return;
}

void MultiLabelTree::SplitNode(CBranch* root_node) {
    const vector<CNode*>& children = root_node->children();

    vector<vector<int>> clusters;
    if (is_pixel_used_)
        SplitOnPixels(children, clusters);
    else
        SplitOnRecords(children, clusters);

    vector<CNode*> new_children(clusters.size());
    for (int i = 0; i < clusters.size(); i++) {
        vector<CNode*> temp_nodes(clusters[i].size());
        for (int j = 0; j < clusters[i].size(); j++)
            temp_nodes[j] = children[clusters[i][j]];
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
    if (mv_dataset_->record_num() < SCALE_NODE_SIZE) {
        BuildOnRecords();
        is_pixel_used_ = false;
    } else {
        BuildOnPixels();
        is_pixel_used_ = true;
    }
}

void MultiLabelTree::BuildOnPixels() {
    double left, right, bottom, top;
    mv_dataset_->GetPosRanges(left, right, bottom, top);

    int imw = LEAF_RESO_SIZE;
    int imh = (top - bottom) / (right - left) * imw;

    vector<CLeaf*> grid_nodes;
    grid_nodes.resize(imw * imh, NULL);

    for (int i = 0; i < mv_dataset_->record_num(); ++i) {
        double* record_pos = mv_dataset_->GetRecordPos(i);

        int x = (int)(record_pos[0] * (imw - 1));
        int y = (int)(record_pos[1] * (imh - 1));
        int grid_index = y * imw + x;
        
        if (grid_nodes[grid_index] == NULL) {
            grid_nodes[grid_index] = new CLeaf(mv_dataset_, i);

            id_node_map_.insert(map<int, CNode*>::value_type(grid_nodes[grid_index]->id(), grid_nodes[grid_index]));
        } else {
            grid_nodes[grid_index]->AppendRecord(i);
        }
    }
    
    leaf_nodes_.clear();
    for (int i = 0; i < grid_nodes.size(); i++)
        if (grid_nodes[i] != NULL) {
            grid_nodes[i]->Update();
            leaf_nodes_.push_back(grid_nodes[i]);
        }
}

void MultiLabelTree::BuildOnRecords() {
    leaf_nodes_.resize(mv_dataset_->record_num());
    for (int i = 0; i < mv_dataset_->record_num(); ++i) {
        double* record_pos = mv_dataset_->GetRecordPos(i);

        leaf_nodes_[i] = new CLeaf(mv_dataset_, i);
        leaf_nodes_[i]->Update();

        id_node_map_.insert(map<int, CNode*>::value_type(leaf_nodes_[i]->id(), leaf_nodes_[i]));
    }
}

void MultiLabelTree::SplitOnPixels(const vector<CNode*>& children, vector<vector<int>>& clusters) {
    double left = 1e10, right = -1e10, bottom = 1e10, top = -1e10;
    for (int i = 0; i < children.size(); i++) {
        const vector<double>& mean_pos = children[i]->mean_pos();
        if (mean_pos[0] < left) left = mean_pos[0];
        if (mean_pos[0] > right) right = mean_pos[0];
        if (mean_pos[1] < bottom) bottom = mean_pos[1];
        if (mean_pos[1] > top) top = mean_pos[1];
    }

    int imw = SLIC_RESO_SIZE;
    int imh = (top - bottom) / (right - left) * imw;
    while (imw * imh < 16 * SLIC_PIXEL_THRESHOLD) {
        imh *= 2;
        imw *= 2;
    }

    vector<vector<vector<double>>> pixel_data;
    vector<vector<int>> pixel_node_count;

    pixel_data.resize(imh);
    pixel_node_count.resize(imh);
    for (int i = 0; i < imh; ++i) {
        pixel_data[i].resize(imw);
        pixel_node_count[i].resize(imw, 0);
        for (int j = 0; j < imw; ++j)
            pixel_data[i][j].resize(mv_dataset_->var_num(), 0);
    }

    float node_width = right - left;
    float node_height = top - bottom;
    for (int i = 0; i < children.size(); i++) {
        const vector<double>& mean_pos = children[i]->mean_pos();
        const vector<double>& mean_values = children[i]->mean_values();
        int xi = (int)((mean_pos[0] - left) / node_width * (imw - 1));
        int yi = (int)((mean_pos[1] - bottom) / node_height * (imh - 1));
        for (int j = 0; j < mv_dataset_->var_num(); j++)
            pixel_data[yi][xi][j] += mean_values[j];

        pixel_node_count[yi][xi]++;
    }

    for (int i = 0; i < imh; ++i)
        for (int j = 0; j < imw; ++j)
            if (pixel_node_count[i][j] != 0) {
                for (int k = 0; k < mv_dataset_->var_num(); ++k)
                    pixel_data[i][j][k] /= pixel_node_count[i][j];
            }

    VectorFieldData* vfd = new VectorFieldData(pixel_data);

    double step = sqrt(imw * imh / (double)SLIC_PIXEL_THRESHOLD);
    int nc = 1.0;

    Slic slic;
    slic.GenerateSuperPixels(vfd, step, nc);
    vector<vector<int>>& spixel_index = slic.GetClusters();
    vector<vector<double>>& spixel_centers = slic.GetCenters();
    vector<int> spixel_count;
    spixel_count.resize(spixel_centers.size(), 0);
    for (int i = 0; i < imh; ++i)
        for (int j = 0; j < imw; ++j) 
            spixel_count[spixel_index[j][i]] += pixel_node_count[i][j];

    vector<vector<double>> node_pos;
    vector<vector<double>> node_value;
    vector<int> spixel_to_node_index;

    float scale_rate = right - left;
    for (int i = 0; i < spixel_count.size(); ++i)
        if (spixel_count[i] != 0) {
            node_pos.push_back(vector<double>{spixel_centers[i][0] / (imw - 1) * scale_rate, spixel_centers[i][1] / (imh - 1) * scale_rate});
            node_value.push_back(vector<double>(spixel_centers[i].begin() + 2, spixel_centers[i].end()));
            spixel_to_node_index.push_back(node_pos.size() - 1);
        }
        else {
            spixel_to_node_index.push_back(-1);
        }


    vector<vector<int>> slic_clusters;
    double child_radius = sqrt((right - left) * (top - bottom) / EXPECTED_CLUSTER_NUM);
    SplitPoints(node_pos, node_value, child_radius, slic_clusters);

    vector<int> node_to_cluster_index(node_pos.size());
    for (int i = 0; i < slic_clusters.size(); i++)
        for (int j = 0; j < slic_clusters[i].size(); j++)
            node_to_cluster_index[slic_clusters[i][j]] = i;

    clusters.clear();
    clusters.resize(slic_clusters.size());

    for (int i = 0; i < children.size(); i++) {
        const vector<double>& mean_pos = children[i]->mean_pos();
        const vector<double>& mean_values = children[i]->mean_values();
        int xi = (int)((mean_pos[0] - left) / node_width * (imw - 1));
        int yi = (int)((mean_pos[1] - bottom) / node_height * (imh - 1));
        int sp = spixel_index[xi][yi];
        int no = spixel_to_node_index[sp];
        int c = node_to_cluster_index[no];
        clusters[c].push_back(i);
    }
}

void MultiLabelTree::SplitOnRecords(const vector<CNode*>& children, vector<vector<int>>& clusters) {
    double left = 1e10, right = -1e10, bottom = 1e10, top = -1e10;
    for (int i = 0; i < children.size(); i++) {
        const vector<double>& mean_pos = children[i]->mean_pos();
        if (mean_pos[0] < left) left = mean_pos[0];
        if (mean_pos[0] > right) right = mean_pos[0];
        if (mean_pos[1] < bottom) bottom = mean_pos[1];
        if (mean_pos[1] > top) top = mean_pos[1];
    }

    vector<vector<double>> node_pos;
	vector<vector<double>> node_value;
    node_pos.resize(children.size());
    node_value.resize(children.size());
	for (int i = 0; i < children.size(); i++) {
        node_pos[i] = children[i]->mean_pos();
        node_value[i] = children[i]->mean_values();
    }

    double child_radius = sqrt((right - left) * (top - bottom) / EXPECTED_CLUSTER_NUM);
    SplitPoints(node_pos, node_value, child_radius, clusters);
}

void MultiLabelTree::SplitPoints(vector<vector<double>>& pos, vector<vector<double>>& value, double radius, vector<vector<int>>& clusters) {
    vector<vector<bool>> connecting;
	Utility::Triangulation(pos, connecting);
    processor_->GenerateCluster(pos, value, mv_dataset_->var_weights(), connecting, radius, clusters);
}