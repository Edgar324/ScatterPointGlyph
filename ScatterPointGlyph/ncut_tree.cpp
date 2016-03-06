#include "ncut_tree.h"
#include "scatter_point_dataset.h"
#include <queue>
#include "utility.h"
#include "LibNCut/libNcut.h"

NCutTree::NCutTree(ScatterPointDataset* data) 
	: TreeCommon(data), un_threshold_(0.1)
{
	libNcutInitialize();
}

NCutTree::~NCutTree()
{
	libNcutTerminate();
}

void NCutTree::SetUncertaintyThreshold(float un_threshold)
{
	un_threshold_ = un_threshold;
}

void NCutTree::GenerateClusters()
{
	std::queue< CNode* > node_queue;
	node_queue.push(root_);
	while (node_queue.size() > 0) {
		CNode* temp_node = node_queue.front();
		node_queue.pop();

		if (temp_node->type() == CNode::BRANCH) {
			CBranch* branch = dynamic_cast<CBranch*>(temp_node);
			/*bool is_un_fit = false;
			for (int i = 0; i < dataset_->var_num; ++i)
				if (branch == root_ || branch->value_variance[i] > un_threshold_) {
					is_un_fit = true;
					break;
				}
			if (is_un_fit) SplitNode(branch);*/

			float accu_un = 0.0;
			for (int i = 0; i < dataset_->var_num; ++i) {
				accu_un += branch->variable_variances[i] * dataset_->var_weights[i];
			}
			if (branch == root_ || accu_un > un_threshold_) SplitNode(branch);

			for (int i = 0; i < branch->linked_nodes.size(); ++i)
				node_queue.push(branch->linked_nodes[i]);
		}
	}
}

void NCutTree::SplitNode(CBranch* node) 
{
	int seg_num = 2;

	std::vector< std::vector< bool > > connect_status;
	Utility::VtkTriangulation(node->linked_nodes, connect_status, this->min_edge_length_);

	mwSize node_num = node->linked_nodes.size();
	mwSize element_num = node_num * node_num - node_num;
	mxArray* weight_array = mxCreateSparse(node_num, node_num, element_num, mxREAL);
	mwIndex *irs, *jcs;
	double *sr;
	sr = mxGetPr(weight_array);
	irs = mxGetIr(weight_array);
	jcs = mxGetJc(weight_array);

	int temp_index = 0;
	for (int j = 0; j < node_num; ++j) {
		jcs[j] = temp_index;

		for (int i = 0; i < node_num; ++i) {
			if (i == j) continue;
			if (!connect_status[i][j]) continue;

			float temp_dis = 0;

			float value_dis = 0;
			for (int k = 0; k < dataset_->var_num; ++k)
				value_dis += abs(node->linked_nodes[i]->average_values[k] - node->linked_nodes[j]->average_values[k]) * dataset_->var_weights[k];

			float pos_dis = sqrt(pow(node->linked_nodes[i]->center_pos[0] - node->linked_nodes[j]->center_pos[0], 2) + pow(node->linked_nodes[i]->center_pos[1] - node->linked_nodes[j]->center_pos[1], 2));

			temp_dis = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;

			irs[temp_index] = i;
			sr[temp_index] = temp_dis;
			temp_index++;
		}
	}
	jcs[node_num] = temp_index;


	mwSize seg_num_arr_length = 1;
	mxArray* arraySegmentsNumber = mxCreateNumericArray(1, &seg_num_arr_length, mxINT32_CLASS, mxREAL);
	*((int *)mxGetPr(arraySegmentsNumber)) = seg_num;

	mxArray *arrayDiscrete = NULL;
	mxArray *arrayEigenvectors = NULL;
	mxArray *arrayEigenvalues = NULL;
	bool b = mlfNcutW(3, &arrayDiscrete, &arrayEigenvectors, &arrayEigenvalues, weight_array, arraySegmentsNumber);
	double* arrayPointer = mxGetPr(arrayDiscrete);
	std::vector< CBranch* > branches;
	branches.resize(seg_num);
	for (int i = 0; i < seg_num; ++i) {
		branches[i] = new CBranch;
		branches[i]->set_level(node->level() + 1);
		branches[i]->parent = node;
		for (int j = 0; j < node_num; ++j)
			if (arrayPointer[i * node_num + j] != 0) branches[i]->linked_nodes.push_back(node->linked_nodes[j]);
		ProgressNodeData(branches[i]);
	}
	node->linked_nodes.clear();
    for (int i = 0; i < branches.size(); ++i) {
        node->linked_nodes.push_back(branches[i]);
        id_node_map_.insert(std::map< int, CNode* >::value_type(branches[i]->id(), branches[i]));
    }

	mxDestroyArray(arraySegmentsNumber);
	mxDestroyArray(arrayEigenvectors);
	mxDestroyArray(arrayEigenvalues);
	mxDestroyArray(arrayDiscrete);
}