#include "ncut_tree.h"
#include "multivariate_dataset.h"
#include <queue>
#include "utility.h"
#include "LibNCut/libNcut.h"

NCutTree::NCutTree(MultivariateDataset* data) 
	: TreeCommon(data) {
	libNcutInitialize();
}

NCutTree::~NCutTree()
{
	libNcutTerminate();
}

void NCutTree::AutoConstructTree(float std_dev_threshold) {
	//if (root_ != NULL) this->SplitNodeRecursively(root_->id(), std_dev_threshold);
}

void NCutTree::SplitNode(CBranch* node) 
{
	//int seg_num = 2;

	//std::vector<std::vector<bool>> connect_status;
 //   // TODO: add triangulation
 //   float min_edge_length;
	//Utility::Triangulation(node->children, connect_status, min_edge_length);

	//mwSize node_num = node->children.size();
	//mwSize element_num = node_num * node_num - node_num;
	//mxArray* weight_array = mxCreateSparse(node_num, node_num, element_num, mxREAL);
	//mwIndex *irs, *jcs;
	//double *sr;
	//sr = mxGetPr(weight_array);
	//irs = mxGetIr(weight_array);
	//jcs = mxGetJc(weight_array);

	//int temp_index = 0;
	//for (int j = 0; j < node_num; ++j) {
	//	jcs[j] = temp_index;

	//	for (int i = 0; i < node_num; ++i) {
	//		if (i == j) continue;
	//		if (!connect_status[i][j]) continue;

	//		float temp_dis = 0;

	//		float value_dis = 0;
	//		for (int k = 0; k < mv_dataset_->var_num; ++k)
	//			value_dis += abs(node->children[i]->mean_values[k] - node->children[j]->mean_values[k]) * mv_dataset_->var_weights[k];

	//		float pos_dis = sqrt(pow(node->children[i]->mean_pos[0] - node->children[j]->mean_pos[0], 2) + pow(node->children[i]->mean_pos[1] - node->children[j]->mean_pos[1], 2));

	//		temp_dis = data_dis_scale_ * value_dis + (1.0 - data_dis_scale_) * pos_dis;

	//		irs[temp_index] = i;
	//		sr[temp_index] = temp_dis;
	//		temp_index++;
	//	}
	//}
	//jcs[node_num] = temp_index;


	//mwSize seg_num_arr_length = 1;
	//mxArray* arraySegmentsNumber = mxCreateNumericArray(1, &seg_num_arr_length, mxINT32_CLASS, mxREAL);
	//*((int *)mxGetPr(arraySegmentsNumber)) = seg_num;

	//mxArray *arrayDiscrete = NULL;
	//mxArray *arrayEigenvectors = NULL;
	//mxArray *arrayEigenvalues = NULL;
	//bool b = mlfNcutW(3, &arrayDiscrete, &arrayEigenvectors, &arrayEigenvalues, weight_array, arraySegmentsNumber);
	//double* arrayPointer = mxGetPr(arrayDiscrete);
	//std::vector<CBranch* > branches;
	//branches.resize(seg_num);
	//for (int i = 0; i < seg_num; ++i) {
	//	branches[i] = new CBranch;
	//	branches[i]->set_level(node->level() + 1);
	//	branches[i]->parent = node;
	//	for (int j = 0; j < node_num; ++j)
	//		if (arrayPointer[i * node_num + j] != 0) branches[i]->children.push_back(node->children[j]);
	//	ProgressNodeData(branches[i]);
	//}
	//node->children.clear();
 //   for (int i = 0; i < branches.size(); ++i) {
 //       node->children.push_back(branches[i]);
 //       id_node_map_.insert(std::map< int, CNode*>::value_type(branches[i]->id(), branches[i]));
 //   }

	//mxDestroyArray(arraySegmentsNumber);
	//mxDestroyArray(arrayEigenvectors);
	//mxDestroyArray(arrayEigenvalues);
	//mxDestroyArray(arrayDiscrete);
}