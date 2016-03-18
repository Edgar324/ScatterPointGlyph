#include "tour_path_generator.h"
#include "path_dataset.h"
#include "transmap_data.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <ctime>
#include "scatter_point_dataset.h"

#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/integer_traits.hpp>
#include <boost/type_traits/detail/ice_not.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/simple_point.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
using namespace boost;
using namespace std;

template<typename VertexListGraph, typename PointContainer,
	typename WeightMap, typename VertexIndexMap>
void connectAllEuclidean(VertexListGraph& g,
						const PointContainer& points,
						WeightMap wmap,            // Property maps passed by value
						VertexIndexMap vmap,       // Property maps passed by value
						int /*sz*/)
{
	using namespace boost;
	using namespace std;
	typedef typename graph_traits<VertexListGraph>::edge_descriptor Edge;
	typedef typename graph_traits<VertexListGraph>::vertex_iterator VItr;

	Edge e;
	bool inserted;

	pair<VItr, VItr> verts(vertices(g));
	for (VItr src(verts.first); src != verts.second; src++)
	{
		for (VItr dest(src); dest != verts.second; dest++)
		{
			if (dest != src)
			{
				double weight(sqrt(pow(
					static_cast<double>(points[vmap[*src]].x -
					points[vmap[*dest]].x), 2.0) +
					pow(static_cast<double>(points[vmap[*dest]].y -
					points[vmap[*src]].y), 2.0)));

				boost::tie(e, inserted) = add_edge(*src, *dest, g);

				wmap[e] = weight;
			}

		}

	}
}

TourPathGenerator::TourPathGenerator() {

}

TourPathGenerator::~TourPathGenerator() {

}

void TourPathGenerator::SetData(TransMapData* data) {
	trans_data_ = data;

	node_dis_.resize(trans_data_->level_one_nodes.size());
	for (int i = 0; i < trans_data_->level_one_nodes.size(); ++i) {
		node_dis_[i].resize(trans_data_->level_one_nodes.size());
		node_dis_[i].assign(trans_data_->level_one_nodes.size(), 0);
	}

	for (int i = 0; i + 1 < node_dis_.size(); ++i)
		for (int j = i + 1; j < node_dis_.size(); ++j) {
			/*for (int k = 0; k < trans_data_->dataset->var_weights.size(); ++k)
				node_dis_[i][j] += abs(trans_data_->level_one_nodes[i]->average_values[k] - trans_data_->level_one_nodes[j]->average_values[k]) * trans_data_->dataset->var_weights[k];*/
			node_dis_[i][j] = sqrt(pow(trans_data_->level_one_nodes[i]->center_pos[0] - trans_data_->level_one_nodes[j]->center_pos[0], 2)
				+ pow(trans_data_->level_one_nodes[i]->center_pos[1] - trans_data_->level_one_nodes[j]->center_pos[1], 2));
			node_dis_[j][i] = node_dis_[i][j];
		}
}

bool TourPathGenerator::GenerateSpanningTree() {
	typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, double > > Graph;
	typedef property_map< Graph, edge_weight_t >::type WeightMap;
	typedef graph_traits< Graph >::edge_descriptor Edge;

	int ncount = trans_data_->level_one_nodes.size();
	std::vector< std::pair< int, int > > node_edges;
	std::vector< float > weights;
	for (int i = 0; i < ncount - 1; ++i)
		for (int j = i + 1; j < ncount; ++j)
			if (trans_data_->node_connecting_status[i][j]) {
				node_edges.push_back(std::pair< int, int >(i, j));
				/*double temp_weight(sqrt(pow(
					static_cast<double>(trans_data_->level_one_nodes[i]->center_pos[0] - trans_data_->level_one_nodes[j]->center_pos[0]), 2.0) +
					pow(static_cast<double>(trans_data_->level_one_nodes[i]->center_pos[1] - trans_data_->level_one_nodes[j]->center_pos[1]), 2.0)));*/
                double temp_weight = 0;
                for (int k = 0; k < trans_data_->var_num; ++k) {
                    temp_weight += abs(trans_data_->level_one_nodes[i]->average_values[k] - trans_data_->level_one_nodes[j]->average_values[k]) * trans_data_->dataset->var_weights[k];
                }
				weights.push_back(temp_weight);
            }
            else {
                node_edges.push_back(std::pair< int, int >(i, j));
				weights.push_back(1000);
            }

	Graph g(ncount);
	WeightMap weight_map(get(edge_weight, g));
	for (int i = 0; i < node_edges.size(); ++i) {
		Edge e; bool inserted;
		tie(e, inserted) = add_edge(node_edges[i].first, node_edges[i].second, g);
		weight_map[e] = weights[i];
	}

	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
	std::vector < Edge > spanning_tree;

	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

	linked_edge_list.clear();
#ifdef DEBUG_ON
	std::cout << "Print the edges in the MST:" << std::endl;
#endif // DEBUG_ON

	
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
#ifdef DEBUG_ON
		std::cout << source(*ei, g) << " <--> " << target(*ei, g)
			<< " with weight of " << weight[*ei]
			<< std::endl;
#endif
		linked_edge_list.push_back(source(*ei, g));
		linked_edge_list.push_back(target(*ei, g));
	}

	return true;
}

bool TourPathGenerator::GenerateVarTrend(int var_index) {
	std::vector< int > node_index;
	std::vector< float > var_values;
	node_index.resize(trans_data_->level_one_nodes.size());
	var_values.resize(trans_data_->level_one_nodes.size());
	for (int i = 0; i < node_index.size(); ++i) {
		node_index[i] = i;
		var_values[i] = trans_data_->level_one_nodes[i]->average_values[var_index];
	}

	for (int i = 0; i < trans_data_->level_one_nodes.size() - 1; ++i)
		for (int j = i + 1; j < trans_data_->level_one_nodes.size(); ++j)
			if (var_values[i] > var_values[j]){
				float temp_value = var_values[j];
				var_values[j] = var_values[i];
				var_values[i] = temp_value;

				int temp_index = node_index[i];
				node_index[i] = node_index[j];
				node_index[j] = temp_index;
			}
	this->trans_edge_list.clear();
	for (int i = 0; i < node_index.size() - 1; ++i) {
		this->trans_edge_list.push_back(node_index[i]);
		this->trans_edge_list.push_back(node_index[i + 1]);
	}

	return true;
}

bool TourPathGenerator::GenerateMinimumPath(int begin, int end, std::vector< int >& tour_list) {
	tour_list.clear();

	std::vector< bool > is_node_used;
	is_node_used.resize(trans_data_->level_one_nodes.size(), false);

	std::vector< float > min_dis;
	min_dis.resize(trans_data_->level_one_nodes.size(), 99999);
	min_dis[begin] = 0;

	std::vector< int > pre_node;
	pre_node.resize(trans_data_->level_one_nodes.size(), -1);


	for (int i = 0; i < is_node_used.size(); ++i) {
		// find the minimum
		float temp_min_dis = 9999;
		float temp_min_index = -1;
		for (int j = 0; j < is_node_used.size(); ++j)
			if (!is_node_used[j] && min_dis[j] < temp_min_dis) {
				temp_min_dis = min_dis[j];
				temp_min_index = j;
			}
		if (temp_min_index == end) break;
		if (temp_min_index == -1) return false;

		is_node_used[temp_min_index] = true;
		// update the other
		for (int j = 0; j < is_node_used.size(); ++j)
			if (!is_node_used[j] && trans_data_->node_connecting_status[temp_min_index][j] && min_dis[temp_min_index] + node_dis_[temp_min_index][j] < min_dis[j]) {
				min_dis[j] = min_dis[temp_min_index] + node_dis_[temp_min_index][j];
				pre_node[j] = temp_min_index;
			}
	}

	tour_list.clear();
	tour_list.push_back(end);
	int cIndex = end;
	while (pre_node[cIndex] != -1) {
		tour_list.push_back(pre_node[cIndex]);
		cIndex = pre_node[cIndex];
	}
	for (int i = 0; i < tour_list.size() / 2; ++i) {
		int temp = tour_list[i];
		tour_list[i] = tour_list[tour_list.size() - i - 1];
		tour_list[tour_list.size() - i - 1] = temp;
	}

	return true;
}

bool TourPathGenerator::GenerateRoundPath() {
	if (trans_data_ == NULL) return false;
	typedef std::vector< simple_point< double > > PositionVec;
	typedef adjacency_matrix < undirectedS, no_property, property < edge_weight_t, double> > Graph;
	typedef graph_traits< Graph >::vertex_descriptor Vertex;

	typedef std::vector< Vertex > Container;
	typedef property_map< Graph, edge_weight_t >::type WeightMap;
	typedef property_map< Graph, vertex_index_t >::type VertexMap;

	int ncount = trans_data_->level_one_nodes.size();

	PositionVec position_vec;
	position_vec.resize(ncount);
	for (int i = 0; i < ncount; ++i) {
		simple_point< double > temp;
		temp.x = trans_data_->level_one_nodes[i]->center_pos[0];
		temp.y = trans_data_->level_one_nodes[i]->center_pos[1];
		position_vec[i] = temp;
	}

	Graph g(ncount);
	WeightMap weight_map(get(edge_weight, g));
	VertexMap v_map = get(vertex_index, g);
	connectAllEuclidean(g, position_vec, weight_map, v_map, ncount);

	Container c;

	double len(0.0);
	try {
		metric_tsp_approx(g, make_tsp_tour_len_visitor(g, back_inserter(c), len, weight_map));
	}
	catch (const bad_graph& e) {
#ifdef DEBUG_ON
		std::cerr << "bad_graph: " << e.what() << endl;
#endif
		return false;
	}

#ifdef DEBUG_ON
	std::cout << "Number of points: " << num_vertices(g) << endl;
	std::cout << "Number of edges: " << num_edges(g) << endl;
	std::cout << "Length of Tour: " << len << endl;
#endif

	tour_node_list.clear();
	for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr){
		tour_node_list.push_back(*itr);
	}

	return true;
}

bool TourPathGenerator::GenerateRoundPath(std::vector< CNode* >& nodes, std::vector< int >& tour_list) {
	if (nodes.size() < 3) {
		for (int i = 0; i < nodes.size(); ++i) tour_list.push_back(i);
		return true;
	}
	typedef std::vector< simple_point< double > > PositionVec;
	typedef adjacency_matrix < undirectedS, no_property, property < edge_weight_t, double> > Graph;
	typedef graph_traits< Graph >::vertex_descriptor Vertex;

	typedef std::vector< Vertex > Container;
	typedef property_map< Graph, edge_weight_t >::type WeightMap;
	typedef property_map< Graph, vertex_index_t >::type VertexMap;

	int ncount = nodes.size();
	int varcount = nodes[0]->average_values.size();

	PositionVec position_vec;
	position_vec.resize(ncount);
	for (int i = 0; i < ncount; ++i) {
		simple_point< double > temp;
		temp.x = nodes[i]->center_pos[0];
		temp.y = nodes[i]->center_pos[1];
		position_vec[i] = temp;
	}

	Graph g(ncount);
	WeightMap weight_map(get(edge_weight, g));
	VertexMap v_map = get(vertex_index, g);
	connectAllEuclidean(g, position_vec, weight_map, v_map, ncount);


	//typedef graph_traits<Graph>::edge_descriptor Edge;
	//typedef graph_traits<Graph>::vertex_iterator VItr;

	//Edge e;
	//bool inserted;

	//pair<VItr, VItr> verts(vertices(g));
	//int i = 0, j = 0;
	//for (VItr src(verts.first); src != verts.second; src++)
	//{
	//	j = 0;
	//	for (VItr dest(src); dest != verts.second; dest++)
	//	{
	//		if (dest != src)
	//		{
	//			double weight = 0;
	//			for (int v = 0; v < varcount; ++v) {
	//				weight += pow(nodes[i]->average_values[v] - nodes[j]->average_values[v], 2);
	//			}
	//			weight = sqrt(weight / varcount);
	//			/*double weight(sqrt(pow(
	//				static_cast<double>(points[vmap[*src]].x -
	//				points[vmap[*dest]].x), 2.0) +
	//				pow(static_cast<double>(points[vmap[*dest]].y -
	//				points[vmap[*src]].y), 2.0)));*/

	//			boost::tie(e, inserted) = add_edge(*src, *dest, g);

	//			weight_map[e] = weight;
	//		}
	//		j++;
	//	}
	//	i++;
	//}

	Container c;

	double len(0.0);
	try {
		metric_tsp_approx(g, make_tsp_tour_len_visitor(g, back_inserter(c), len, weight_map));
	}
	catch (const bad_graph& e) {
#ifdef DEBUG_ON
		std::cerr << "bad_graph: " << e.what() << endl;
#endif
		return false;
	}
	
#ifdef DEBUG_ON
	std::cout << "Number of points: " << num_vertices(g) << endl;
	std::cout << "Number of edges: " << num_edges(g) << endl;
	std::cout << "Length of Tour: " << len << endl;
#endif

	std::vector< int > temp_list;
	for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr){
		temp_list.push_back(*itr);
	}
	
	int max_index = 0;
	float max_dis = -1e10;
	for (int i = 0; i < temp_list.size() - 1; ++i) {
		float temp_dis = sqrt(pow(nodes[temp_list[i]]->center_pos[0] - nodes[temp_list[i + 1]]->center_pos[0], 2) + pow(nodes[temp_list[i]]->center_pos[1] - nodes[temp_list[i + 1]]->center_pos[1], 2));
		if (temp_dis > max_dis) {
			max_dis = temp_dis;
			max_index = i;
		}
	}
	tour_list.clear();
	for (int i = max_index + 1; i < temp_list.size() - 1; ++i)
		tour_list.push_back(temp_list[i]);
	for (int i = 0; i <= max_index; ++i)
		tour_list.push_back(temp_list[i]);

	return true;
}

bool TourPathGenerator::GenerateRoundPath(std::vector< std::vector< float > >& node_dis, std::vector< int >& tour_list) {
	typedef std::vector< simple_point< double > > PositionVec;
	typedef adjacency_matrix < undirectedS, no_property, property < edge_weight_t, double> > Graph;
	typedef graph_traits< Graph >::vertex_descriptor Vertex;

	typedef std::vector< Vertex > Container;
	typedef property_map< Graph, edge_weight_t >::type WeightMap;
	typedef property_map< Graph, vertex_index_t >::type VertexMap;

	int ncount = node_dis.size();

	PositionVec position_vec;
	position_vec.resize(ncount);
	for (int i = 0; i < ncount; ++i) {
		simple_point< double > temp;
		temp.x = 0;
		temp.y = 0;
		position_vec[i] = temp;
	}

	Graph g(ncount);
	WeightMap weight_map(get(edge_weight, g));
	VertexMap v_map = get(vertex_index, g);
	//connectAllEuclidean(g, position_vec, weight_map, v_map, ncount);

	typedef graph_traits< Graph >::vertex_iterator VItr;
	typedef graph_traits< Graph >::edge_descriptor Edge;

	Edge e;
	bool inserted;

	pair<VItr, VItr> verts(vertices(g));
	int x = 0;
	for (VItr src(verts.first); src != verts.second; src++)
	{	
		int y = 0;
		for (VItr dest(src); dest != verts.second; dest++)
		{
			if (dest != src)
			{
				boost::tie(e, inserted) = add_edge(*src, *dest, g);

				weight_map[e] = node_dis[x][y];
			}
			y++;
		}
		x++;
	}

	Container c;

	double len(0.0);
	try {
		metric_tsp_approx(g, make_tsp_tour_len_visitor(g, back_inserter(c), len, weight_map));
	}
	catch (const bad_graph& e) {
#ifdef DEBUG_ON
		std::cerr << "bad_graph: " << e.what() << endl;
#endif
		return false;
	}

#ifdef DEBUG_ON
	std::cout << "Number of points: " << num_vertices(g) << endl;
	std::cout << "Number of edges: " << num_edges(g) << endl;
	std::cout << "Length of Tour: " << len << endl;
#endif

	std::vector< int > temp_list;
	for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr){
		temp_list.push_back(*itr);
	}

	int max_index = 0;
	float max_dis = -1e10;
	for (int i = 0; i < temp_list.size() - 1; ++i) {
		float temp_dis = node_dis[temp_list[i]][temp_list[i + 1]];
		if (temp_dis > max_dis) {
			max_dis = temp_dis;
			max_index = i;
		}
	}
	tour_list.clear();
	for (int i = max_index + 1; i < temp_list.size() - 1; ++i)
		tour_list.push_back(temp_list[i]);
	for (int i = 0; i <= max_index; ++i)
		tour_list.push_back(temp_list[i]);

	return true;
}

bool TourPathGenerator::GeneratePath(int begin_index /* = -1 */, int end_index /* = -1 */) {
	if (begin_index == -1 || end_index == -1) return false;
	return true;
}

PathRecord* TourPathGenerator::GetPath() {
	if (trans_data_ == 0) return 0;

	PathRecord* record = new PathRecord;
	record->item_values.resize(tour_node_list.size());
	record->item_color.resize(tour_node_list.size());
	for (int i = 0; i < tour_node_list.size(); ++i) {
		record->item_values[i] = trans_data_->level_one_nodes[tour_node_list[i]]->average_values;
		record->item_color[i] = trans_data_->level_one_nodes[tour_node_list[i]]->color;
	}
	record->change_values.resize(tour_node_list.size() - 1);
	for (int i = 0; i < tour_node_list.size() - 1; ++i) {
		CNode* current_node = trans_data_->level_one_nodes[tour_node_list[i]];
		CNode* next_node = trans_data_->level_one_nodes[tour_node_list[i + 1]];
		record->change_values[i].resize(current_node->average_values.size());
		for (int j = 0; j < current_node->average_values.size(); ++j)
			record->change_values[i][j] = next_node->average_values[j] - current_node->average_values[j];
	}
	for (int k = 0; k < trans_data_->level_one_nodes[0]->average_values.size(); ++k) 
		record->var_names.push_back("Test");

	return record;
}