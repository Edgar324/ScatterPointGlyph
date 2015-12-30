#include "tour_path_generator.h"
#include "path_dataset.h"
#include "transmap_data.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <ctime>

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
			if (/*trans_data_->node_connecting_status[i][j]*/true) {
				node_edges.push_back(std::pair< int, int >(i, j));
				double temp_weight(sqrt(pow(
					static_cast<double>(trans_data_->level_one_nodes[i]->center_pos[0] - trans_data_->level_one_nodes[j]->center_pos[0]), 2.0) +
					pow(static_cast<double>(trans_data_->level_one_nodes[i]->center_pos[1] - trans_data_->level_one_nodes[j]->center_pos[1]), 2.0)));
				weights.push_back(temp_weight);
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

	edge_list.clear();
	std::cout << "Print the edges in the MST:" << std::endl;
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
		std::cout << source(*ei, g) << " <--> " << target(*ei, g)
			<< " with weight of " << weight[*ei]
			<< std::endl;
		edge_list.push_back(source(*ei, g));
		edge_list.push_back(target(*ei, g));
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
		std::cerr << "bad_graph: " << e.what() << endl;
		return -1;
	}

	std::cout << "Number of points: " << num_vertices(g) << endl;
	std::cout << "Number of edges: " << num_edges(g) << endl;
	std::cout << "Length of Tour: " << len << endl;

	tour_node_list.clear();
	for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr){
		tour_node_list.push_back(*itr);
	}

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
		record->item_color[i] = QColor(trans_data_->level_one_colors[3 * tour_node_list[i]], trans_data_->level_one_colors[3 * tour_node_list[i] + 1], trans_data_->level_one_colors[3 * tour_node_list[i] + 2]);
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