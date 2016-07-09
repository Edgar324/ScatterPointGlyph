#include "tour_path_generator.h"
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


bool TourPathGenerator::GenerateRoundPath(std::vector<CNode*>& nodes, std::vector<int>& tour_list) {
	if (nodes.size() < 3) {
		for (int i = 0; i < nodes.size(); ++i) tour_list.push_back(i);
		return true;
	}
	typedef std::vector<simple_point< double>> PositionVec;
	typedef adjacency_matrix < undirectedS, no_property, property < edge_weight_t, double> > Graph;
	typedef graph_traits< Graph >::vertex_descriptor Vertex;

	typedef std::vector<Vertex > Container;
	typedef property_map< Graph, edge_weight_t >::type WeightMap;
	typedef property_map< Graph, vertex_index_t >::type VertexMap;

	int ncount = nodes.size();
	int varcount = nodes[0]->mean_values.size();

	PositionVec position_vec;
	position_vec.resize(ncount);
	for (int i = 0; i < ncount; ++i) {
		simple_point< double > temp;
		temp.x = nodes[i]->mean_pos[0];
		temp.y = nodes[i]->mean_pos[1];
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

	std::vector<int> temp_list;
	for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr){
		temp_list.push_back(*itr);
	}
	
	int max_index = 0;
	float max_dis = -1e10;
	for (int i = 0; i < temp_list.size() - 1; ++i) {
		float temp_dis = sqrt(pow(nodes[temp_list[i]]->mean_pos[0] - nodes[temp_list[i + 1]]->mean_pos[0], 2) + pow(nodes[temp_list[i]]->mean_pos[1] - nodes[temp_list[i + 1]]->mean_pos[1], 2));
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

bool TourPathGenerator::GenerateRoundPath(std::vector<std::vector<float>>& node_dis, std::vector<int>& tour_list) {
	typedef std::vector<simple_point< double>> PositionVec;
	typedef adjacency_matrix < undirectedS, no_property, property < edge_weight_t, double> > Graph;
	typedef graph_traits< Graph >::vertex_descriptor Vertex;

	typedef std::vector<Vertex > Container;
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

	std::vector<int> temp_list;
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