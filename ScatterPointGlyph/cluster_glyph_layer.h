#ifndef HIER_CLUSTER_GLYPH_LAYER_H_
#define HIER_CLUSTER_GLYPH_LAYER_H_

#include "basic_glyph_layer.h"

class ClusterGlyphLayer : public BasicGlyphLayer
{
public:
	ClusterGlyphLayer();
	~ClusterGlyphLayer();

	void SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& value);
	void SetHighlighCluster(int cluster_one, int cluster_two);
	void AddGestaltGlyph(std::vector< int >& point_index);
	void ClearGlyph();

private:
	std::vector< std::vector< float > > glyph_pos_;
	std::vector< std::vector< float > > glyph_values_;

	std::vector< float > line_paras_;

	int highlight_cluster_[2];

	void UpdateGlyphActor();
};

#endif