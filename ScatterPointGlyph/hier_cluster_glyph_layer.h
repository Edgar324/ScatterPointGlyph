#ifndef HIER_CLUSTER_GLYPH_LAYER_H_
#define HIER_CLUSTER_GLYPH_LAYER_H_

#include "basic_glyph_layer.h"

class HierClusterGlyphLayer : public BasicGlyphLayer
{
public:
	HierClusterGlyphLayer();
	~HierClusterGlyphLayer();

	void SetData(std::vector< float >& pos, std::vector< std::vector< float > >& value);
	void SetHighlighCluster(int cluster_one, int cluster_two);

private:
	std::vector< float > glyph_pos_;
	std::vector< std::vector< float > > glyph_values_;

	int highlight_cluster_[2];

	void UpdateGlyphActor();
};

#endif