#ifndef HIER_CLUSTER_GLYPH_LAYER_H_
#define HIER_CLUSTER_GLYPH_LAYER_H_

#include "basic_glyph_layer.h"

class ScatterPointDataset;

class ClusterGlyphLayer : public BasicGlyphLayer
{
public:
	ClusterGlyphLayer();
	~ClusterGlyphLayer();

	void SetData(ScatterPointDataset* data);
	void SetRadiusRange(float maxR, float minR);

	void SetClusterIndex(int cluster_count, std::vector< int >& point_index);
	void Select(float x, float y);
	std::vector< bool >& GetClusterSelection() { return is_cluster_selected_;  }

	void ClearGlyph();

private:
	ScatterPointDataset* dataset_;
	float radius_range_[2];
	int current_cluster_count_;

	std::vector< float > center_x_, center_y_;
	std::vector< bool > is_cluster_selected_;

	void InitGlyphActor();
};

#endif