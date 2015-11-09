#ifndef GESTALT_PROCESSOR_H_
#define GESTALT_PROCESSOR_H_

#include <vector>

class GestaltProcessor
{
public:
	GestaltProcessor();
	~GestaltProcessor();

	void SetData(std::vector< std::vector< float > >& pos, std::vector< std::vector< float > >& values);

private:
	void ExtractGestalts();
	void LabelOptimization();
};

#endif