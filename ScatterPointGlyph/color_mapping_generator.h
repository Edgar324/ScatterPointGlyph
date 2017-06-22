#ifndef COLOR_MAPPING_GENERATOR_H_
#define	COLOR_MAPPING_GENERATOR_H_

#include <QtGui/QColor>
#include <vector>

enum ColorMappingType{
	GRAY_MAPPING = 0,
    UNCERTAINTY_MAPPING,
	METEO_RAIN_MAPPING,
    METEO_TMP_MAPPING,
    RMS_MAPPING,
    PROBABILISTIC_MAPPING,
	METEO_PRES_MAPPING,
	SIMILARITY_MAPPING,
    VARIANCE_MAPPING
};

class ColorMappingGenerator
{
public:
	static ColorMappingGenerator* GetInstance();
	static bool DeleteInstance();

    typedef struct{
        float value_index;
        QColor color;
    } ColorIndex;
	
	QColor GetColor(ColorMappingType color_type, float current_value, float min_value = 0, float max_value = 1);
    void GetColorIndex(ColorMappingType color_type, std::vector<ColorIndex>& index, bool& is_linear);
    void GetQualitativeColors(int color_num, std::vector<QColor>& colors);

	void SetRmsMapping(float min_value, float max_value);
	void GetRmsMapping(float& min_value, float& max_value);

protected:
	ColorMappingGenerator();
	~ColorMappingGenerator();

private:
	static ColorMappingGenerator* instance_;

    std::vector<QColor> qualitative_colors_;
	std::vector<std::vector<ColorIndex>>color_index_vec_;
    std::vector<bool> is_linear_;
	void ConstructColorList();
};

#endif