#include "color_mapping_generator.h"

ColorMappingGenerator* ColorMappingGenerator::instance_ = NULL;

ColorMappingGenerator* ColorMappingGenerator::GetInstance(){
    if ( instance_ == NULL ){
        instance_ = new ColorMappingGenerator;
    }

    return instance_;
}

bool ColorMappingGenerator::DeleteInstance(){
    if ( instance_ != NULL ){
        delete instance_;
        return true;
    } else 
        return false;
}

ColorMappingGenerator::ColorMappingGenerator(){
    ConstructColorList();
}

ColorMappingGenerator::~ColorMappingGenerator(){

}

QColor ColorMappingGenerator::GetColor(ColorMappingType color_type, float current_value, float min_value, float max_value){
    float value = current_value;
    if ( value < min_value ) value = min_value;
    if ( value > max_value ) value = max_value;

	if ( value < color_index_vec_[color_type][0].value_index ) value = color_index_vec_[color_type][0].value_index;
	if ( value > color_index_vec_[color_type][color_index_vec_[color_type].size() - 1].value_index ) value = color_index_vec_[color_type][color_index_vec_[color_type].size() - 1].value_index;

    int current_index = 0;
    while ( current_index < color_index_vec_[color_type].size() && value > color_index_vec_[color_type][current_index].value_index ) current_index++;
    if ( current_index == 0 )
        return color_index_vec_[color_type][0].color;
    else {
        QColor color;
        if ( is_linear_[color_type] ){
            float scale = (color_index_vec_[color_type][current_index].value_index - value) / (color_index_vec_[color_type][current_index].value_index - color_index_vec_[color_type][current_index - 1].value_index);
            color.setRed(scale * color_index_vec_[color_type][current_index - 1].color.red() + (1.0 - scale) * color_index_vec_[color_type][current_index].color.red());
            color.setGreen(scale * color_index_vec_[color_type][current_index - 1].color.green() + (1.0 - scale) * color_index_vec_[color_type][current_index].color.green());
            color.setBlue(scale * color_index_vec_[color_type][current_index - 1].color.blue() + (1.0 - scale) * color_index_vec_[color_type][current_index].color.blue());
        } else {
            color = color_index_vec_[color_type][current_index - 1].color;
        }

        return color;
    }
}

void ColorMappingGenerator::GetColorIndex(ColorMappingType color_type, std::vector<ColorIndex >& index, bool& is_linear){
    index = color_index_vec_[color_type];
    is_linear = is_linear_[color_type];
}

void ColorMappingGenerator::ConstructColorList(){
    // gray 0 - 1
    ColorIndex gray_mapping[2];
    gray_mapping[0].value_index = 0.0;
    gray_mapping[0].color = QColor(0, 0, 0);
    gray_mapping[1].value_index = 1.0;
    gray_mapping[1].color = QColor(255, 255, 255);
    std::vector<ColorIndex > gray_index_vec;
    gray_index_vec.push_back(gray_mapping[0]);
    gray_index_vec.push_back(gray_mapping[1]);
    color_index_vec_.push_back(gray_index_vec);
    is_linear_.push_back(true);

    // uncertainty 0 - 1
    ColorIndex uncertainty_mapping[2];
    uncertainty_mapping[0].value_index = 0;
    uncertainty_mapping[0].color = QColor(255,255,255);
    uncertainty_mapping[1].value_index = 1.0;
    uncertainty_mapping[1].color = QColor(0,0,0);
    std::vector<ColorIndex > uncertainty_index_vec;
    for ( int i = 0; i < 2; ++i ) uncertainty_index_vec.push_back(uncertainty_mapping[i]);
    color_index_vec_.push_back(uncertainty_index_vec);
    is_linear_.push_back(true);

    // rain 0-100
    /*ColorIndex rain_mapping[7];
    rain_mapping[0].value_index = 0;
    rain_mapping[0].color = QColor(255, 255, 255);
    rain_mapping[1].value_index = 5;
    rain_mapping[1].color = QColor(170, 220, 170);
    rain_mapping[2].value_index = 10;
    rain_mapping[2].color = QColor(50, 190, 30);
    rain_mapping[3].value_index = 20;
    rain_mapping[3].color = QColor(20, 200, 240);
    rain_mapping[4].value_index = 30;
    rain_mapping[4].color = QColor(20, 30, 240);
    rain_mapping[5].value_index = 50;
    rain_mapping[5].color = QColor(220, 15, 240);
    rain_mapping[6].value_index = 200;
    rain_mapping[6].color = QColor(220, 15, 240);*/
    ColorIndex rain_mapping[7];
    rain_mapping[0].value_index = 0;
    rain_mapping[0].color = QColor(255, 255, 255);
    rain_mapping[1].value_index = 1;
    rain_mapping[1].color = QColor(170, 220, 170);
    rain_mapping[2].value_index = 5;
    rain_mapping[2].color = QColor(50, 190, 30);
    rain_mapping[3].value_index = 10;
    rain_mapping[3].color = QColor(20, 200, 240);
    rain_mapping[4].value_index = 25;
    rain_mapping[4].color = QColor(20, 30, 240);
    rain_mapping[5].value_index = 50;
    rain_mapping[5].color = QColor(220, 15, 240);
    rain_mapping[6].value_index = 100;
    rain_mapping[6].color = QColor(220, 15, 240);
    std::vector<ColorIndex > rain_index_vec;
    for ( int i = 0; i < 7; ++i ) rain_index_vec.push_back(rain_mapping[i]);
    color_index_vec_.push_back(rain_index_vec);
    is_linear_.push_back(false);

    ColorIndex tmp_mapping[10];
    tmp_mapping[0].value_index = -40 + 273.15;
    tmp_mapping[0].color = QColor(4,90,141);
    tmp_mapping[1].value_index = -30 + 273.15;
    tmp_mapping[1].color = QColor(43,140,190);
    tmp_mapping[2].value_index = -20 + 273.15;
    tmp_mapping[2].color = QColor(116,169,207);
    tmp_mapping[3].value_index = -10 + 273.15;
    tmp_mapping[3].color = QColor(189,201,225);
    tmp_mapping[4].value_index = 0 + 273.15;
    tmp_mapping[4].color = QColor(241,238,246);
    tmp_mapping[5].value_index = 10 + 273.15;
    tmp_mapping[5].color = QColor(254,240,217);
    tmp_mapping[6].value_index = 20 + 273.15;
    tmp_mapping[6].color = QColor(253,204,138);
    tmp_mapping[7].value_index = 30 + 273.15;
    tmp_mapping[7].color = QColor(252,141,89);
    tmp_mapping[8].value_index = 40 + 273.15;
    tmp_mapping[8].color = QColor(227,74,51);
    tmp_mapping[9].value_index = 50 + 273.15;
    tmp_mapping[9].color = QColor(179,0,0);
    std::vector<ColorIndex > tmp_index_vec;
    for ( int i = 0; i < 10; ++i ) tmp_index_vec.push_back(tmp_mapping[i]);
    color_index_vec_.push_back(tmp_index_vec);
    is_linear_.push_back(false);

    ColorIndex rms_mapping[9];
    /*rms_mapping[0].value_index = 0;
    rms_mapping[0].color = QColor(255,247,236);
    rms_mapping[1].value_index = 0.125;
    rms_mapping[1].color = QColor(254,232,200);
    rms_mapping[2].value_index = 0.25;
    rms_mapping[2].color = QColor(253,212,158);
    rms_mapping[3].value_index = 0.375;
    rms_mapping[3].color = QColor(253,187,132);
    rms_mapping[4].value_index = 0.5;
    rms_mapping[4].color = QColor(252,141,89);
    rms_mapping[5].value_index = 0.625;
    rms_mapping[5].color = QColor(239,101,72);
    rms_mapping[6].value_index = 0.75;
    rms_mapping[6].color = QColor(215,48,31);
    rms_mapping[7].value_index = 0.875;
    rms_mapping[7].color = QColor(153,0,0);
    rms_mapping[8].value_index = 1.0;
    rms_mapping[8].color = QColor(153,0,0);*/

    /*rms_mapping[0].value_index = 0;
    rms_mapping[0].color = QColor(237,248,233);
    rms_mapping[1].value_index = 0.2;
    rms_mapping[1].color = QColor(199,233,192);
    rms_mapping[2].value_index = 0.4;
    rms_mapping[2].color = QColor(161,217,155);
    rms_mapping[3].value_index = 0.6;
    rms_mapping[3].color = QColor(116,196,118);
    rms_mapping[4].value_index = 0.8;
    rms_mapping[4].color = QColor(65,171,93);
    rms_mapping[5].value_index = 1.0;
    rms_mapping[5].color = QColor(35,139,69);*/

    rms_mapping[0].value_index = 0;
    rms_mapping[0].color = QColor(239,243,255);
    rms_mapping[1].value_index = 2;
    rms_mapping[1].color = QColor(198,219,239);
    rms_mapping[2].value_index = 4;
    rms_mapping[2].color = QColor(158,202,225);
    rms_mapping[3].value_index = 6;
    rms_mapping[3].color = QColor(107,174,214);
    rms_mapping[4].value_index = 8;
    rms_mapping[4].color = QColor(66,146,198);
    rms_mapping[5].value_index = 10;
    rms_mapping[5].color = QColor(33,113,181);

    std::vector<ColorIndex > rms_index_vec;
    for ( int i = 0; i < 6; ++i ) rms_index_vec.push_back(rms_mapping[i]);
    color_index_vec_.push_back(rms_index_vec);
    is_linear_.push_back(true);

    ColorIndex prob_mapping[13];
    prob_mapping[0].value_index = 0;
    prob_mapping[0].color = QColor(255,255,255);
    prob_mapping[1].value_index = 0.1;
    prob_mapping[1].color = QColor(228,255,255);
    prob_mapping[2].value_index = 0.2;
    prob_mapping[2].color = QColor(196,232,255);
    prob_mapping[3].value_index = 0.3;
    prob_mapping[3].color = QColor(143,179,255);
    prob_mapping[4].value_index = 0.4;
    prob_mapping[4].color = QColor(216,249,216);
    prob_mapping[5].value_index = 0.5;
    prob_mapping[5].color = QColor(166,236,166);
    prob_mapping[6].value_index = 0.6;
    prob_mapping[6].color = QColor(66,247,66);
    prob_mapping[7].value_index = 0.7;
    prob_mapping[7].color = QColor(255,255,0);
    //prob_mapping[8].value_index = 0.7;
    //prob_mapping[8].color = QColor(255,215,0);
    prob_mapping[8].value_index = 0.8;
    prob_mapping[8].color = QColor(255,165,0);
    prob_mapping[9].value_index = 0.9;
    prob_mapping[9].color = QColor(246,163,174);
    prob_mapping[10].value_index = 0.95;
    prob_mapping[10].color = QColor(250,82,87);
    prob_mapping[11].value_index = 1.000;
    prob_mapping[11].color = QColor(255,0,0);
    std::vector<ColorIndex > prob_index_vec;
    for ( int i = 0; i < 12; ++i ) prob_index_vec.push_back(prob_mapping[i]);
    color_index_vec_.push_back(prob_index_vec);
    is_linear_.push_back(false);

	ColorIndex msl_pres_mapping[10];
	msl_pres_mapping[0].value_index = 100000;
	msl_pres_mapping[0].color = QColor(4,90,141);
	msl_pres_mapping[1].value_index = 400 + 100000;
	msl_pres_mapping[1].color = QColor(43,140,190);
	msl_pres_mapping[2].value_index = 800 + 100000;
	msl_pres_mapping[2].color = QColor(116,169,207);
	msl_pres_mapping[3].value_index = 1200 + 100000;
	msl_pres_mapping[3].color = QColor(189,201,225);
	msl_pres_mapping[4].value_index = 1600 + 100000;
	msl_pres_mapping[4].color = QColor(241,238,246);
	msl_pres_mapping[5].value_index = 2000 + 100000;
	msl_pres_mapping[5].color = QColor(254,240,217);
	msl_pres_mapping[6].value_index = 2400 + 100000;
	msl_pres_mapping[6].color = QColor(253,204,138);
	msl_pres_mapping[7].value_index = 2800 + 100000;
	msl_pres_mapping[7].color = QColor(252,141,89);
	msl_pres_mapping[8].value_index = 3200 + 100000;
	msl_pres_mapping[8].color = QColor(227,74,51);
	msl_pres_mapping[9].value_index = 3600 + 100000;
	msl_pres_mapping[9].color = QColor(179,0,0);
	std::vector<ColorIndex > msl_pres_index_vec;
	for ( int i = 0; i < 10; ++i ) msl_pres_index_vec.push_back(msl_pres_mapping[i]);
	color_index_vec_.push_back(msl_pres_index_vec);
	is_linear_.push_back(false);

	ColorIndex similarity_mapping[9];

    similarity_mapping[0].value_index = 0;
    similarity_mapping[0].color = QColor(239,243,255);
    similarity_mapping[1].value_index = 0.2;
    similarity_mapping[1].color = QColor(198,219,239);
    similarity_mapping[2].value_index = 0.4;
    similarity_mapping[2].color = QColor(158,202,225);
    similarity_mapping[3].value_index = 0.6;
    similarity_mapping[3].color = QColor(107,174,214);
    similarity_mapping[4].value_index = 0.8;
    similarity_mapping[4].color = QColor(66,146,198);
    similarity_mapping[5].value_index = 1.0;
    similarity_mapping[5].color = QColor(33,113,181);

    std::vector<ColorIndex > similarity_index_vec;
    for ( int i = 0; i < 6; ++i ) similarity_index_vec.push_back(similarity_mapping[i]);
    color_index_vec_.push_back(similarity_index_vec);
    is_linear_.push_back(true);

    ColorIndex variance_mapping[3];
    variance_mapping[0].value_index = 0.0;
    variance_mapping[0].color = QColor(229,245,224);
    variance_mapping[1].value_index = 0.5;
    variance_mapping[1].color = QColor(161,217,155);
    variance_mapping[2].value_index = 1.0;
    variance_mapping[2].color = QColor(49,163,84);
    std::vector<ColorIndex > variacen_index_vec;
    for ( int i = 0; i < 3; ++i ) variacen_index_vec.push_back(variance_mapping[i]);
    color_index_vec_.push_back(variacen_index_vec);
    is_linear_.push_back(true);


    /*qualitative_colors_.push_back(QColor(190, 174, 212));
    qualitative_colors_.push_back(QColor(178, 223, 138));
    qualitative_colors_.push_back(QColor(31, 107, 211));
    qualitative_colors_.push_back(QColor(255, 127, 0));*/

    qualitative_colors_.push_back(QColor(166,206,227));
    qualitative_colors_.push_back(QColor(31,120,180));
    qualitative_colors_.push_back(QColor(178,223,138));
    qualitative_colors_.push_back(QColor(51,160,44));
    qualitative_colors_.push_back(QColor(251,154,153));
    qualitative_colors_.push_back(QColor(227,26,28));
    qualitative_colors_.push_back(QColor(253,191,111));
    qualitative_colors_.push_back(QColor(255,127,0));
    qualitative_colors_.push_back(QColor(202,178,214));
    qualitative_colors_.push_back(QColor(106,61,154));
    qualitative_colors_.push_back(QColor(255,255,153));
    qualitative_colors_.push_back(QColor(177,89,40));
    qualitative_colors_.push_back(QColor(141,211,199));
    qualitative_colors_.push_back(QColor(204,235,197));
    qualitative_colors_.push_back(QColor(253,180,98));
}

void ColorMappingGenerator::GetQualitativeColors(int color_num, std::vector<QColor >& colors){
    colors.clear();
    for ( int i = 0; i < color_num && i < qualitative_colors_.size(); ++i )
        colors.push_back(qualitative_colors_[i]);
    for (int i = qualitative_colors_.size(); i < color_num; ++i)
        colors.push_back(QColor(128, 128, 128));
}

void ColorMappingGenerator::SetRmsMapping(float min_value, float max_value){
	for ( int i = 0; i < color_index_vec_[RMS_MAPPING].size(); ++i ){
		color_index_vec_[RMS_MAPPING][i].value_index = min_value + (max_value - min_value) * i / (color_index_vec_[RMS_MAPPING].size() - 1);
	}
}

void ColorMappingGenerator::GetRmsMapping(float& min_value, float& max_value){
	min_value = color_index_vec_[RMS_MAPPING][0].value_index;
	max_value = color_index_vec_[RMS_MAPPING][1].value_index;
}