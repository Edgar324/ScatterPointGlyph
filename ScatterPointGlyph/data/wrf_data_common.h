#ifndef WRF_DATA_COMMON_H_
#define WRF_DATA_COMMON_H_

#include <string>
#include <vector>

enum WrfModelType{
    WRF_UNKNOWN_MODEL = 0x0,
    WRF_NCEP_ENSEMBLES,
    WRF_NCEP_ENSEMBLE_MEAN,
	WRF_NCEP_ENSEMBLE_SPREAD,
    WRF_T639,
    WRF_EC_FINE,
    WRF_JAPAN_GSM,
    WRF_NCEP,
    WRF_REANALYSIS,
    WRF_ELEMENT_UNCERTAINTY,
    WRF_ELEMENT_PROBABILISTIC
};

enum WrfElementType{
    WRF_UNKNOWN_ELEMENT = 0x0,
    WRF_ACCUMULATED_PRECIPITATION,
    WRF_PRECIPITABLE_WATER,
    WRF_CAPE,
    WRF_CIN,
    WRF_T2M,
    WRF_MSLP,
	WRF_PRESSURE_SURFACE,
    WRF_TEMPERATURE_500HPA,
    WRF_HEIGHT_500HPA,
    WRF_RELATIVE_HUMIDITY_500HPA,
    WRF_TEMPERATURE_850HPA,
    WRF_UNCERTAINTY,
    WRF_PROBABILISTIC,
    WRF_RMS
};

enum WrfForecastingType{
    WRF_ENSEMBLER_PQPF = 0x0,
    WRF_REFORECAST_PQPF,
    WRF_ENSEMBLE_AVERAGE,
    WRF_ADAPTED_REFORECAST_PQPF,
    WRF_REALISTIC
};

struct MapRange{
    float start_x, start_y, end_x, end_y;
    float y_grid_space, x_grid_space;
    int y_grid_number, x_grid_number;
};

struct ProbabilisticPara{
    int time;
    WrfElementType element;
    float thresh;
    std::string para_name;
	int analog_number;
};

struct RmsUnit{
    float lon;
    float lat;
    int grid_point_index;
    std::vector< float > values;
};

class WrfValueMap
{
public:
    WrfValueMap() {
        model_type = WRF_UNKNOWN_MODEL;
        element_type = WRF_UNKNOWN_ELEMENT;
        value_per_grid = 1;
        is_updated = true;
    }
    ~WrfValueMap() {}

    enum ValueMapType{
        UNKNOWN_VALUE_MAP = 0x0,
        DISCRETE_VALUE_MAP,
        GRID_VALUE_MAP,
        LAMBCC_VALUE_MAP
    };

    virtual ValueMapType type(){
        return UNKNOWN_VALUE_MAP;
    }

    WrfModelType model_type;
    WrfElementType element_type;
    std::string title;
    int value_per_grid;
    bool is_updated;
};

class WrfDiscreteValueMap : public WrfValueMap
{
public:
    WrfDiscreteValueMap() 
        : dimension(0){
    }
    ~WrfDiscreteValueMap() {}

    int dimension;
    // x, y
    std::vector< float > coor;
    std::vector< float > values;

    virtual WrfValueMap::ValueMapType type() {
        return WrfValueMap::DISCRETE_VALUE_MAP;
    }
};

class WrfGridValueMap : public WrfValueMap
{
public:
    WrfGridValueMap(){
        map_range.x_grid_space = 0.0f;
        map_range.y_grid_space = 0.0f;
        map_range.start_x =  -1.0f;
        map_range.end_x = -1.0f;
        map_range.start_y = -1.0f;
        map_range.end_y = -1.0f;
        map_range.y_grid_number = 0;
        map_range.x_grid_space = 0;
        values = NULL;
    }

    WrfGridValueMap(const WrfGridValueMap& value_map){
        this->model_type = value_map.model_type;
        this->element_type = value_map.element_type;
        this->title = value_map.title;
        this->value_per_grid = value_map.value_per_grid;
        this->is_updated = value_map.is_updated;
        this->map_range = value_map.map_range;
        this->values = new float[value_map.value_per_grid * value_map.map_range.x_grid_number * value_map.map_range.y_grid_number];
        memcpy(this->values, value_map.values, value_map.value_per_grid * value_map.map_range.x_grid_number * value_map.map_range.y_grid_number * sizeof(float));
    }

    ~WrfGridValueMap() {}

    virtual WrfValueMap::ValueMapType type() {
        return WrfValueMap::GRID_VALUE_MAP;
    }

    MapRange map_range;
    float* values;

    const WrfGridValueMap& operator = (WrfGridValueMap& value_map){
        this->model_type = value_map.model_type;
        this->element_type = value_map.element_type;
        this->title = value_map.title;
        this->value_per_grid = value_map.value_per_grid;
        this->is_updated = value_map.is_updated;
        this->map_range = value_map.map_range;
        this->values = new float[value_map.value_per_grid * value_map.map_range.x_grid_number * value_map.map_range.y_grid_number];
        memcpy(this->values, value_map.values, value_map.value_per_grid * value_map.map_range.x_grid_number * value_map.map_range.y_grid_number * sizeof(float));
        return *this;
    }

    WrfGridValueMap* Copy(){
        WrfGridValueMap* map = new WrfGridValueMap;
        map->model_type = this->model_type;
        map->element_type = this->element_type;
        map->title = this->title;
        map->value_per_grid = this->value_per_grid;
        map->map_range = this->map_range;
        map->is_updated = this->is_updated;
        map->values = this->values;

        return map;
    }

    WrfGridValueMap* DeepCopy(){
        WrfGridValueMap* map = new WrfGridValueMap;
        map->model_type = this->model_type;
        map->element_type = this->element_type;
        map->title = this->title;
        map->value_per_grid = this->value_per_grid;
        map->map_range = this->map_range;
        map->is_updated = this->is_updated;
        map->values = new float[this->value_per_grid * this->map_range.x_grid_number * this->map_range.y_grid_number];
        memcpy(map->values, this->values, this->value_per_grid * this->map_range.x_grid_number * this->map_range.y_grid_number * sizeof(float));

        return map;
    }
};

class WrfLambCcValueMap : public WrfGridValueMap
{
public:
    WrfLambCcValueMap() {}
    ~WrfLambCcValueMap() {}

    virtual WrfValueMap::ValueMapType type() {
        return WrfValueMap::LAMBCC_VALUE_MAP;
    }
    
    std::vector< float > lon;
    std::vector< float > lat;

    double standard_parallel[2];
    double longitude_of_central_meridian;
    double latitude_of_projection_origin;
    double false_easting, false_northing;

    WrfLambCcValueMap* Copy(){
        WrfLambCcValueMap* map = new WrfLambCcValueMap;
        map->model_type = this->model_type;
        map->element_type = this->element_type;
        map->title = this->title;
        map->value_per_grid = this->value_per_grid;
        map->is_updated = this->is_updated;
        map->map_range = this->map_range;
        map->values = this->values;

        map->lon = this->lon;
        map->lat = this->lat;
        map->standard_parallel[0] = this->standard_parallel[0];
        map->standard_parallel[1] = this->standard_parallel[1];
        map->longitude_of_central_meridian = this->longitude_of_central_meridian;
        map->latitude_of_projection_origin = this->latitude_of_projection_origin;
        map->false_easting = this->false_easting;
        map->false_northing = this->false_northing;

        return map;
    }

    WrfLambCcValueMap* DeepCopy(){
        WrfLambCcValueMap* map = new WrfLambCcValueMap;
        map->model_type = this->model_type;
        map->element_type = this->element_type;
        map->title = this->title;
        map->value_per_grid = this->value_per_grid;
        map->is_updated = this->is_updated;
        map->map_range = this->map_range;
        map->values = new float[this->value_per_grid * this->map_range.x_grid_number * this->map_range.y_grid_number];
        memcpy(map->values, this->values, this->value_per_grid * this->map_range.x_grid_number * this->map_range.y_grid_number * sizeof(float));

        map->lon = this->lon;
        map->lat = this->lat;
        map->standard_parallel[0] = this->standard_parallel[0];
        map->standard_parallel[1] = this->standard_parallel[1];
        map->longitude_of_central_meridian = this->longitude_of_central_meridian;
        map->latitude_of_projection_origin = this->latitude_of_projection_origin;
        map->false_easting = this->false_easting;
        map->false_northing = this->false_northing;

        return map;
    }
};

class AxisData
{
public:
    AxisData(){
        name = "";
        analog_number = 0;
        r = 0;
        g = 0;
        b = 0;
        is_axis = false;
    }
    ~AxisData(){}

    std::string name;
    int analog_number;
    int r, g, b;
    bool is_axis;

    std::vector< WrfElementType > variables;
    std::vector< float > weights;
    std::vector< float > normalize_values;

    std::vector< float > rms_values;
    std::vector< int > sort_date_index;
    std::vector< int > date_sort_index;
};

#endif