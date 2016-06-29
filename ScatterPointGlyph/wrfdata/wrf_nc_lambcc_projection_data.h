#ifndef WRF_NC_LAMB_CC_PROJECTION_DATA_H_
#define WRF_NC_LAMB_CC_PROJECTION_DATA_H_

#include <vector>
#include <map>
#include "wrf_nc_map_data.h"
#include "wrf_data_common.h"

class WrfNcLambCcProjectionData : public WrfNcMapData{
public:
    WrfNcLambCcProjectionData();
    WrfNcLambCcProjectionData(WrfNcLambCcProjectionData& data);
    WrfNcLambCcProjectionData(std::string& file_name);
    ~WrfNcLambCcProjectionData();

    virtual bool CheckExisting(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);
    virtual bool CheckExisting(WrfModelType model);

    virtual float* GetData(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);
    WrfLambCcValueMap* GetMap(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);

    virtual void GetMaps(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<WrfGridValueMap*>& maps);
    virtual void GetData(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<float*>& data);

    virtual NcMapDataType DataType() { return WrfNcMapData::LAMBCC_PROJECTION_DATA; }

    void SaveData(const char* file_name);
    WrfNcLambCcProjectionData* Convert2Range(float start_lon, float end_lon, float start_lat, float end_lat);

    WrfModelType model_;
    WrfElementType element_;
    std::vector<float> lon_;
    std::vector<float> lat_;
    std::vector<float> x_dist_;
    std::vector<float> y_dist_;
    double standard_parallel_[2];
    double longitude_of_central_meridian_;
    double latitude_of_projection_origin_;
    double false_easting_, false_northing_;

    std::map<int, int> time_map_;
    size_t time_len_, y_len_, x_len_, nbnds_len_;

    // [time][lat][lon]
    float* data_values_;

    MapRange map_range_;

private:
    WrfLambCcValueMap* map_;

    void LoadData(const char* file_name);
    void LoadTotalData(const char* file_name);
    void LoadConvertedData(const char* file_name);
};

#endif