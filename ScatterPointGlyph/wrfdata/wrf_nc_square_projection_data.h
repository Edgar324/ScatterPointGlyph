#ifndef WRF_NC_SQUARE_PROJECTION_DATA_H_
#define WRF_NC_SQUARE_PROJECTION_DATA_H_

#include "wrf_nc_map_data.h"
#include <vector>
#include <map>
#include "wrf_data_common.h"

class WrfNcSquareProjectionData : public WrfNcMapData{
public:
    WrfNcSquareProjectionData();
    WrfNcSquareProjectionData(std::string& file_name, WrfModelType model, WrfElementType element);
    WrfNcSquareProjectionData(std::string& file_name, WrfModelType model, WrfElementType element, int year);
    ~WrfNcSquareProjectionData();

    virtual bool CheckExisting(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);
    virtual bool CheckExisting(WrfModelType model);

    float* GetData(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);
    WrfGridValueMap* GetMap(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0);

    virtual void GetMaps(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<WrfGridValueMap*>& maps);
    virtual void GetData(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<float*>& data);

    virtual NcMapDataType DataType() { return WrfNcMapData::SQUARE_PROJECTION_DATA; }

    MapRange map_range_;

    std::map<int, int> time_map_;
    std::map<int, int> fhour_map_;
    WrfModelType model_;
    WrfElementType element_;
    std::vector<float> lat_;
    std::vector<float> lon_;
    size_t time_len_, lat_len_, lon_len_, ens_len_, fhour_len_;

    // [time][ens][lat][lon]
    float* data_values_;

private:
    void LoadData(const char* file_name);
};

#endif