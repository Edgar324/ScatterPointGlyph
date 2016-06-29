#ifndef WRF_NC_MAP_DATA_H_
#define WRF_NC_MAP_DATA_H_

#include "wrf_data_common.h"
#include <string>

class WrfNcMapData{
public:
    WrfNcMapData() {}
    ~WrfNcMapData() {}

    enum NcMapDataType{
        SQUARE_PROJECTION_DATA = 0x0,
        LAMBCC_PROJECTION_DATA
    };

    // @param time int time since 1980-1-1 0-0-0
    virtual bool CheckExisting(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num = 0) = 0;
    virtual bool CheckExisting(WrfModelType model) = 0;
    virtual void GetMaps(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<WrfGridValueMap*>& maps) = 0;
    virtual void GetData(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<float*>& data) = 0;
    virtual NcMapDataType DataType() = 0;

protected:
    std::string file_name_;
};

#endif