#ifndef WRF_DATA_MANAGER_H_
#define WRF_DATA_MANAGER_H_

#include <map>
#include <QtCore/QDateTime>
#include "wrf_data_common.h"
#include "wrf_nc_map_data.h"

class WrfNcLambCcProjectionData;
class WrfNcSquareProjectionData;

class WrfDataManager
{
public:
    static WrfDataManager* GetInstance();
    static bool DeleteInstance();

    void LoadDefaultData();

    void LoadEnsembleData(WrfElementType element, std::string& file_name);
    void LoadEnsembleMeanData(WrfElementType element, std::string& file_name);
    void LoadNarrData(WrfElementType element, std::string& file_name);
    void LoadChinaReanalysisData(WrfElementType element, std::string& file_name);
    void LoadCombinedData(std::string& file_name);

    void GetModelMapRange(WrfModelType model_type, MapRange& range);
    void GetEnsembleMapRange(MapRange& range);
    bool GetGridValueMap(QDateTime& datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector< WrfGridValueMap* >& maps);
    bool GetGridValueMap(int datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector< WrfGridValueMap* >& maps);
    bool GetGridData(QDateTime& datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector< float* >& data);
    bool GetGridData(int datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector< float* >& data);

    void GetEnsembleElements(std::vector< WrfElementType >& elements);
    void GetEnsembleMeanElements(std::vector< WrfElementType >& elements);
    void GetObservationElements(std::vector< WrfElementType >& elements);

    std::vector< std::vector< float > > GetMapInformation() { return border_polygons_; }

protected:
    WrfDataManager();
    ~WrfDataManager();

    static WrfDataManager* instance_;

private:
    std::vector< WrfNcMapData* > map_data_vec_;

    std::vector< WrfElementType > ensemble_elements_;
    std::vector< WrfElementType > ensemble_mean_elements_;
    std::vector< WrfElementType > observation_elements_;

    MapRange ens_map_range_;

    // map information
    std::vector< std::vector< float > > border_polygons_;

    void LoadMapInformation(const char* file_name, float start_x, float end_x, float start_y, float end_y);
};

#endif