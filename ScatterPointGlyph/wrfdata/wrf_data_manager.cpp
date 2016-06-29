#include "wrf_data_manager.h"
#include "wrf_data_common.h"
#include "netcdf.h"
#include <iostream>
#include <fstream>
#include <QtWidgets/QMessageBox>
#include "wrf_nc_square_projection_data.h"
#include "wrf_nc_lambcc_projection_data.h"
#include "wrf_utility.h"

WrfDataManager* WrfDataManager::instance_ = NULL;

WrfDataManager* WrfDataManager::GetInstance(){
    if ( instance_ == NULL ){
        instance_ = new WrfDataManager;
    }

    return instance_;
}

bool WrfDataManager::DeleteInstance(){
    if ( instance_ != NULL ){
        delete instance_;
        return true;
    } else {
        return false;
    }
}

WrfDataManager::WrfDataManager(){
}

WrfDataManager::~WrfDataManager(){
    // clear all the data existed
}

void WrfDataManager::LoadDefaultData(){
    //WrfNcSquareProjectionData* data = new WrfNcSquareProjectionData(std::string("E:\\Weather Forcast\\Data\\apcp_sfc_latlon_all_19841201_20041201_liaoD0wwUX.nc"));
    //ShowNcInfo("E:\\Weather Forcast\\Data\\apcp_sfc_latlon_all_19841201_20041201_liaoD0wwUX.nc");
    //ShowNcInfo("E:\\Weather Forcast\\Data\\apcp.1979.nc");
    //WrfNcLambCcProjectionData* lambcc_data = new WrfNcLambCcProjectionData(std::string("E:\\Weather Forcast\\Data\\apcp.1979.nc"));
    //ShowNcInfo("E:\\Weather Forcast\\Data\\border_lines\\gshhg-gmt-2.3.3\\binned_border_f.nc");
    //LoadMapInformation("E:\\Weather Forcast\\Data\\border_lines\\gshhg-gmt-2.3.3\\binned_border_f.nc");
    //ShowNcInfo("E:\\Weather Forcast\\Data\\tmp_2m_latlon_all_19841201_20041201_liaochoRjy.nc");
    //ShowNcInfo("E:/Weather Forcast/Data/apcp_sfc_latlon_all_20141101_20141105.nc");
    //QDateTime temp(QDate(1800, 1,1), QTime(0, 0, 0));
    //qint64 hour = temp.msecsTo(QDateTime(QDate(1, 1, 1), QTime(0, 0, 0))) / 3600000;
    //ShowNcInfo("./Data/apcpusdaily_1979010100-2005013100.nc");
    //LoadData(WRF_NCEP_ENSEMBLE0, WRF_RAIN_24_HOURS, std::string(".\\Data\\apcp_sfc_latlon_all_19841201_20041201_liaoD0wwUX.nc"));
    //LoadData(WRF_NCEP_ENSEMBLE0, WRF_TEMPERATURE_500HPA, std::string(".\\Data\\tmp_2m_latlon_all_19841201_20041201_liaochoRjy.nc"));
    //ShowNcInfo("./Data/from_tom/refcstv2_precip_ccpav2_000_to_012.nc");
    //LoadData(WRF_NCEP_ENSEMBLE0, WRF_RAIN, std::string(".\\Data\\china\\apcp_sfc_latlon_all_19850101_20141116_liaoHVLxK3.nc"));
    //LoadData(WRF_NCEP_MEAN, WRF_RAIN_24_HOURS, std::string(".\\Data\\china\\apcp_sfc_latlon_mean_19850101_20141116_liaoCntTBY.nc"));
    //LoadData(WRF_NCEP_MEAN, WRF_TEMPERATURE_500HPA, std::string(".\\Data\\china\\pwat_eatm_latlon_mean_19850101_20141116_liaoIEhsjq.nc"));
    //WrfNcLambCcProjectionData* lambcc_data = new WrfNcLambCcProjectionData(std::string("E:\\Weather Forcast\\Data\\apcp.1979.nc"));
    //WrfNcLambCcProjectionData* new_data = lambcc_data->Convert2Range(-127, -65, 23, 51);
    //new_data->SaveData("test.nc");
    //LoadData(WRF_OBSERVATION, WRF_RAIN_24_HOURS, std::string("./Data/apcpusdaily_1979010100-2005013100.nc"));

    this->LoadEnsembleData(WRF_ACCUMULATED_PRECIPITATION, std::string(".\\Data\\china\\apcp_sfc_latlon_all_19850101_20141116_liaoHVLxK3.nc"));
    //this->LoadEnsembleMeanData(WRF_ACCUMULATED_PRECIPITATION, std::string(".\\Data\\china\\apcp_sfc_latlon_mean_19850101_20141116_liaoCntTBY.nc"));
    //this->LoadCombinedData(std::string("./Data/from_tom/refcstv2_precip_ccpav2_012_to_024.nc"));
    //this->LoadNarrData(WRF_ACCUMULATED_PRECIPITATION, std::string("./Data/apcpusdaily_1979010100-2005013100.nc"));
    //this->LoadNarrData(WRF_RAIN, std::string("E:\\Weather Forcast\\Data\\apcp.1979.nc"));
    for ( int i = 2000; i <= 2013; ++i ){
        QString str = QString("./Data/china_reanalysis/SURF_CLI_CHN_PRE_DAY_GRID_0.25-%0").arg(i);
        this->LoadChinaReanalysisData(WRF_ACCUMULATED_PRECIPITATION, std::string(str.toLocal8Bit().data()));
    }
}

void WrfDataManager::LoadEnsembleData(WrfElementType element, std::string& file_name){
    WrfNcSquareProjectionData* square_data = new WrfNcSquareProjectionData(file_name, WRF_NCEP_ENSEMBLES, element);
    this->map_data_vec_.push_back(square_data);
    ensemble_elements_.push_back(element);

    if ( ensemble_elements_.size() == 1 && ensemble_mean_elements_.size() == 0 ){
        ens_map_range_ = square_data->map_range_;
        LoadMapInformation("./Resources/border.txt", ens_map_range_.start_x, ens_map_range_.end_x, ens_map_range_.start_y, ens_map_range_.end_y);
    }
}

void WrfDataManager::LoadEnsembleMeanData(WrfElementType element, std::string& file_name){
    WrfNcSquareProjectionData* square_data = new WrfNcSquareProjectionData(file_name, WRF_NCEP_ENSEMBLE_MEAN, element);
    this->map_data_vec_.push_back(square_data);
    ensemble_mean_elements_.push_back(element);

    if ( ensemble_mean_elements_.size() == 1 && ensemble_elements_.size() == 0 ){
        ens_map_range_ = square_data->map_range_;
        LoadMapInformation("./Resources/border.txt", ens_map_range_.start_x, ens_map_range_.end_x, ens_map_range_.start_y, ens_map_range_.end_y);
    }
}

void WrfDataManager::LoadNarrData(WrfElementType element, std::string& file_name){
     WrfNcLambCcProjectionData* data = new WrfNcLambCcProjectionData(file_name);
     this->map_data_vec_.push_back(data);
     observation_elements_.push_back(element);
}

void WrfDataManager::LoadChinaReanalysisData(WrfElementType element, std::string& file_name){
    std::string year_string = file_name.substr(file_name.length() - 4, 4);
    WrfNcSquareProjectionData* data = new WrfNcSquareProjectionData(file_name, WRF_REANALYSIS, element, atoi(year_string.c_str()));
    this->map_data_vec_.push_back(data);

    bool is_exist = false;
    for ( int i = 0; i <observation_elements_.size(); ++i )
        if ( data->element_ == observation_elements_[i] ){
            is_exist = true;
            break;
        }
    if ( !is_exist ) observation_elements_.push_back(data->element_);
}

void WrfDataManager::GetEnsembleElements(std::vector<WrfElementType>& elements){
    elements = ensemble_elements_;
}

void WrfDataManager::GetEnsembleMeanElements(std::vector<WrfElementType>& elements){
    elements = ensemble_mean_elements_;
}

void WrfDataManager::GetObservationElements(std::vector<WrfElementType>& elements){
    elements = observation_elements_;
}

bool WrfDataManager::GetGridValueMap(QDateTime& datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector<WrfGridValueMap*>& maps){
    bool is_exist = false;

    qint64 time = datetime.msecsTo(QDateTime(QDate(1800, 1, 1), QTime(0, 0, 0))) / 3600000 * -1;
    for ( int i = 0; i <map_data_vec_.size(); ++i )
        if ( map_data_vec_[i]->CheckExisting((int)time, model_type, element_type, fhour) ){
            map_data_vec_[i]->GetMaps(time, model_type, element_type, fhour, maps);
            is_exist = true;

            break;
        }

    return is_exist;
}

bool WrfDataManager::GetGridValueMap(int datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector<WrfGridValueMap*>& maps){
    bool is_exist = false;

    for ( int i = 0; i <map_data_vec_.size(); ++i )
        if ( map_data_vec_[i]->CheckExisting(datetime, model_type, element_type, fhour) ){
            map_data_vec_[i]->GetMaps(datetime, model_type, element_type, fhour, maps);
            is_exist = true;

            break;
        }

    return is_exist;
}

bool WrfDataManager::GetGridData(QDateTime& datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector<float*>& data){
    bool is_exist = false;

    qint64 time = datetime.msecsTo(QDateTime(QDate(1800, 1, 1), QTime(0, 0, 0))) / 3600000 * -1;
    for ( int i = 0; i <map_data_vec_.size(); ++i )
        if ( map_data_vec_[i]->CheckExisting((int)time, model_type, element_type, fhour) ){
            map_data_vec_[i]->GetData(time, model_type, element_type, fhour, data);
            is_exist = true;

            break;
        }

    return is_exist;
}

bool WrfDataManager::GetGridData(int datetime, WrfModelType model_type, WrfElementType element_type, int fhour, std::vector<float*>& data){
    bool is_exist = false;

    for ( int i = 0; i <map_data_vec_.size(); ++i )
        if ( map_data_vec_[i]->CheckExisting(datetime, model_type, element_type, fhour) ){
            map_data_vec_[i]->GetData(datetime, model_type, element_type, fhour, data);
            is_exist = true;

            break;
        }

    return is_exist;
}

void WrfDataManager::GetModelMapRange(WrfModelType model_type, MapRange& range){
    for ( int i = 0; i <map_data_vec_.size(); ++i )
        if ( map_data_vec_[i]->CheckExisting(model_type) ){
            switch ( map_data_vec_[i]->DataType() ){
            case WrfNcMapData::SQUARE_PROJECTION_DATA:{
                WrfNcSquareProjectionData* nc_data = dynamic_cast<WrfNcSquareProjectionData*>(map_data_vec_[i]);
                range = nc_data->map_range_;
                break;
                                                      }
            case WrfNcMapData::LAMBCC_PROJECTION_DATA:{
                WrfNcLambCcProjectionData* nc_data = dynamic_cast<WrfNcLambCcProjectionData*>(map_data_vec_[i]);
                range = nc_data->map_range_;
                break;
                                                      }
            default:{
                std::cout <<"No map range for the selected model_type: " <<enum_model_to_string(model_type) <<std::endl; 
                exit(-1);
                    }                 
            }
            break;
        }
}

void WrfDataManager::GetEnsembleMapRange(MapRange& range){
    range = ens_map_range_;
}

void WrfDataManager::LoadMapInformation(const char* file_name, float start_x, float end_x, float start_y, float end_y){
    border_polygons_.clear();

    std::ifstream border_file(file_name);
    if ( border_file.good() ){
        while ( !border_file.eof() ){
            int poly_size;
            float x, y;
            std::vector<float> temp_poly;
            border_file>> poly_size;
            temp_poly.resize(poly_size * 2);
            bool is_negative = false;
            for ( int i = 0; i <poly_size; ++i ){
                border_file>> temp_poly[2 * i]>> temp_poly[2 * i + 1];
                if ( temp_poly[2 * i] <0 ) is_negative = true;
            }
            if ( is_negative )
                for ( int i = 0; i <temp_poly.size() / 2; ++i ) temp_poly[2 * i] += 360;
            bool is_used = false;
            for ( int i = 0; i <temp_poly.size() / 2; ++i )
                if ( (temp_poly[2 * i] - start_x) * (temp_poly[2 * i] - end_x) <= 0 && (temp_poly[2 * i +1] - start_y) * (temp_poly[2 * i + 1] - end_y) <= 0 ){
                    is_used = true;
                    break;
                }
            if ( is_used ) border_polygons_.push_back(temp_poly);
        }
    }
}

void WrfDataManager::LoadCombinedData(std::string& file_name){
    int ncid, dim_id, var_id;
    int status;

    status = nc_open(file_name.c_str(), NC_NOWRITE, &ncid);

    if ( status != NC_NOERR ) return;

    size_t xa_len, ya_len, xf_len, yf_len, time_len, ens_len;
    std::vector<float> xa, ya, xf, yf, time;
    std::vector<int> ensv;
    std::vector<float> lons_anal, lats_anal, lons_fcst, lats_fcst;
    std::vector<int> init_time, fcst_time_begin, fcst_time_end;

    nc_inq_dimid(ncid, "xa", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &xa_len);

    nc_inq_dimid(ncid, "ya", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &ya_len);

    nc_inq_dimid(ncid, "xf", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &xf_len);

    nc_inq_dimid(ncid, "yf", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &yf_len);

    nc_inq_dimid(ncid, "time", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &time_len);

    nc_inq_dimid(ncid, "ens", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &ens_len);

    nc_inq_varid(ncid, "xa", &var_id);
    xa.resize(xa_len);
    nc_get_var_float(ncid, var_id, xa.data());

    nc_inq_varid(ncid, "ya", &var_id);
    ya.resize(ya_len);
    nc_get_var_float(ncid, var_id, ya.data());

    nc_inq_varid(ncid, "xf", &var_id);
    xf.resize(xf_len);
    nc_get_var_float(ncid, var_id, xf.data());

    nc_inq_varid(ncid, "yf", &var_id);
    yf.resize(yf_len);
    nc_get_var_float(ncid, var_id, yf.data());

    nc_inq_varid(ncid, "time", &var_id);
    time.resize(time_len);
    nc_get_var_float(ncid, var_id, time.data());

    qint64 temp_time = QDateTime(QDate(1, 1, 1), QTime(0, 0, 0)).msecsTo(QDateTime(QDate(1800, 1, 1), QTime(0, 0, 0))) / 3600000;
    for ( int i = 0; i <time.size(); ++i ) time[i] -= temp_time;

    nc_inq_varid(ncid, "ensv", &var_id);
    ensv.resize(ens_len);
    nc_get_var_int(ncid, var_id, ensv.data());

    nc_inq_varid(ncid, "lons_anal", &var_id);
    lons_anal.resize(ya_len * xa_len);
    nc_get_var_float(ncid, var_id, lons_anal.data());

    nc_inq_varid(ncid, "lats_anal", &var_id);
    lats_anal.resize(ya_len * xa_len);
    nc_get_var_float(ncid, var_id, lats_anal.data());

    nc_inq_varid(ncid, "lons_fcst", &var_id);
    lons_fcst.resize(yf_len * xf_len);
    nc_get_var_float(ncid, var_id, lons_fcst.data());

    nc_inq_varid(ncid, "lats_fcst", &var_id);
    lats_fcst.resize(yf_len * xf_len);
    nc_get_var_float(ncid, var_id, lats_fcst.data());

    nc_inq_varid(ncid, "yyyymmddhh_init", &var_id);
    init_time.resize(time_len);
    nc_get_var_int(ncid, var_id, init_time.data());

    nc_inq_varid(ncid, "yyyymmddhh_fcstb", &var_id);
    fcst_time_begin.resize(time_len);
    nc_get_var_int(ncid, var_id, fcst_time_begin.data());

    nc_inq_varid(ncid, "yyyymmddhh_fcste", &var_id);
    fcst_time_end.resize(time_len);
    nc_get_var_int(ncid, var_id, fcst_time_end.data());

    // pre process data
    std::map<int, int> ens_time_map, anal_time_map;
    std::map<int, int> fhour_map;
    std::vector<float> lat_anal_vec, lon_anal_vec, lat_fcst_vec, lon_fcst_vec;
    MapRange anal_map_range, fcst_map_range;

    int bias_hour = (fcst_time_end[0] - init_time[0]) / 100 * 24 + (fcst_time_end[0] - init_time[0]) % 100;
    for ( int i = 0; i <time.size(); ++i ) {
        ens_time_map.insert(std::map<int, int>::value_type((int)(time[i] - bias_hour), i));
    }
    fhour_map.insert(std::map<int, int>::value_type(bias_hour, 0));
    lat_anal_vec.resize(ya_len);
    for ( int i = 0; i <ya_len; ++i ) lat_anal_vec[i] = lats_anal[xa_len * i];
    lon_anal_vec.resize(xa_len);
    for ( int i = 0; i <xa_len; ++i ) lon_anal_vec[i] = lons_anal[i] + 360;
    lat_fcst_vec.resize(yf_len);
    for ( int i = 0; i <yf_len; ++i ) lat_fcst_vec[i] = lats_fcst[xf_len * i];
    lon_fcst_vec.resize(xf_len);
    for ( int i = 0; i <xf_len; ++i ) lon_fcst_vec[i] = lons_fcst[i] + 360;

    anal_map_range.start_x = lon_anal_vec[0];
    anal_map_range.end_x = lon_anal_vec[lon_anal_vec.size() - 1];
    anal_map_range.start_y = lat_anal_vec[0];
    anal_map_range.end_y = lat_anal_vec[lat_anal_vec.size() - 1];
    anal_map_range.x_grid_space = (lon_anal_vec[lon_anal_vec.size() - 1] - lon_anal_vec[0]) / (lon_anal_vec.size() - 1);
    anal_map_range.y_grid_space = (lat_anal_vec[lat_anal_vec.size() - 1] - lat_anal_vec[0]) / (lat_anal_vec.size() - 1);
    anal_map_range.x_grid_number = xa_len;
    anal_map_range.y_grid_number = ya_len;

    fcst_map_range.start_x = lon_fcst_vec[0];
    fcst_map_range.end_x = lon_fcst_vec[lon_fcst_vec.size() - 1];
    fcst_map_range.start_y = lat_fcst_vec[0];
    fcst_map_range.end_y = lat_fcst_vec[lat_fcst_vec.size() - 1];
    fcst_map_range.x_grid_space = (lon_fcst_vec[lon_fcst_vec.size() - 1] - lon_fcst_vec[0]) / (lon_fcst_vec.size() - 1);
    fcst_map_range.y_grid_space = (lat_fcst_vec[lat_fcst_vec.size() - 1] - lat_fcst_vec[0]) / (lat_fcst_vec.size() - 1);
    fcst_map_range.x_grid_number = xf_len;
    fcst_map_range.y_grid_number = yf_len;

	std::cout <<"X Space: "<<fcst_map_range.x_grid_space <<"  Y Space: "<<fcst_map_range.y_grid_space <<std::endl;

    ens_map_range_ = fcst_map_range;
    LoadMapInformation("./Resources/border.txt", ens_map_range_.start_x, ens_map_range_.end_x, ens_map_range_.start_y, ens_map_range_.end_y);

    // apcp_anal
    WrfNcSquareProjectionData* apcp_anal_data = new WrfNcSquareProjectionData;
    apcp_anal_data->map_range_ = anal_map_range;
    apcp_anal_data->time_map_ = ens_time_map;
    apcp_anal_data->fhour_map_ = fhour_map;
    apcp_anal_data->lat_ = lat_anal_vec;
    apcp_anal_data->lon_ = lon_anal_vec;
    apcp_anal_data->time_len_ = time_len;
    apcp_anal_data->lat_len_ = ya_len;
    apcp_anal_data->lon_len_ = xa_len;
    apcp_anal_data->ens_len_ = 1;
    apcp_anal_data->fhour_len_ = 1;
    apcp_anal_data->model_ = WRF_REANALYSIS;
    apcp_anal_data->element_ = WRF_ACCUMULATED_PRECIPITATION;
    apcp_anal_data->data_values_ = new float[time_len * ya_len * xa_len];
    nc_inq_varid(ncid, "apcp_anal", &var_id);
    nc_get_var_float(ncid, var_id, apcp_anal_data->data_values_);
    for ( int i = 0; i <time_len * ya_len * xa_len; ++i )
        if ( apcp_anal_data->data_values_[i] <-1 ) apcp_anal_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(apcp_anal_data);
    observation_elements_.push_back(WRF_ACCUMULATED_PRECIPITATION);

    // apcp_fcst_ens
    WrfNcSquareProjectionData* apcp_fcst_ens_data = new WrfNcSquareProjectionData;
    apcp_fcst_ens_data->map_range_ = fcst_map_range;
    apcp_fcst_ens_data->time_map_ = ens_time_map;
    apcp_fcst_ens_data->fhour_map_ = fhour_map;
    apcp_fcst_ens_data->lat_ = lat_fcst_vec;
    apcp_fcst_ens_data->lon_ = lon_fcst_vec;
    apcp_fcst_ens_data->time_len_ = time_len;
    apcp_fcst_ens_data->lat_len_ = yf_len;
    apcp_fcst_ens_data->lon_len_ = xf_len;
    apcp_fcst_ens_data->ens_len_ = ens_len;
    apcp_fcst_ens_data->fhour_len_ = 1;
    apcp_fcst_ens_data->model_ = WRF_NCEP_ENSEMBLES;
    apcp_fcst_ens_data->element_ = WRF_ACCUMULATED_PRECIPITATION;
    apcp_fcst_ens_data->data_values_ = new float[time_len * ens_len * yf_len * xf_len];
    nc_inq_varid(ncid, "apcp_fcst_ens", &var_id);
    nc_get_var_float(ncid, var_id, apcp_fcst_ens_data->data_values_);
    for ( int i = 0; i <time_len * ens_len * yf_len * xf_len; ++i )
        if ( apcp_fcst_ens_data->data_values_[i] <-1 ) apcp_fcst_ens_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(apcp_fcst_ens_data);
    ensemble_elements_.push_back(WRF_ACCUMULATED_PRECIPITATION);

    // T2m_fcst_mean
    WrfNcSquareProjectionData* t2m_fcst_mean_data = new WrfNcSquareProjectionData;
    t2m_fcst_mean_data->map_range_ = fcst_map_range;
    t2m_fcst_mean_data->time_map_ = ens_time_map;
    t2m_fcst_mean_data->fhour_map_ = fhour_map;
    t2m_fcst_mean_data->lat_ = lat_fcst_vec;
    t2m_fcst_mean_data->lon_ = lon_fcst_vec;
    t2m_fcst_mean_data->time_len_ = time_len;
    t2m_fcst_mean_data->lat_len_ = yf_len;
    t2m_fcst_mean_data->lon_len_ = xf_len;
    t2m_fcst_mean_data->ens_len_ = 1;
    t2m_fcst_mean_data->fhour_len_ = 1;
    t2m_fcst_mean_data->model_ = WRF_NCEP_ENSEMBLE_MEAN;
    t2m_fcst_mean_data->element_ = WRF_T2M;
    t2m_fcst_mean_data->data_values_ = new float[time_len * yf_len * xf_len];
    nc_inq_varid(ncid, "T2m_fcst_mean", &var_id);
    nc_get_var_float(ncid, var_id, t2m_fcst_mean_data->data_values_);
    for ( int i = 0; i <time_len * yf_len * xf_len; ++i )
        if ( t2m_fcst_mean_data->data_values_[i] <-1 ) t2m_fcst_mean_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(t2m_fcst_mean_data);
    ensemble_mean_elements_.push_back(WRF_T2M);

    // pwat_fcst_mean
    WrfNcSquareProjectionData* pwat_fcst_mean_data = new WrfNcSquareProjectionData;
    pwat_fcst_mean_data->map_range_ = fcst_map_range;
    pwat_fcst_mean_data->time_map_ = ens_time_map;
    pwat_fcst_mean_data->fhour_map_ = fhour_map;
    pwat_fcst_mean_data->lat_ = lat_fcst_vec;
    pwat_fcst_mean_data->lon_ = lon_fcst_vec;
    pwat_fcst_mean_data->time_len_ = time_len;
    pwat_fcst_mean_data->lat_len_ = yf_len;
    pwat_fcst_mean_data->lon_len_ = xf_len;
    pwat_fcst_mean_data->ens_len_ = 1;
    pwat_fcst_mean_data->fhour_len_ = 1;
    pwat_fcst_mean_data->model_ = WRF_NCEP_ENSEMBLE_MEAN;
    pwat_fcst_mean_data->element_ = WRF_PRECIPITABLE_WATER;
    pwat_fcst_mean_data->data_values_ = new float[time_len * yf_len * xf_len];
    nc_inq_varid(ncid, "pwat_fcst_mean", &var_id);
    nc_get_var_float(ncid, var_id, pwat_fcst_mean_data->data_values_);
    for ( int i = 0; i <time_len * yf_len * xf_len; ++i )
        if ( pwat_fcst_mean_data->data_values_[i] <-1 ) pwat_fcst_mean_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(pwat_fcst_mean_data);
    ensemble_mean_elements_.push_back(WRF_PRECIPITABLE_WATER);

    // cape_fcst_mean
    /*WrfNcSquareProjectionData* cape_fcst_mean_data = new WrfNcSquareProjectionData;
    cape_fcst_mean_data->map_range_ = fcst_map_range;
    cape_fcst_mean_data->time_map_ = ens_time_map;
    cape_fcst_mean_data->fhour_map_ = fhour_map;
    cape_fcst_mean_data->lat_ = lat_fcst_vec;
    cape_fcst_mean_data->lon_ = lon_fcst_vec;
    cape_fcst_mean_data->time_len_ = time_len;
    cape_fcst_mean_data->lat_len_ = yf_len;
    cape_fcst_mean_data->lon_len_ = xf_len;
    cape_fcst_mean_data->ens_len_ = 1;
    cape_fcst_mean_data->fhour_len_ = 1;
    cape_fcst_mean_data->model_ = WRF_NCEP_ENSEMBLE_MEAN;
    cape_fcst_mean_data->element_ = WRF_CAPE;
    cape_fcst_mean_data->data_values_ = new float[time_len * yf_len * xf_len];
    nc_inq_varid(ncid, "cape_fcst_mean", &var_id);
    nc_get_var_float(ncid, var_id, cape_fcst_mean_data->data_values_);
    for ( int i = 0; i <time_len * yf_len * xf_len; ++i )
        if ( cape_fcst_mean_data->data_values_[i] <-1 ) cape_fcst_mean_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(cape_fcst_mean_data);
    ensemble_mean_elements_.push_back(WRF_CAPE);*/

    // cin_fcst_mean
    /*WrfNcSquareProjectionData* cin_fcst_mean_data = new WrfNcSquareProjectionData;
    cin_fcst_mean_data->map_range_ = fcst_map_range;
    cin_fcst_mean_data->time_map_ = ens_time_map;
    cin_fcst_mean_data->fhour_map_ = fhour_map;
    cin_fcst_mean_data->lat_ = lat_fcst_vec;
    cin_fcst_mean_data->lon_ = lon_fcst_vec;
    cin_fcst_mean_data->time_len_ = time_len;
    cin_fcst_mean_data->lat_len_ = yf_len;
    cin_fcst_mean_data->lon_len_ = xf_len;
    cin_fcst_mean_data->ens_len_ = 1;
    cin_fcst_mean_data->fhour_len_ = 1;
    cin_fcst_mean_data->model_ = WRF_NCEP_ENSEMBLE_MEAN;
    cin_fcst_mean_data->element_ = WRF_CIN;
    cin_fcst_mean_data->data_values_ = new float[time_len * yf_len * xf_len];
    nc_inq_varid(ncid, "cin_fcst_mean", &var_id);
    nc_get_var_float(ncid, var_id, cin_fcst_mean_data->data_values_);
    for ( int i = 0; i <time_len * yf_len * xf_len; ++i )
        if ( cin_fcst_mean_data->data_values_[i] <-4000 ) cin_fcst_mean_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(cin_fcst_mean_data);
    ensemble_mean_elements_.push_back(WRF_CIN);*/

    //// MSLP_fcst_mean
    WrfNcSquareProjectionData* mslp_fcst_mean_data = new WrfNcSquareProjectionData;
    mslp_fcst_mean_data->map_range_ = fcst_map_range;
    mslp_fcst_mean_data->time_map_ = ens_time_map;
    mslp_fcst_mean_data->fhour_map_ = fhour_map;
    mslp_fcst_mean_data->lat_ = lat_fcst_vec;
    mslp_fcst_mean_data->lon_ = lon_fcst_vec;
    mslp_fcst_mean_data->time_len_ = time_len;
    mslp_fcst_mean_data->lat_len_ = yf_len;
    mslp_fcst_mean_data->lon_len_ = xf_len;
    mslp_fcst_mean_data->ens_len_ = 1;
    mslp_fcst_mean_data->fhour_len_ = 1;
    mslp_fcst_mean_data->model_ = WRF_NCEP_ENSEMBLE_MEAN;
    mslp_fcst_mean_data->element_ = WRF_MSLP;
    mslp_fcst_mean_data->data_values_ = new float[time_len * yf_len * xf_len];
    nc_inq_varid(ncid, "MSLP_fcst_mean", &var_id);
    nc_get_var_float(ncid, var_id, mslp_fcst_mean_data->data_values_);
    for ( int i = 0; i <time_len * yf_len * xf_len; ++i )
        if ( mslp_fcst_mean_data->data_values_[i] <-1 ) mslp_fcst_mean_data->data_values_[i] = -10001;
    this->map_data_vec_.push_back(mslp_fcst_mean_data);
    ensemble_mean_elements_.push_back(WRF_MSLP);

    nc_close(ncid);

    QMessageBox::information(0, QString("Time range"), QString("Time Range: From %0 to %1").arg(init_time[0]).arg(init_time[init_time.size() - 1]));
}