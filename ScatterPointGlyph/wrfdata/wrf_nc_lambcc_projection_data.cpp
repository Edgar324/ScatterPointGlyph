#include "wrf_nc_lambcc_projection_data.h"
#include "netcdf.h"
#include <iostream>

WrfNcLambCcProjectionData::WrfNcLambCcProjectionData()
    : data_values_(NULL), map_(NULL) {

}

WrfNcLambCcProjectionData::WrfNcLambCcProjectionData(std::string& file_name)
    : data_values_(NULL), map_(NULL){
    file_name_ = file_name;

    //LoadConvertedData("test.nc");

    LoadData(file_name.c_str());
    //LoadTotalData(file_name.c_str());
}

WrfNcLambCcProjectionData::WrfNcLambCcProjectionData(WrfNcLambCcProjectionData& data){
    this->model_ = data.model_;
    this->element_ = data.element_;
    this->lon_ = data.lon_;
    this->lat_ = data.lat_;
    this->x_dist_ = data.x_dist_;
    this->y_dist_ = data.y_dist_;
    this->standard_parallel_[0] = data.standard_parallel_[0];
    this->standard_parallel_[1] = data.standard_parallel_[1];
    this->longitude_of_central_meridian_ = data.longitude_of_central_meridian_;
    this->latitude_of_projection_origin_ = data.latitude_of_projection_origin_;
    this->false_easting_ = data.false_easting_;
    this->false_northing_ = data.false_northing_;
    this->time_map_ = data.time_map_;
    this->time_len_ = data.time_len_;
    this->y_len_ = data.y_len_;
    this->x_len_ = data.x_len_;
    this->nbnds_len_ = data.nbnds_len_;
    this->map_range_ = data.map_range_;
    this->data_values_ = new float[y_len_ * x_len_ * time_len_];
    memcpy(this->data_values_, data.data_values_, sizeof(float) * y_len_ * x_len_ * time_len_);
}

WrfNcLambCcProjectionData::~WrfNcLambCcProjectionData(){

}

void WrfNcLambCcProjectionData::LoadData(const char* file_name){
    int ncid, dim_id, var_id;
    int status;
    double scale_factor, add_offset;
    nc_type att_type;

    nc_open(file_name, NC_NOWRITE, &ncid);

    nc_inq_dimid(ncid, "time", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &time_len_);

    nc_inq_dimid(ncid, "y", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &y_len_);

    nc_inq_dimid(ncid, "x", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &x_len_);

    nc_inq_dimid(ncid, "nbnds", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &nbnds_len_);

    // load analysis time
    nc_inq_varid(ncid, "time", &var_id);
    std::vector<double> temp_time;
    temp_time.resize(time_len_);
    nc_get_var_double(ncid, var_id, temp_time.data());
    for ( int i = 0; i <temp_time.size(); ++i ) {
        time_map_.insert(std::map<int, int>::value_type((int)temp_time[i], i));
    }

    nc_inq_varid(ncid, "lat", &var_id);
    lat_.resize(y_len_ * x_len_);
    nc_get_var_float(ncid, var_id, lat_.data());

    nc_inq_varid(ncid, "lon", &var_id);
    lon_.resize(y_len_ * x_len_);
    nc_get_var_float(ncid, var_id, lon_.data());

    nc_inq_varid(ncid, "x", &var_id);
    x_dist_.resize(x_len_);
    nc_get_var_float(ncid, var_id, x_dist_.data());

    nc_inq_varid(ncid, "y", &var_id);
    y_dist_.resize(y_len_);
    nc_get_var_float(ncid, var_id, y_dist_.data());

    nc_inq_varid(ncid, "Lambert_Conformal", &var_id);
    nc_get_att(ncid, var_id, "standard_parallel", standard_parallel_);
    nc_get_att(ncid, var_id, "longitude_of_central_meridian", &longitude_of_central_meridian_);
    nc_get_att(ncid, var_id, "latitude_of_projection_origin", &latitude_of_projection_origin_);
    nc_get_att(ncid, var_id, "false_northing", &false_northing_);
    nc_get_att(ncid, var_id, "false_easting", &false_easting_);

    std::vector<double> time_bnds;
    time_bnds.resize(time_len_ * nbnds_len_);
    nc_inq_varid(ncid, "time_bnds", &var_id);
    nc_get_var_double(ncid, var_id, time_bnds.data());

    nc_inq_varid(ncid, "apcp", &var_id);
    nc_get_att(ncid, var_id, "scale_factor", &scale_factor);
    nc_get_att(ncid, var_id, "add_offset", &add_offset);
    size_t data_length = y_len_ * x_len_ * time_len_;    
    data_values_ = new float[data_length];
    if ( data_values_ == NULL ) {
        std::cout <<"Out of memory!" <<std::endl;
        return;
    }

    short* temp_data_values = new short[data_length];
    if ( temp_data_values == NULL ){
        std::cout <<"Out of memory!" <<std::endl;
        delete []data_values_;
        data_values_ = NULL;
        return;
    }
    nc_get_var_short(ncid, var_id, temp_data_values);
    for ( int i = 0; i <data_length; ++i ){
        data_values_[i] = temp_data_values[i] * scale_factor + add_offset;
    }
    delete []temp_data_values;

    nc_close(ncid);

    model_ = WRF_REANALYSIS;
    element_ = WRF_ACCUMULATED_PRECIPITATION;
    map_range_.start_x = x_dist_[0];
    map_range_.end_x = x_dist_[x_dist_.size() - 1];
    map_range_.start_y = y_dist_[0];
    map_range_.end_y = y_dist_[y_dist_.size() - 1];
    map_range_.x_grid_space = x_dist_[1] - x_dist_[0];
    map_range_.y_grid_space = y_dist_[1] - y_dist_[0];
    map_range_.x_grid_number = x_len_;
    map_range_.y_grid_number = y_len_;
}

void WrfNcLambCcProjectionData::LoadConvertedData(const char* file_name){
    int ncid, dim_id, var_id;
    int status;
    double scale_factor, add_offset;
    nc_type att_type;

    nc_open(file_name, NC_NOWRITE, &ncid);

    nc_inq_dimid(ncid, "time", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &time_len_);

    nc_inq_dimid(ncid, "y", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &y_len_);

    nc_inq_dimid(ncid, "x", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &x_len_);

    nc_inq_dimid(ncid, "nbnds", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &nbnds_len_);

    // load analysis time
    nc_inq_varid(ncid, "time", &var_id);
    std::vector<int> temp_time;
    temp_time.resize(time_len_);
    nc_get_var_int(ncid, var_id, temp_time.data());
    for ( int i = 0; i <temp_time.size(); ++i ) {
        time_map_.insert(std::map<int, int>::value_type(temp_time[i], i));
    }

    nc_inq_varid(ncid, "lat", &var_id);
    lat_.resize(y_len_ * x_len_);
    nc_get_var_float(ncid, var_id, lat_.data());

    nc_inq_varid(ncid, "lon", &var_id);
    lon_.resize(y_len_ * x_len_);
    nc_get_var_float(ncid, var_id, lon_.data());

    nc_inq_varid(ncid, "x", &var_id);
    x_dist_.resize(x_len_);
    nc_get_var_float(ncid, var_id, x_dist_.data());

    nc_inq_varid(ncid, "y", &var_id);
    y_dist_.resize(y_len_);
    nc_get_var_float(ncid, var_id, y_dist_.data());

    nc_inq_varid(ncid, "Lambert_Conformal", &var_id);
    nc_get_att(ncid, var_id, "standard_parallel", standard_parallel_);
    nc_get_att(ncid, var_id, "longitude_of_central_meridian", &longitude_of_central_meridian_);
    nc_get_att(ncid, var_id, "latitude_of_projection_origin", &latitude_of_projection_origin_);
    nc_get_att(ncid, var_id, "false_northing", &false_northing_);
    nc_get_att(ncid, var_id, "false_easting", &false_easting_);

    std::vector<double> time_bnds;
    time_bnds.resize(time_len_ * nbnds_len_);
    nc_inq_varid(ncid, "time_bnds", &var_id);
    nc_get_var_double(ncid, var_id, time_bnds.data());

    nc_inq_varid(ncid, "apcp", &var_id);
    size_t data_length = y_len_ * x_len_ * time_len_;
    data_values_ = new float[data_length];
    if ( data_values_ == NULL ) {
        std::cout <<"Out of memory!" <<std::endl;
        return;
    }
    nc_get_var_float(ncid, var_id, data_values_);

    nc_close(ncid);

    model_ = WRF_REANALYSIS;
    element_ = WRF_ACCUMULATED_PRECIPITATION;
    map_range_.start_x = x_dist_[0];
    map_range_.end_x = x_dist_[x_dist_.size() - 1];
    map_range_.start_y = y_dist_[0];
    map_range_.end_y = y_dist_[y_dist_.size() - 1];
    map_range_.x_grid_space = x_dist_[1] - x_dist_[0];
    map_range_.y_grid_space = y_dist_[1] - y_dist_[0];
    map_range_.x_grid_number = x_len_;
    map_range_.y_grid_number = y_len_;
}

void WrfNcLambCcProjectionData::LoadTotalData(const char* file_name){
    int ncid, dim_id, var_id;
    int status;

    nc_open(file_name, NC_NOWRITE, &ncid);

    status = nc_inq_dimid(ncid, "time", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &time_len_);

    nc_inq_dimid(ncid, "ydist", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &y_len_);

    nc_inq_dimid(ncid, "xdist", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &x_len_);

    // load analysis time
    nc_inq_varid(ncid, "time", &var_id);
    std::vector<int> temp_time;
    temp_time.resize(time_len_);
    nc_get_var_int(ncid, var_id, temp_time.data());
    for ( int i = 0; i <temp_time.size(); ++i ) {
        time_map_.insert(std::map<int, int>::value_type(temp_time[i] - 15769752, i));
    }

    nc_inq_varid(ncid, "lon", &var_id);
    lon_.resize(x_len_ * y_len_);
    nc_get_var_float(ncid, var_id, lon_.data());

    nc_inq_varid(ncid, "lat", &var_id);
    lat_.resize(x_len_ * y_len_);
    nc_get_var_float(ncid, var_id, lat_.data());

    nc_inq_varid(ncid, "xdist", &var_id);
    x_dist_.resize(x_len_);
    nc_get_var_float(ncid, var_id, x_dist_.data());

    nc_inq_varid(ncid, "ydist", &var_id);
    y_dist_.resize(y_len_);
    nc_get_var_float(ncid, var_id, y_dist_.data());

    nc_inq_varid(ncid, "apcp", &var_id);
    data_values_ = new float[y_len_ * x_len_ * time_len_];
    nc_get_var_float(ncid, var_id, data_values_);

    nc_close(ncid);

    model_ = WRF_REANALYSIS;
    element_ = WRF_ACCUMULATED_PRECIPITATION;
    map_range_.start_x = x_dist_[0];
    map_range_.end_x = x_dist_[x_dist_.size() - 1];
    map_range_.start_y = y_dist_[0];
    map_range_.end_y = y_dist_[y_dist_.size() - 1];
    map_range_.x_grid_space = x_dist_[1] - x_dist_[0];
    map_range_.y_grid_space = y_dist_[1] - y_dist_[0];
    map_range_.x_grid_number = x_len_;
    map_range_.y_grid_number = y_len_;
}

void WrfNcLambCcProjectionData::SaveData(const char* file_name){
    int ncid, dim_id, var_id;
    int status;
    double scale_factor, add_offset;
    nc_type att_type;
    int pdim[4];
    int temp_ids[5];
    int pvar[8];

    nc_create(file_name, NC_WRITE, &ncid);

    nc_def_dim(ncid, "time", time_len_, &pdim[0]); 
    nc_def_dim(ncid, "y", y_len_, &pdim[1]);
    nc_def_dim(ncid, "x", x_len_, &pdim[2]);
    nc_def_dim(ncid, "nbnds", nbnds_len_, &pdim[3]);

    temp_ids[0] = pdim[0];
    nc_def_var(ncid, "time", NC_INT, 1, temp_ids, &pvar[0]);
    temp_ids[0] = pdim[1];
    temp_ids[1] = pdim[2];
    nc_def_var(ncid, "lat", NC_FLOAT, 2, temp_ids, &pvar[1]);
    nc_def_var(ncid, "lon", NC_FLOAT, 2, temp_ids, &pvar[2]);
    temp_ids[0] =pdim[2];
    nc_def_var(ncid, "x", NC_FLOAT, 1, temp_ids, &pvar[3]);
    temp_ids[0] =pdim[1];
    nc_def_var(ncid, "y", NC_FLOAT, 1, temp_ids, &pvar[4]);

    nc_def_var(ncid, "Lambert_Conformal", NC_CHAR, 0, temp_ids, &pvar[5]);
    nc_put_att(ncid, pvar[5], "standard_parallel", NC_DOUBLE, 2, standard_parallel_);
    nc_put_att(ncid, pvar[5], "longitude_of_central_meridian", NC_DOUBLE, 1, &longitude_of_central_meridian_);
    nc_put_att(ncid, pvar[5], "latitude_of_projection_origin", NC_DOUBLE, 1, &latitude_of_projection_origin_);
    nc_put_att(ncid, pvar[5], "false_northing", NC_DOUBLE, 1, &false_northing_);
    nc_put_att(ncid, pvar[5], "false_easting", NC_DOUBLE, 1, &false_easting_);

    temp_ids[0] = pdim[0];
    temp_ids[1] = pdim[1];
    temp_ids[2] = pdim[2];
    nc_def_var(ncid, "apcp", NC_FLOAT, 3, temp_ids, &pvar[6]);

    nc_enddef(ncid);

    std::vector<int> temp_time;
    std::map<int, int>::iterator iter = time_map_.begin();
    while ( iter != time_map_.end() ){
        temp_time.push_back(iter->first);
        iter++;
    }
    nc_put_var(ncid, pvar[0], temp_time.data());
    nc_put_var(ncid, pvar[1], lat_.data());
    nc_put_var(ncid, pvar[2], lon_.data());
    nc_put_var(ncid, pvar[3], x_dist_.data());
    nc_put_var(ncid, pvar[4], y_dist_.data());
    nc_put_var(ncid, pvar[6], data_values_);

    nc_close(ncid);
}

bool WrfNcLambCcProjectionData::CheckExisting(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num){
    if ( model != model_ || element != element_ || time_map_.find(time) == time_map_.end() ) return false;
    else return true;
}

bool WrfNcLambCcProjectionData::CheckExisting(WrfModelType model){
    if ( model != model_ ) return false;
    else return true;
}

float* WrfNcLambCcProjectionData::GetData(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num){
    if ( data_values_ == NULL ) return NULL;

    long index = time_map_[time] * x_len_ * y_len_;

    return data_values_ + index;
}

WrfLambCcValueMap* WrfNcLambCcProjectionData::GetMap(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num){
    if ( data_values_ == NULL ) return NULL;

    long index = time_map_[time] * x_len_ * y_len_;
    if ( map_ == NULL ) {
        map_ = new WrfLambCcValueMap;
    }
    map_->model_type = model;
    map_->element_type = element;
    map_->map_range = this->map_range_;
    map_->values = this->data_values_ + index;
    map_->lon = this->lon_;
    map_->lat = this->lat_;
    map_->standard_parallel[0] = this->standard_parallel_[0];
    map_->standard_parallel[1] = this->standard_parallel_[1];
    map_->longitude_of_central_meridian = this->longitude_of_central_meridian_;
    map_->latitude_of_projection_origin = this->latitude_of_projection_origin_;
    map_->false_easting = this->false_easting_;
    map_->false_northing = this->false_northing_;

    return map_;
}

void WrfNcLambCcProjectionData::GetMaps(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<WrfGridValueMap*>& maps){

}

void WrfNcLambCcProjectionData::GetData(int time, WrfModelType model, WrfElementType element, int fhour, std::vector<float*>& data){

}

WrfNcLambCcProjectionData* WrfNcLambCcProjectionData::Convert2Range(float start_lon, float end_lon, float start_lat, float end_lat){
    WrfNcLambCcProjectionData* lambcc_data = new WrfNcLambCcProjectionData;

    int left = -1, bottom = -1, right = -1, top = -1;
    float min_dis = 1e10;
    std::vector<bool> x_exist, y_exist;
    x_exist.resize(this->x_len_, false);
    y_exist.resize(this->y_len_, false);
    for ( int i = 0; i <y_len_; ++i ){
        for ( int j = 0; j <x_len_; ++j ){
            float temp_lon = this->lon_[i * x_len_ + j];
            float temp_lat = this->lat_[i * x_len_ + j];
            if ( (temp_lat - start_lat) * (temp_lat - end_lat) <= 0 && (temp_lon - start_lon) * (temp_lon - end_lon) <= 0 ){
                
            }
        }
    }
    for ( int i = 0; i <x_len_; ++i ){
        if ( left == -1 && x_exist[i] ){
            left = i;
        }
        if ( x_exist[i] ) right = i;
    }
    for ( int i = 0; i <y_len_; ++i ){
        if ( bottom == -1 & y_exist[i] ){
            bottom = i;
        }
        if ( y_exist[i] ) top = i;
    }

    lambcc_data->time_len_ = this->time_len_;
    lambcc_data->y_len_ = top - bottom + 1;
    lambcc_data->x_len_ = right - left + 1;
    lambcc_data->nbnds_len_ = this->nbnds_len_;

    lambcc_data->model_ = this->model_;
    lambcc_data->element_ = this->element_;
    lambcc_data->x_dist_.resize(right - left + 1);
    for ( int i = 0; i <lambcc_data->x_len_; ++i ) lambcc_data->x_dist_[i] = this->x_dist_[left + i];
    lambcc_data->y_dist_.resize(top - bottom + 1);
    for ( int i = 0; i <lambcc_data->y_len_; ++i ) lambcc_data->y_dist_[i] = this->y_dist_[bottom + i];
    lambcc_data->lon_.resize(lambcc_data->x_len_ * lambcc_data->y_len_);
    lambcc_data->lat_.resize(lambcc_data->x_len_ * lambcc_data->y_len_);
    for ( int i = 0; i <lambcc_data->y_len_; ++i )
        for ( int j = 0; j <lambcc_data->x_len_; ++j ){
            lambcc_data->lon_[i * lambcc_data->x_len_ + j] = this->lon_[(i + bottom) * this->x_len_ + j + left];
            lambcc_data->lat_[i * lambcc_data->x_len_ + j] = this->lat_[(i + bottom) * this->x_len_ + j + left];
        }
    lambcc_data->standard_parallel_[0] = this->standard_parallel_[0];
    lambcc_data->standard_parallel_[1] = this->standard_parallel_[1];
    lambcc_data->longitude_of_central_meridian_ = this->longitude_of_central_meridian_;
    lambcc_data->latitude_of_projection_origin_ = this->latitude_of_projection_origin_;
    lambcc_data->false_easting_ = this->false_easting_;
    lambcc_data->false_northing_ = this->false_northing_;
    lambcc_data->time_map_ = this->time_map_;
    lambcc_data->map_range_ = this->map_range_;
    lambcc_data->data_values_ = new float[lambcc_data->time_len_ * lambcc_data->x_len_ * lambcc_data->y_len_];
    for ( int t = 0; t <time_len_; ++t )
        for ( int y = 0; y <lambcc_data->y_len_; ++y )
            for ( int x = 0; x <lambcc_data->x_len_; ++x ){
                int pre_index = t * y_len_ * x_len_ + (y + bottom ) * x_len_ + x + left;
                int current_index = t * lambcc_data->y_len_ * lambcc_data->x_len_ + y * lambcc_data->x_len_ + x;
                lambcc_data->data_values_[current_index] = this->data_values_[pre_index];
            }

    return lambcc_data;
}