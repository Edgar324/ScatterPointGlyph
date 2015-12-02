#include "wrf_nc_square_projection_data.h"
#include <QtCore/QDateTime>
#include "netcdf.h"
#include <iostream>
#include <fstream>

WrfNcSquareProjectionData::WrfNcSquareProjectionData(){

}

WrfNcSquareProjectionData::WrfNcSquareProjectionData(std::string& file_name, WrfModelType model, WrfElementType element)
    : data_values_(NULL), element_(element), model_(model){
    file_name_ = file_name;

    LoadData(file_name_.c_str());
}

WrfNcSquareProjectionData::WrfNcSquareProjectionData(std::string& file_name, WrfModelType model, WrfElementType element, int year)
    : data_values_(NULL), element_(element), model_(model){
    time_len_ = 365;
    if ( (year % 4 == 0 && year % 100 != 0) || (year % 100 == 0 && year % 400 == 0) ) time_len_ = 366;
    QDateTime temp_time = QDateTime(QDate(year, 1, 1), QTime(0, 0));
    int time = temp_time.msecsTo(QDateTime(QDate(1800, 1, 1), QTime(0, 0, 0))) / 3600000 * -1;
    time -= 24;
    for ( int i = 0; i < time_len_; ++i ){
        time_map_.insert(std::map< int, int >::value_type(time, i));
        time += 24;
    }

    fhour_len_ = 1;
    fhour_map_.insert(std::map< int, int >::value_type(24, 0));

    lon_len_ = 280;
    lat_len_ = 180;
    map_range_.start_x = 70.125;
    map_range_.x_grid_space = 0.25;
    map_range_.x_grid_number = 280;
    map_range_.end_x = map_range_.start_x + map_range_.x_grid_space * (map_range_.x_grid_number - 1);
    map_range_.start_y = 15.125;
    map_range_.y_grid_space = 0.25;
    map_range_.y_grid_number = 180;
    map_range_.end_y = map_range_.start_y + map_range_.y_grid_space * (map_range_.y_grid_number - 1);

    lon_.resize(lon_len_);
    for ( int i = 0; i < lon_len_; ++i )
        lon_[i] = map_range_.start_x + i * map_range_.x_grid_space;
    lat_.resize(lat_len_);
    for ( int i = 0; i < lat_len_; ++i )
        lat_[i] = map_range_.start_y + i * map_range_.y_grid_space;

    ens_len_ = 1;

    data_values_ = new float[time_len_ * lon_len_ * lat_len_];

    std::vector< float > temp_values;
    temp_values.resize(time_len_ * lon_len_ * lat_len_ * 2);
    std::ifstream input_stream(file_name.c_str(), std::ios::binary);
    if ( input_stream.good() ){
        input_stream.read((char*)temp_values.data(), time_len_ * lon_len_ * lat_len_ * sizeof(float) * 2);
    }
    for ( int i = 0; i < time_len_; ++i ){
        memcpy(data_values_ + lon_len_ * lat_len_ * i, (char*)(&temp_values[i * lon_len_ * lat_len_ * 2]), lon_len_ * lat_len_ * sizeof(float));
    }
    input_stream.close();
}

WrfNcSquareProjectionData::~WrfNcSquareProjectionData(){
    if ( data_values_ != NULL ) delete []data_values_;
}

void WrfNcSquareProjectionData::LoadData(const char* file_name){
    int ncid, dim_id, var_id;
    int status;

    nc_open(file_name, NC_NOWRITE, &ncid);

    nc_inq_dimid(ncid, "time", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &time_len_);

    nc_inq_dimid(ncid, "lat", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &lat_len_);

    nc_inq_dimid(ncid, "lon", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &lon_len_);

    status = nc_inq_dimid(ncid, "ens", &dim_id);
    if ( status == NC_NOERR ){
        nc_inq_dimlen(ncid, dim_id, &ens_len_);
    } else {
        ens_len_ = 1;
    }
    

    nc_inq_dimid(ncid, "fhour", &dim_id);
    nc_inq_dimlen(ncid, dim_id, &fhour_len_);

    // load analysis time
    nc_inq_varid(ncid, "time", &var_id);
    std::vector< double > temp_time;
    temp_time.resize(time_len_ * fhour_len_);
    nc_get_var_double(ncid, var_id, temp_time.data());
    for ( int i = 0; i < temp_time.size(); ++i ) {
        time_map_.insert(std::map< int, int >::value_type((int)temp_time[i], i));
    }

    nc_inq_varid(ncid, "intValidTime", &var_id);
    std::vector< int > intvalidtime;
    intvalidtime.resize(time_len_ * fhour_len_);
    nc_get_var_int(ncid, var_id, intvalidtime.data());

    nc_inq_varid(ncid, "intTime", &var_id);
    std::vector< int > inttime;
    inttime.resize(time_len_ * fhour_len_);
    nc_get_var_int(ncid, var_id, inttime.data());

    nc_inq_varid(ncid, "lat", &var_id);
    lat_.resize(lat_len_);
    nc_get_var_float(ncid, var_id, lat_.data());

    nc_inq_varid(ncid, "lon", &var_id);
    lon_.resize(lon_len_);
    nc_get_var_float(ncid, var_id, lon_.data());
    

    nc_inq_varid(ncid, "fhour", &var_id);
    std::vector< int > temp_fhour;
    temp_fhour.resize(fhour_len_);
    nc_get_var_int(ncid, var_id, temp_fhour.data());
    for ( int i = 0; i < temp_fhour.size(); ++i )
        fhour_map_.insert(std::map< int, int >::value_type(temp_fhour[i], i));

    switch ( element_ ){
    case WRF_ACCUMULATED_PRECIPITATION:
        nc_inq_varid(ncid, "Total_precipitation", &var_id);
        break;
	case WRF_PRECIPITABLE_WATER:
		nc_inq_varid(ncid, "Precipitable_water", &var_id);
		break;
    case WRF_T2M:
        nc_inq_varid(ncid, "Temperature_height_above_ground", &var_id);
        break;
	case WRF_PRESSURE_SURFACE:
		nc_inq_varid(ncid, "Pressure_surface", &var_id);
		break;
	case WRF_MSLP:
		nc_inq_varid(ncid, "Pressure", &var_id);
		break;
    default:
        return;
    }
    
    data_values_ = new float[time_len_ * ens_len_ * fhour_len_ * lat_len_ * lon_len_];
    if ( data_values_ == NULL ){
        std::cout << "Out of memory" << std::endl;
        return;
    }
    nc_get_var_float(ncid, var_id, data_values_);

    nc_close(ncid);

    map_range_.start_x = lon_[0];
    map_range_.end_x = lon_[lon_.size() - 1];
    map_range_.start_y = lat_[0];
    map_range_.end_y = lat_[lat_.size() - 1];
    map_range_.x_grid_space = lon_[1] - lon_[0];
    map_range_.y_grid_space = lat_[1] - lat_[0];
    map_range_.x_grid_number = lon_len_;
    map_range_.y_grid_number = lat_len_;
}

bool WrfNcSquareProjectionData::CheckExisting(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num){
    if ( element != element_ || model != model_ || fhour_map_.find(fhour) == fhour_map_.end() 
        || ens_num >= ens_len_ || time_map_.find(time) == time_map_.end() ) 
        return false;
    
    return true;
}

bool WrfNcSquareProjectionData::CheckExisting(WrfModelType model){
    if ( model == model_ ) return true;
    else return false;
}
float* WrfNcSquareProjectionData::GetData(int time, WrfModelType model, WrfElementType element, int fhour, int ens_num){
    if ( data_values_ == NULL ) return NULL;

    long index = time_map_[time];
    index *= lat_len_ * lon_len_ * fhour_len_ * ens_len_;
    index += ens_num * lat_len_ * lon_len_ * fhour_len_;
    index += lat_len_ * lon_len_ * fhour_map_[fhour];

    return data_values_ + index;
}

WrfGridValueMap* WrfNcSquareProjectionData::GetMap(int time, WrfModelType model, WrfElementType element, int ens_num, int fhour){
    if ( data_values_ == NULL ) return NULL;

    WrfGridValueMap* temp_map = new WrfGridValueMap;

    long index = time_map_[time];
    index *= lat_len_ * lon_len_ * fhour_len_ * ens_len_;
    index += ens_num * lat_len_ * lon_len_ * fhour_len_;
    index += lat_len_ * lon_len_ * fhour_map_[fhour];

    temp_map->model_type = model;
    temp_map->element_type = element;
    temp_map->map_range = this->map_range_;
    temp_map->values = data_values_ + index;

    return temp_map;
}

void WrfNcSquareProjectionData::GetMaps(int time, WrfModelType model, WrfElementType element, int fhour, std::vector< WrfGridValueMap* >& maps){
    if ( data_values_ == NULL ) return;

    maps.clear();
    for ( int ens_num = 0; ens_num < ens_len_; ++ens_num ){
        WrfGridValueMap* temp_map = new WrfGridValueMap;

        long index = time_map_[time];
        index *= lat_len_ * lon_len_ * fhour_len_ * ens_len_;
        index += ens_num * lat_len_ * lon_len_ * fhour_len_;
        index += lat_len_ * lon_len_ * fhour_map_[fhour];

        temp_map->model_type = model;
        temp_map->element_type = element;
        temp_map->map_range = this->map_range_;
        temp_map->values = data_values_ + index;

        maps.push_back(temp_map);
    }
}

void WrfNcSquareProjectionData::GetData(int time, WrfModelType model, WrfElementType element, int fhour, std::vector< float* >& data){
    if ( data_values_ == NULL ) return;

    data.clear();
    for ( int ens_num = 0; ens_num < ens_len_; ++ens_num ){
        long index = time_map_[time];
        index *= lat_len_ * lon_len_ * fhour_len_ * ens_len_;
        index += ens_num * lat_len_ * lon_len_ * fhour_len_;
        index += lat_len_ * lon_len_ * fhour_map_[fhour];

        data.push_back(data_values_ + index);
    }
}