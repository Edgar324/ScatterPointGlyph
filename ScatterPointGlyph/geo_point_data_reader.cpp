/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "geo_point_data_reader.h"
#include <string>
#include <fstream>
using namespace std;
#include <QString>
#include "scatter_point_dataset.h"

GeoPointDataReader::GeoPointDataReader() {

}

GeoPointDataReader::~GeoPointDataReader() {

}

ScatterPointDataset* GeoPointDataReader::LoadFile(const char* file_name) {
    ScatterPointDataset* dataset = new ScatterPointDataset;

    std::ifstream input_file(file_name);
	char char_str[1000];
	input_file.getline(char_str, 1000);
	QString value_str = QString::fromLocal8Bit(char_str);
	QStringList value_list = value_str.split(' ');

	int record_num, var_num;
	record_num = value_list.at(0).toInt();
	var_num = value_list.at(1).toInt();

	input_file.getline(char_str, 1000);
	value_str = QString::fromLocal8Bit(char_str);
	value_list = value_str.split(' ');
	for (int i = 0; i < value_list.size(); ++i) 
        dataset->var_names.push_back(value_list.at(i));

	dataset->original_point_pos.resize(2);
    dataset->original_point_pos[0].resize(record_num);
    dataset->original_point_pos[1].resize(record_num);
	dataset->original_point_values.resize(var_num);
    for (int i = 0; i < var_num; ++i)
        dataset->original_point_values[i].resize(record_num);

	for (int i = 0; i < record_num; ++i) {
		input_file.getline(char_str, 1000);
		value_str = QString::fromLocal8Bit(char_str);
		value_list = value_str.split(',');

		for (int j = 0; j < var_num; ++j)
			dataset->original_point_values[j][i] = value_list.at(2 + j).toFloat();

		dataset->original_point_pos[0][i] = value_list.at(0).toFloat();
		dataset->original_point_pos[1][i] = value_list.at(1).toFloat();
	}
	input_file.close();

    dataset->DirectConstruct();

    return dataset;
}