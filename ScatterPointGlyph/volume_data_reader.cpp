/*
 Copyright (C) Hongsen Liao
 Institute of Computer Graphics and Computer Aided Design
 School of Software, Tsinghua University
 
 Contact: liaohs082@gmail.com
 
 All rights reserved.
*/

#include "volume_data_reader.h"
#include <vector>
#include <string>
using namespace std;

#include "scatter_volume_dataset.h"
#include "pca_projector.h"

#define BigtoLittle32(A) ((((unsigned int)(A) & 0xff000000) >> 24) | (((unsigned int)(A) & 0x00ff0000) >> 8) | \
             (((unsigned int)(A) & 0x0000ff00) << 8) | (((unsigned int)(A) & 0x000000ff) << 24))

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

VolumeDataReader::VolumeDataReader() {

}

VolumeDataReader::~VolumeDataReader() {

}

ScatterPointDataset* VolumeDataReader::LoadFile(const char* file_name) {
    ScatterVolumeDataset* volume = new ScatterVolumeDataset(500, 500, 100);
    string path_names[] = {"CLOUDf25", "Pf25", "PRECIPf25", "QCLOUDf25", "QGRAUPf25", "QRAINf25", "QSNOWf25", "QVAPORf25", "TCf25", "Uf25", "Vf25", "Wf25"};

    //vector<int> temp_index{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    vector<int> temp_index{0, 1, 2, 7};
    int var_num = temp_index.size();
    for (int i = 0; i < var_num; ++i)
        volume->var_names.push_back(QString(path_names[i].c_str()));

    for (int i = 0; i < 12; ++i)
        path_names[i] = string("./TestData/hurricane/") + path_names[i] + string(".bin");
    float low_range[] = {0.0, -5472, 0, 0, 0, 0, 0, 0, -84, -80, -77, -10};
    float high_range[] = {0.004, 3226, 0.02, 0.004, 0.02, 0.02, 0.002, 0.03, 32, 86, 83, 29};

    int voxel_size = 500 * 500 * 100;
    
	volume->original_point_values.resize(var_num);

    vector<float> temp_data;
    
    temp_data.resize(500 * 500 * 100);
    for (int i = 0; i < var_num; ++i) {
        int index = temp_index[i];
        FILE *pfile = fopen(path_names[index].c_str(), "rb");
        fread(temp_data.data(), sizeof(float) * 500 * 500 * 100, 1, pfile);
        fclose(pfile);

        volume->original_point_values[i].resize(voxel_size);

        float* data_ptr = temp_data.data();
        for (int j = 0; j < temp_data.size(); ++j) {
            temp_data[j] = ReverseFloat(temp_data[j]);
            volume->original_point_values[i][j] = temp_data[j];
        }
    }

    volume->original_point_pos.resize(2);
    volume->original_point_pos[0].resize(voxel_size);
    volume->original_point_pos[1].resize(voxel_size);
    for (int i = 0; i < voxel_size; ++i) {
        volume->original_point_pos[0][i] = volume->original_point_values[2][i];
        volume->original_point_pos[1][i] = volume->original_point_values[3][i];
    }

    vector<vector<float>> value_ranges;
    value_ranges.resize(var_num);
    for (int i = 0; i < var_num; ++i) {
        value_ranges[i].resize(2);
        value_ranges[i][0] = low_range[temp_index[i]];
        value_ranges[i][1] = high_range[temp_index[i]];
    }
    volume->ConstructWithValueRanges(value_ranges);

    int step = 10;
    int index = 0;
    vector<vector<float>> vals;
    vals.resize(var_num);
    for (int z = step / 2; z < 100 / step - 1; ++z)
        for (int y = step / 2; y < 500 / step - 1; ++y)
            for (int x = step / 2; x < 500 / step - 1; ++x) {
                int index = x + y * 500 + z * 500 * 500;
                if (volume->is_valid[index]) {
                    for (int i = 0; i < var_num; ++i)
                        vals[i].push_back(volume->normalized_point_values[i][index]);
                }
            }

    PcaProjector projector;
    projector.Evaluate(vals);
    projector.Apply(volume->normalized_point_values, volume->original_point_pos);

    volume->NormalizePos();

    return volume;
}