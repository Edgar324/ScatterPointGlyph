#ifndef WRF_UTILITY_H_
#define WRF_UTILITY_H_

#include "wrf_data_common.h"

extern "C" {
    char* enum_model_to_string(WrfModelType type);
    char* enum_element_to_string(WrfElementType type);
}
#endif