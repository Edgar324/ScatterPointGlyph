#include "wrf_utility.h"

char* enum_model_to_string(WrfModelType type){
    switch ( type ){
    case WRF_UNKNOWN_MODEL:
        return "Unknown Model";
    case WRF_NCEP_ENSEMBLES:
        return "Ensemble";
    case WRF_NCEP_ENSEMBLE_MEAN:
        return "Ensemble Mean";
    case WRF_REANALYSIS:
        return "Observation";
    case WRF_T639:
        return "T639";
    case WRF_EC_FINE:
        return "Ec Fine";
    case WRF_JAPAN_GSM:
        return "Japan Gsm";
    case WRF_NCEP:
        return "NCEP";
    case WRF_ELEMENT_UNCERTAINTY:
        return "Uncertainty";
    case WRF_ELEMENT_PROBABILISTIC:
        return "Probabilistic";
    default:
        return "";
    }
};

char* enum_element_to_string(WrfElementType type){
    switch ( type ){
    case WRF_UNKNOWN_ELEMENT:
        return "Unknown Element";
    case WRF_ACCUMULATED_PRECIPITATION:
        return "APCP";
    case WRF_PRECIPITABLE_WATER:
        return "PrecW";
    case WRF_CAPE:
        return "CAPE";
    case WRF_CIN:
        return "CIN";
    case WRF_T2M:
        return "T2M";
    case WRF_MSLP:
        return "MSLP";
	case WRF_PRESSURE_SURFACE:
		return "Surface Pressure";
    case WRF_TEMPERATURE_500HPA:
        return "temper500";
    case WRF_HEIGHT_500HPA:
        return "height500";
    case WRF_RELATIVE_HUMIDITY_500HPA:
        return "rh500";
    case WRF_UNCERTAINTY:
        return "Uncertainty";
    case WRF_PROBABILISTIC:
        return "Probabilistic";
    default:
        return "";
    }
};