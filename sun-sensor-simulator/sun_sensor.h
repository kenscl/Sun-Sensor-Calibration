#include "math/matlib.h"
#ifndef __SUN_SENSOR 
#define __SUN_SENSOR 
class Sun_Sensor {
public:
    Sun_Sensor();

    Vector_3D generate_response_v1(double angle_x, double angle_y);
    Vector_3D generate_response_v2(double angle_x, double angle_y);
    Vector_3D generate_response_v3(double angle_x, double angle_y);

};

#endif
