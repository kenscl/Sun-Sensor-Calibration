#include "sun_sensor.h"

Sun_Sensor::Sun_Sensor(){}

Vector_3D Sun_Sensor::generate_response_v1(double angle_x, double angle_y) {
    double response_x, response_y;
    Vector_3D result;
    response_x = 0.3 * sin (3./2 * angle_x);
    response_y = 0.3 * sin (M_PI / (2 * 60) * angle_y);
    result.x = response_x;
    result.y = response_y;

    return result;
}

Vector_3D Sun_Sensor::generate_response_v2(double angle_x, double angle_y) {
    double response_x, response_y;
    Vector_3D result;
    double sign_x = angle_x / abs(angle_x);
    response_x = 
        0.34 * angle_x -
        sign_x * 0.044 * pow(abs(angle_x), 0.8) * sin(7.4 * pow(abs(angle_x), 0.66)) 
        + 0.001 * sin (60 * angle_x);
        response_y = 0.315 * sin (M_PI / (2 * 60) * angle_y) 
        + 0.002 * angle_y * sin(M_PI / 50 * angle_y) * sin(M_PI / 50 * angle_y);
    result.x = response_x;
    result.y = response_y;

    return result;
}

Vector_3D Sun_Sensor::generate_response_v3(double angle_x, double angle_y) {
    double response_x, response_y;
    Vector_3D result;
    response_x = 0.3 * sin (3./2 * angle_x) - 0.015 * sin(10 * angle_x) + 0.001 * sin (60 * angle_x);
    response_y = 0.315 * sin (M_PI / (2 * 60) * angle_y) 
            + 0.002 * angle_y * sin(M_PI / 50 * angle_y) * sin(M_PI / 50 * angle_y)
            + 0.01 * sin (M_PI / (20) * angle_y) * sin (M_PI / (20) * angle_y) * sin (M_PI / (20) * angle_y);
            result.x = response_x;
            result.y = response_y;

    return result;
}

