#pragma once
#include <../math/matlib.h>
#include <string>
#include <vector>
class Vasco {
public:
    // filter
    Quaternion<double> q_i_c, q_bk ;
    Vector<double> x, x_hat;
    Vector<double> y;
    Vector<double> w_x;
    Vector<double> w_y;

    Matrix<double> F;
    Matrix<double> A;
    Matrix<double> H;
    Matrix<double> Q;
    Matrix<double> R;
    Matrix<double> P, P_bar;
    Matrix<double> K;

    double convergence_threshold;

    //correction variables
    
    Vector<double> bias; 
    Matrix<double> mag_misalignment;
    Vector<double> mag_bias;
    std::vector<Matrix<double>> sun_sensor_misalignment;
    std::vector<Vector<double>> sun_sensor_bias;
    std::vector<Vector<double>> cov_matrix_diag;

    // record keeping
    
    std::vector<Quaternion<double>> attitude_estimates;
    std::vector<Vector<double>> residuals;
    std::vector<double> elapsed_time;
    std::vector<Quaternion<double>> mag_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_xp_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_xm_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_yp_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_ym_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_zp_misalignment_estimates;
    std::vector<Quaternion<double>> sun_sensor_zm_misalignment_estimates;
public:
    Vasco();
    void determine_misalignments(Vector<double> sun_refrence, Vector<double> mag_refrence, std::string sun_measured, std::string mag_measured, std::string gyro_measured);

    Quaternion<double> triad(Vector<double> sun_refrence, Vector<double> sun_measured, Vector<double> mag_refrence, Vector<double> mag_measured);

    // filter stuff
    void propagate(Vector<double> w, double dt);
    void predict(double dt, Vector<double> w, Vector<double> mag_model, Vector<double> sun_model, int sun_sensor_id);
    void update(Vector<double> mag_measured, Vector<double> sun_measured, int sun_sensor_id, Vector<double> mag_refrence, Vector<double> sun_refrence );

    void set_attitude_error_variance(double variance);
    void set_gyro_error_variance(double variance);
    void set_mag_error_variance(double variance);
    void set_sun_error_variance(double variance);

    void set_mag_measurement_variance(double variance);
    void set_sun_measurement_variance(double variance);

    void set_attitude_process_variance(double variance);
    void set_gyro_process_variance(double variance);
    void set_mag_process_variance(double variance);
    void set_sun_process_variance(double variance);

    bool has_converged();
    void feed_forward();
    void full_feedback(double elapsed_time);
    void partial_feedback();

    void store_data(double time);

    Vector<double> correct_gyro(Vector<double> w);
    Vector<double> correct_mag(Vector<double> mag);
    Vector<double> correct_sun(Vector<double> sun, int id);

    void correct_matricies();

    void run(std::vector<double> dt, std::vector<Vector<double>> w, std::vector<Vector<double>> mag_model, std::vector<Vector<double>> mag_measured, std::vector<Vector<double>> sun_model, std::vector<Vector<double>> sun_measured, std::vector<int> sun_sensor_id, double threashold,std::vector<Quaternion<double>> attitude); 
    void print(double time);

    Matrix<double> b_to_c();
};
