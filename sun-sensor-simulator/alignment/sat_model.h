#include "../math/matlib.h"
class Sat_Model {
    public: 
        Matrix<double> mag_misalignment;
        Matrix<double> gyro_misalignments;
        std::vector<Matrix<double>> ss_misaligments;
        Vector<double> sun_vector;
        Vector<double> mag_vector;

        Matrix<double> dcm_SB_mag;
        Matrix<double> dcm_SB_gyro;
        std::vector<Matrix<double>> dcm_SB_ss;

        double angle_z;
        double angle_y;

        Quaternion<double> attitude;

        Vector<double> gyro_measurment;

    public: 
        Sat_Model(Matrix<double> mag_misalignment, Matrix<double> gyro_misalignments, std::vector<Matrix<double>> ss_misaligments, Vector<double> sun_vector, Vector<double> mag_vector);
        std::vector<Vector<double>> get_sun_vectors();
        Vector<double> get_mag_vector();
        Vector<double> get_gyro_vector();
        void rotate(double angle, Vector<double> axis);
};
