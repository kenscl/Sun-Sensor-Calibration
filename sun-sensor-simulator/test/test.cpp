#include "test.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <cstdio>
#include <vector>

void split_vector_odd_even(const std::vector<double> init, std::vector<double>& even, std::vector<double>& odd) {
    for (uint i = 0; i < init.size(); ++i) {
        if (i % 2 == 0) {
            even.push_back(init.at(i));
        } else {
            odd.push_back(init.at(i));
        }
    }
}


void polyfit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;
    int err_min_index_x = -1;
    int err_min_index_y = -1;
    double err_min_x = MAXFLOAT;
    double err_min_y = MAXFLOAT;


    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);
    for (uint i = 0; i < 17; ++i) {
        // fit data
        Polynomial_Fit polyfit_x(i, gt_x_even, measurement_x_even);
        meas_x = polyfit_x.calc(measurement_x_odd);
        convert_to_deg(meas_x);


        Polynomial_Fit polyfit_y(i, gt_y_even, measurement_y_even);
        meas_y = polyfit_y.calc(measurement_y_odd);
        convert_to_deg(meas_y);

        // error
        double err_x = rmse(gt_x_odd, meas_x);
        double err_y = rmse(gt_y_odd, meas_y);

        index.push_back(i);
        error_x.push_back(err_x);
        error_y.push_back(err_y);

        error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
        error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

        if (err_x < err_min_x) {
            err_min_x = err_x;
            err_min_index_x = i;
        }

        if (err_y < err_min_y) {
            err_min_y = err_y;
            err_min_index_y = i;
        }

        // size 
        size.push_back(sizeof(std::vector<double>) + sizeof(double) * polyfit_x.weights.size());

        // time
        auto start = timer::high_resolution_clock::now();
        for (int j = 0; j < 1000; ++j) {
            polyfit_x.calc(measurement_y);
        }

        auto end = timer::high_resolution_clock::now();
        time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);


        printf("Degree %d, error x %f, error y %f, err45 x %f, err45 y %f, size %f, time %f\n", i, err_x, err_y, error_45_x.back(), error_45_y.back(), size.back(), time.back());

    }
    printf("X: minimum error at degree %d with error %f\n", err_min_index_x, err_min_x);
    printf("Y: minimum error at degree %d with error %f\n", err_min_index_y, err_min_y);

    //error plot
    //plt::figure_size(800, 600);
    plt::title("POLYNOMIAL FIT ERROR BASED ON DEGREE");
    plt::xlabel("DEGREE OF POLYNOMIAL FIT");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE", index, error_x);
    plt::named_plot("Y-RMSE", index, error_y);
    plt::legend();
    std::string filename = "../plots/Polyfit_Count.eps";
    plt::save(filename);

    print_vector(error_45_x);
    print_vector(error_45_y);
    plt::cla();
    plt::clf();
    plt::title("POLYNOMIAL FIT ERROR BASED ON DEGREE (INNER 45 DEG)");
    plt::xlabel("DEGREE OF POLYNOMIAL FIT");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE 45", index, error_45_x);
    plt::named_plot("Y-RMSE 45", index, error_45_y);
    plt::legend();
    filename = "../plots/Polyfit_45.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("POLYNOMIAL FIT TIME BASED ON DEGREE");
    plt::xlabel("DEGREE OF POLYNOMIAL FIT");
    plt::ylabel("TIME [µs]");
    plt::grid(true);
    plt::plot(index, time);
    filename = "../plots/Polyfit_time.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();
    plt::title("POLYNOMIAL FIT SIZE BASED ON DEGREE");
    plt::xlabel("DEGREE OF POLYNOMIAL FIT");
    plt::ylabel("SIZE IN BYTES");
    plt::grid(true);
    plt::plot(index, size);
    filename = "../plots/Polyfit_size.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();

    // fit data
    Polynomial_Fit polyfit_x(7, gt_x_even, measurement_x_even);
    meas_x = polyfit_x.calc(measurement_x_odd);
    convert_to_deg(meas_x);

    Polynomial_Fit polyfit_y(7, gt_y_even, measurement_y_even);
    meas_y = polyfit_y.calc(measurement_y_odd);
    convert_to_deg(meas_y);

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
    }
    for (uint i = 0; i < meas_y.size(); ++i) {
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("POLYNOMIAL FIT");
    plt::grid(true);
    plt::named_plot("polynomial fit X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("polynomial fit Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    filename = "../plots/Polyfit.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();
}

void sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;
    int err_min_index_x = -1;
    int err_min_index_y = -1;
    double err_min_x = MAXFLOAT;
    double err_min_y = MAXFLOAT;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);
    for (uint i = 3; i < 5; ++i) {
        Sin_Fit fit_x(i, 5000 ,1e-9, gt_x, measurement_x);
        meas_x = fit_x.calc(measurement_x_odd, 1e-3);
        convert_to_deg(meas_x);

        Sin_Fit fit_y(i, 5000 ,1e-9, gt_y, measurement_y);
        meas_y = fit_y.calc(measurement_y_odd, 1e-3);
        convert_to_deg(meas_y);

        // error
        double err_x = rmse(gt_x_odd, meas_x);
        double err_y = rmse(gt_y_odd, meas_y);

        index.push_back(i);
        error_x.push_back(err_x);
        error_y.push_back(err_y);
        error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
        error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

        if (err_x < err_min_x) {
            err_min_x = err_x;
            err_min_index_x = i;
        }

        if (err_y < err_min_y) {
            err_min_y = err_y;
            err_min_index_y = i;
        }

        // size 
        size.push_back(sizeof(Vector<double>) + sizeof(double) * fit_x.parameters.size);

        // time
        auto start = timer::high_resolution_clock::now();
        for (int j = 0; j < 1000; ++j) {
            fit_x.calc(measurement_y, 1e-3);
        }

        auto end = timer::high_resolution_clock::now();
        time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

        //printf("Degree %d, error %f, size %f, time %f\n", i, err, size.back(), time.back());
        printf("Degree %d, error x %f, error y %f, err45 x %f, err45 y %f, size %f, time %f\n", i, err_x, err_y, error_45_x.back(), error_45_y.back(), size.back(), time.back());

    }
    printf("X: minimum error at degree %d with error %f\n", err_min_index_x, err_min_x);
    printf("Y: minimum error at degree %d with error %f\n", err_min_index_y, err_min_y);

    //error plot
    //plt::figure_size(800, 600);
    plt::title("SIN FIT ERROR BASED ON DEGREE");
    plt::xlabel("DEGREE OF SIN FIT");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE", index, error_x);
    plt::named_plot("Y-RMSE", index, error_y);
    plt::legend();
    std::string filename = "../plots/Sin_fit_Count.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("SIN FIT ERROR BASED ON DEGREE (INNER 45 DEG)");
    plt::xlabel("DEGREE OF SIN FIT");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE 45", index, error_45_x);
    plt::named_plot("Y-RMSE 45", index, error_45_y);
    plt::legend();
    filename = "../plots/Sin_fit_45.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("SIN FIT TIME BASED ON DEGREE");
    plt::xlabel("DEGREE OF SIN FIT");
    plt::ylabel("TIME [µs]");
    plt::grid(true);
    plt::plot(index, time);
    filename = "../plots/Sin_fit_time.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();
    plt::title("SIN FIT SIZE BASED ON DEGREE");
    plt::xlabel("DEGREE OF SIN FIT");
    plt::ylabel("SIZE IN BYTES");
    plt::grid(true);
    plt::plot(index, size);
    filename = "../plots/Sin_fit_size.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();


    Sin_Fit fit_x(3, 5000 ,1e-6, gt_x, measurement_x);
    meas_x = fit_x.calc(measurement_x_odd, 1e-3);
    convert_to_deg(meas_x);

    Sin_Fit fit_y(3, 5000 ,1e-6, gt_y, measurement_y);
    meas_y = fit_y.calc(measurement_y_odd, 1e-3);
    convert_to_deg(meas_y);

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
    }
    for (uint i = 0; i < meas_y.size(); ++i) {
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("SIN FIT");
    plt::grid(true);
    plt::named_plot("Sin fit X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("Sin fit Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    filename = "../plots/Sin_fit.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();

}


void e_sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);
    Extended_sin_fit fit_x(5000 ,1e-6, gt_x_even, measurement_x_even);
    Extended_sin_fit fit_y(5000 ,1e-6, gt_y_even, measurement_y_even);

    meas_x = fit_x.calc(measurement_x_odd, 1e-3);
    convert_to_deg(meas_x);
    double err_x = rmse(gt_x_odd, meas_x);


    meas_y = fit_y.calc(measurement_y_odd, 1e-3);
    convert_to_deg(meas_y);
    double err_y = rmse(gt_y_odd, meas_y);

    error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
    error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

    size.push_back(sizeof(Vector<double>) + sizeof(double) * fit_x.parameters.size);

    auto start = timer::high_resolution_clock::now();
    for (int j = 0; j < 1000; ++j) {
        fit_x.calc(measurement_y, 1e-3);
    }

    auto end = timer::high_resolution_clock::now();
    time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

    printf("X: error %f, size %f, time %f\n", err_x, size.back(), time.back());
    printf("Y: error %f, size %f, time %f\n", err_y, size.back(), time.back());
    printf("X: error 45 %f\n", error_45_x.back());
    printf("Y: error 45 %f\n", error_45_y.back());

    plt::cla();
    plt::clf();

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("MODEL 2");
    plt::grid(true);
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    plt::named_plot("fit X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("fit Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    std::string filename = "../plots/e_Sin_fit.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
}

void alt_sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);

    std::vector<double> test;
    for (float i = -180 * D2R; i < 180 * D2R; i += 0.05) {
        test.push_back(i);
    }

    Alt_sin_fit fit_x(5000 ,1e-9, gt_x, measurement_x);
    // print the parameters
    for (uint i = 0; i < fit_x.parameters.size; ++i) {
        printf("parameter %d: %f\n", i, fit_x.parameters.data.at(i));
    }
    
    Alt_sin_fit fit_y(5000 ,1e-9, gt_y, measurement_y);
    meas_x = fit_x.calc(test, 1e-3);
    meas_y = fit_y.calc(measurement_y_odd, 1e-3);
    convert_to_deg(meas_x);
    convert_to_deg(meas_y);
    convert_to_deg(test);
    //double err_x = rmse(gt_x_odd, meas_x);
    double err_y = rmse(gt_y_odd, meas_y);
    error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
    error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

    size.push_back(sizeof(Vector<double>) + sizeof(double) * fit_x.parameters.size);

    auto start = timer::high_resolution_clock::now();
    for (int j = 0; j < 1000; ++j) {
        fit_x.calc(measurement_y, 1e-3);
    }

    auto end = timer::high_resolution_clock::now();
    time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

    //printf("X: error %f, size %f, time %f\n", err_x, size.back(), time.back());
    printf("Y: error %f, size %f, time %f\n", err_y, size.back(), time.back());
    printf("X: error 45 %f\n", error_45_x.back());
    printf("Y: error 45 %f\n", error_45_y.back());

    plt::cla();
    plt::clf();


    for (uint i = 0; i < meas_x.size(); ++i) {
        //meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("MODEL 3");
    plt::grid(true);
    plt::legend();
    plt::named_plot("fit X-Axis", test, meas_x);
    //plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    //plt::named_plot("fit Y-Axis", gt_y_odd, meas_y);
    //plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    plt::legend();
    plt::ylim(-65,65);
    std::string filename = "../plots/alt_Sin_fit.eps";
    plt::save(filename);
    plt::show();



    plt::cla();
    plt::clf();


}

void lut_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;
    int err_min_index_x = -1;
    int err_min_index_y = -1;
    double err_min_x = MAXFLOAT;
    double err_min_y = MAXFLOAT;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);

    //for (uint i = 5; i < 120; i += 1) {
    //    LUT fit_x(gt_x_even, measurement_x_even, i);
    //    meas_x = fit_x.calc(measurement_x_odd);
    //    convert_to_deg(meas_x);

    //    LUT fit_y(gt_y_even, measurement_y_even, i);
    //    meas_y = fit_y.calc(measurement_y_odd);
    //    convert_to_deg(meas_y);

    //    // error
    //    double err_x = rmse(gt_x_odd, meas_x);
    //    double err_y = rmse(gt_y_odd, meas_y);

    //    index.push_back(i);
    //    error_x.push_back(err_x);
    //    error_y.push_back(err_y);
    //    error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
    //    error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

    //    if (err_x < err_min_x) {
    //        err_min_x = err_x;
    //        err_min_index_x = i;
    //    }

    //    if (err_y < err_min_y) {
    //        err_min_y = err_y;
    //        err_min_index_y = i;
    //    }

    //    // size 
    //    size.push_back((sizeof(std::vector<double>) + sizeof(double) * fit_x.parameters.size()) * 2) ;

    //    // time
    //    auto start = timer::high_resolution_clock::now();
    //    for (int j = 0; j < 1000; ++j) {
    //        fit_x.calc(measurement_y);
    //    }

    //    auto end = timer::high_resolution_clock::now();
    //    time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

    //    printf("Degree %d, error x %f, error y %f, err45 x %f, err45 y %f, size %f, time %f\n", i, err_x, err_y, error_45_x.back(), error_45_y.back(), size.back(), time.back());

    //}
    //printf("X: minimum error at degree %d with error %f\n", err_min_index_x, err_min_x);
    //printf("Y: minimum error at degree %d with error %f\n", err_min_index_y, err_min_y);

    ////error plot
    ////plt::figure_size(800, 600);
    //plt::title("LINEAR INTERPOLATION ERROR BASED ON DEGREE");
    //plt::xlabel("DEGREE OF LINEAR INTERPOLATION");
    //plt::ylabel("RMSE (DEG)");
    //plt::named_plot("X-RMSE", index, error_x);
    //plt::named_plot("Y-RMSE", index, error_y);
    //plt::grid(true);
    //plt::legend();
    //plt::ylim(0., 3.);
    std::string filename = "../plots/LERP_Count.eps";
    //plt::save(filename);

    //plt::cla();
    //plt::clf();
    //plt::title("LINEAR INTERPOLATION ERROR BASED ON DEGREE (INNER 45 DEG)");
    //plt::xlabel("DEGREE OF LINEAR INTERPOLATION");
    //plt::ylabel("RMSE (DEG)");
    //plt::named_plot("X-RMSE 45", index, error_45_x);
    //plt::named_plot("Y-RMSE 45", index, error_45_y);
    //plt::grid(true);
    //plt::legend();
    //plt::ylim(0., 0.05);
    //filename = "../plots/LERP_45.eps";
    //plt::save(filename);

    //plt::cla();
    //plt::clf();
    //plt::title("LINEAR INTERPOLATION TIME BASED ON DEGREE");
    //plt::xlabel("DEGREE OF LINEAR INTERPOLATION");
    //plt::ylabel("TIME [µs]");
    //plt::grid(true);
    //plt::plot(index, time);
    //filename = "../plots/LERP_time.eps";
    //plt::save(filename);


    //plt::cla();
    //plt::clf();
    //plt::title("LINEAR INTERPOLATION SIZE BASED ON DEGREE");
    //plt::xlabel("DEGREE OF LINEAR INTERPOLATION");
    //plt::ylabel("SIZE IN BYTES");
    //plt::grid(true);
    //plt::plot(index, size);
    //filename = "../plots/LERP_size.eps";
    //plt::save(filename);

    plt::cla();
    plt::clf();

    LUT fit_x(gt_x_even, measurement_x_even, 7);
    meas_x = fit_x.calc(measurement_x_odd);
    convert_to_deg(meas_x);

    LUT fit_y(gt_y_even, measurement_y_even, 7);
    meas_y = fit_y.calc(measurement_y_odd);
    convert_to_deg(meas_y);

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
    }
    for (uint i = 0; i < meas_y.size(); ++i) {
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("LINEAR INTERPOLATION");
    plt::grid(true);
    plt::named_plot("Linear interpolation X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("Linear interpolation Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    filename = "../plots/LERP.eps";
    plt::show();
    //plt::save(filename);

    plt::cla();
    plt::clf();
}

void cspline_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {

    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;
    int err_min_index_x = -1;
    int err_min_index_y = -1;
    double err_min_x = MAXFLOAT;
    double err_min_y = MAXFLOAT;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);

    for (uint i = 5; i < 120; i += 1) {
        CSpline fit_x(gt_x, measurement_x, i);
        meas_x = fit_x.calc(measurement_x_odd);
        convert_to_deg(meas_x);

        CSpline fit_y(gt_y_even, measurement_y_even, i);
        meas_y = fit_y.calc(measurement_y_odd);
        convert_to_deg(meas_y);

        // error
        double err_x = rmse(gt_x_odd, meas_x);
        double err_y = rmse(gt_y_odd, meas_y);

        index.push_back(i);
        error_x.push_back(err_x);
        error_y.push_back(err_y);
        error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
        error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

        if (err_x < err_min_x) {
            err_min_x = err_x;
            err_min_index_x = i;
        }

        if (err_y < err_min_y) {
            err_min_y = err_y;
            err_min_index_y = i;
        }

        // size 
        size.push_back(sizeof(std::vector<double>) + sizeof(double) * fit_x.y.size() + sizeof(std::vector<std::vector<double>>) * fit_x.coefficients.size() * 4);

        // time
        auto start = timer::high_resolution_clock::now();
        for (int j = 0; j < 1000; ++j) {
            fit_x.calc(measurement_y);
        }

        auto end = timer::high_resolution_clock::now();
        time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

        //printf("Degree %d, error %f, size %f, time %f\n", i, err, size.back(), time.back());
        printf("Degree %d, error x %f, error y %f, err45 x %f, err45 y %f, size %f, time %f\n", i, err_x, err_y, error_45_x.back(), error_45_y.back(), size.back(), time.back());

    }
    printf("X: minimum error at degree %d with error %f\n", err_min_index_x, err_min_x);
    printf("Y: minimum error at degree %d with error %f\n", err_min_index_y, err_min_y);

    //error plot
    //plt::figure_size(800, 600);
    plt::title("C-SPLINE ERROR BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE", index, error_x);
    plt::named_plot("Y-RMSE", index, error_y);
    plt::legend();
    plt::ylim(0., 5.);
    std::string filename = "../plots/CSpline_Count.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("C-SPLINE ERROR BASED ON DEGREE (INNER 45 DEG)");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE 45", index, error_45_x);
    plt::named_plot("Y-RMSE 45", index, error_45_y);
    plt::legend();
    plt::ylim(0., .1);
    filename = "../plots/CSpline_45.eps";
    plt::save(filename);
    plt::cla();
    plt::clf();

    plt::cla();
    plt::clf();
    plt::title("C-SPLINE TIME BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("TIME [µs]");
    plt::grid(true);
    plt::plot(index, time);
    filename = "../plots/CSpline_time.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();
    plt::title("C-SPLINE SIZE BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("SIZE IN BYTES");
    plt::grid(true);
    plt::plot(index, size);
    filename = "../plots/CSpline_size.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();

    CSpline fit_x(gt_x, measurement_x, 50);
    meas_x = fit_x.calc(measurement_x_odd);
    convert_to_deg(meas_x);

    CSpline fit_y(gt_y_even, measurement_y_even, 50);
    meas_y = fit_y.calc(measurement_y_odd);
    convert_to_deg(meas_y);

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
    }

    for (uint i = 0; i < meas_y.size(); ++i) {
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("C-SPLINE");
    plt::grid(true);
    plt::named_plot("C-Spline X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("C-Spline Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    filename = "../plots/CSpline.eps";
    plt::save(filename);
    plt::show();


    plt::cla();
    plt::clf();
}

void bspline_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y) {
    std::vector<double> index, error_x, error_y, time, size, meas_x, meas_y, meas_err_x, meas_err_y, error_45_x, error_45_y;
    int err_min_index_x = -1;
    int err_min_index_y = -1;
    double err_min_x = MAXFLOAT;
    double err_min_y = MAXFLOAT;

    std::vector<double> gt_x_even, gt_x_odd, measurement_x_even, measurement_x_odd, gt_y_even, gt_y_odd, measurement_y_even, measurement_y_odd;

    split_vector_odd_even(gt_x, gt_x_even, gt_x_odd);
    split_vector_odd_even(measurement_x, measurement_x_even, measurement_x_odd);

    split_vector_odd_even(gt_y, gt_y_even, gt_y_odd);
    split_vector_odd_even(measurement_y, measurement_y_even, measurement_y_odd);

    convert_to_deg(gt_x_odd);
    convert_to_deg(gt_y_odd);

    for (uint i = 75; i < 120; i += 1) {
        BSpline fit_x(gt_x_even, measurement_x_even, i);
        meas_x = fit_x.calc(measurement_x_odd);
        convert_to_deg(meas_x);

        BSpline fit_y(gt_y_even, measurement_y_even, i);
        meas_y = fit_y.calc(measurement_y_odd);
        convert_to_deg(meas_y);

        // error
        double err_x = rmse(gt_x_odd, meas_x);
        double err_y = rmse(gt_y_odd, meas_y);

        index.push_back(i);
        error_x.push_back(err_x);
        error_y.push_back(err_y);
        error_45_x.push_back(rmse_45(meas_x, gt_x_odd));
        error_45_y.push_back(rmse_45(meas_y, gt_y_odd));

        if (err_x < err_min_x) {
            err_min_x = err_x;
            err_min_index_x = i;
        }

        if (err_y < err_min_y) {
            err_min_y = err_y;
            err_min_index_y = i;
        }

        // size 
        size.push_back((sizeof(std::vector<double>) + sizeof(double) * fit_x.y.size()) * 2);

        // time
        auto start = timer::high_resolution_clock::now();
        for (int j = 0; j < 1000; ++j) {
            //fit_x.calc(measurement_y);
        }

        auto end = timer::high_resolution_clock::now();
        time.push_back(timer::duration_cast<timer::microseconds>(end - start).count() / 1000.0);

        printf("Degree %d, error x %f, error y %f, err45 x %f, err45 y %f, size %f, time %f\n", i, err_x, err_y, error_45_x.back(), error_45_y.back(), size.back(), time.back());
    plt::title("B-SPLINE");
    plt::grid(true);
    plt::named_plot("B-Spline X-Axis", gt_x_odd, meas_x);
    plt::named_plot("B-Spline Y-Axis", gt_y_odd, meas_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    plt::save("bsp_err_75.eps");
    plt::show();
    }
    printf("X: minimum error at degree %d with error %f\n", err_min_index_x, err_min_x);
    printf("Y: minimum error at degree %d with error %f\n", err_min_index_y, err_min_y);

    //error plot
    //plt::figure_size(800, 600);
    plt::title("B-SPLINE ERROR BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE", index, error_x);
    plt::named_plot("Y-RMSE", index, error_y);
    plt::ylim(0, 10);
    std::string filename = "../plots/BSpline_Count.eps";
    plt::ylim(0., 5.);
    plt::legend();
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("B-SPLINE ERROR BASED ON DEGREE (INNER 45 DEG)");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("RMSE (DEG)");
    plt::grid(true);
    plt::named_plot("X-RMSE 45", index, error_45_x);
    plt::named_plot("Y-RMSE 45", index, error_45_y);
    plt::legend();
    filename = "../plots/BSpline_45.eps";
    plt::save(filename);

    plt::cla();
    plt::clf();
    plt::title("B-SPLINE TIME BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("TIME [µs]");
    plt::grid(true);
    plt::plot(index, time);
    filename = "../plots/BSpline_time.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();
    plt::title("B-SPLINE SIZE BASED ON DEGREE");
    plt::xlabel("NUMBER OF KNOTS");
    plt::ylabel("SIZE IN BYTES");
    plt::grid(true);
    plt::plot(index, size);
    filename = "../plots/BSpline_size.eps";
    plt::save(filename);


    plt::cla();
    plt::clf();

    BSpline fit_x(gt_x_even, measurement_x_even, 50);
    meas_x = fit_x.calc(measurement_x_odd);
    convert_to_deg(meas_x);

    BSpline fit_y(gt_y_even, measurement_y_even, 50);
    meas_y = fit_y.calc(measurement_y_odd);
    convert_to_deg(meas_y);

    for (uint i = 0; i < meas_x.size(); ++i) {
        meas_err_x.push_back(gt_x_odd.at(i) - meas_x.at(i));
    }

    for (uint i = 0; i < meas_y.size(); ++i) {
        meas_err_y.push_back(gt_y_odd.at(i) - meas_y.at(i));
    }

    plt::title("B-SPLINE");
    plt::grid(true);
    plt::named_plot("B-Spline X-Axis", gt_x_odd, meas_x);
    plt::named_plot("error X-Axis", gt_x_odd, meas_err_x);
    plt::named_plot("B-Spline Y-Axis", gt_y_odd, meas_y);
    plt::named_plot("error Y-Axis", gt_y_odd, meas_err_y);
    plt::legend();
    plt::xlabel("INCIDENT LIGHT ANGLE [deg]");
    plt::ylabel("ANGULAR RESPONSE [deg]");
    filename = "../plots/BSpline.eps";
    plt::save(filename);
    plt::show();


    plt::cla();
    plt::clf();
}
