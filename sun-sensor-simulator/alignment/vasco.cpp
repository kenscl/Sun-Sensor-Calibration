#include "vasco.h"
#include <cstdio>
#include <vector>

Vasco::Vasco() {
  x = Vector<double>(48);
  for (uint i = 0; i < 48; ++i) {
    x.data.at(i) = .0;
  }

  P = Matrix<double>().Identity(48);
  Q = Matrix<double>().Identity(48);
  R = Matrix<double>().Identity(21);

  set_attitude_error_variance(1e1); // P
  set_gyro_error_variance(1e0);
  set_mag_error_variance(1e0);
  set_sun_error_variance(1e0);

  set_mag_measurement_variance(1e-5); // R
  set_sun_measurement_variance(1e-3);

  set_attitude_process_variance(1e-9); // Q
  set_gyro_process_variance(1e-16);
  set_mag_process_variance(1e-8);
  set_sun_process_variance(1e-8);


  sun_sensor_misalignment = std::vector<Matrix<double>>(6);
  sun_sensor_bias = std::vector<Vector<double>>(6);

  for (uint i = 0; i < sun_sensor_misalignment.size(); ++i) {
    sun_sensor_misalignment.at(i) = Matrix<double>().Identity(3);
    sun_sensor_bias.at(i) = Vector<double>(3);
  }

  mag_misalignment = Matrix<double>().Identity(3);
  mag_misalignment = Matrix<double>().Identity(3);
  mag_bias = Vector<double>(3);
  bias = Vector<double>(3);

  for (int i = 0; i< 6; ++i) {
    sun_sensor_misalignment.at(i) = Matrix<double>().Identity(3);
  }

}

Quaternion<double> Vasco::triad(Vector<double> sun_refrence, Vector<double> sun_measured, Vector<double> mag_refrence, Vector<double> mag_measured) {
  // triad algorithm
  Vector<double> S_hat = sun_refrence.normalize();
  Vector<double> s_hat = sun_measured.normalize();

  Vector<double> M_hat = (sun_refrence.dot(mag_refrence)).normalize();
  Vector<double> m_hat = (sun_measured.dot(mag_measured)).normalize();

  Matrix<double> refrence(3,3);
  refrence.set_col(0, S_hat);
  refrence.set_col(1, M_hat);
  refrence.set_col(2, S_hat.dot(M_hat));

  Matrix<double> measured(3,3);
  measured.set_col(0, s_hat);
  measured.set_col(1, m_hat);
  measured.set_col(2, s_hat.dot(m_hat));

  Matrix<double> A_hat = refrence * measured.transpose();
  return Quaternion<double>(A_hat);
}

void Vasco::propagate(Vector<double> w, double dt) {
  Matrix<double> omega_v(4,4);
  omega_v.data[0][1] = -w.data[0];
  omega_v.data[0][2] = -w.data[1];
  omega_v.data[0][3] = -w.data[2];

  omega_v.data[1][0] = w.data[0];
  omega_v.data[1][2] = w.data[2];
  omega_v.data[1][3] = -w.data[1];

  omega_v.data[2][0] = w.data[1];
  omega_v.data[2][1] = -w.data[2];
  omega_v.data[2][3] = w.data[0];

  omega_v.data[3][0] = w.data[2];
  omega_v.data[3][1] = w.data[1];
  omega_v.data[3][2] = -w.data[0];

  Matrix<double> identity = Matrix<double>().Identity(4);
  q_i_c = q_i_c + q_i_c * (omega_v * (0.5 * dt));
}


void Vasco::predict(double dt, Vector<double> w, Vector<double> mag_model, Vector<double> sun_model, int sun_sensor_id) {
  // determine F (A)
  /*
    * A: 
    * -[w]^x I3  | 0_6x(9*N)
    *  0_3 λ*I_3 |
    *  __________|_________
    *  0_(9*N)x6 | I_9*N
    *
    *  6 ss + 1 mag = 7;
    *  -> 69x69
    */

  // transform the sun and mag model to the local frame
  mag_model = mag_model.normalize();
  sun_model = sun_model.normalize();

  Quaternion<double> mag_q(mag_model);
  Quaternion<double> sun_q(sun_model);
  mag_q = q_i_c * mag_q * q_i_c.conjugate();
  sun_q = q_i_c * sun_q * q_i_c.conjugate();

  mag_model.data.at(0) = mag_q.i;
  mag_model.data.at(1) = mag_q.j;
  mag_model.data.at(2) = mag_q.k;

  sun_model.data.at(0) = sun_q.i;
  sun_model.data.at(1) = sun_q.j;
  sun_model.data.at(2) = sun_q.k;

  A = Matrix<double>(48, 48);
  A.data.at(0).at(1) = w.data.at(2);
  A.data.at(0).at(2) = -w.data.at(1);

  A.data.at(1).at(0) = -w.data.at(2);
  A.data.at(1).at(2) = w.data.at(0);

  A.data.at(2).at(0) = w.data.at(1);
  A.data.at(2).at(1) = -w.data.at(0);

  A.data.at(0).at(3) = 1.;
  A.data.at(1).at(4) = 1.;
  A.data.at(2).at(5) = 1.;

  double lamda = -.1; // gyro dependent, set to 1 here

  A.data.at(3).at(3) = lamda;
  A.data.at(4).at(4) = lamda;
  A.data.at(5).at(5) = lamda;

  //w.print();
  //for (int i = 0; i < 6; ++i) {
  //    for (int j = 0; j < 6; ++j) {
  //        printf("%f ", A.data.at(i).at(j));
  //    }
  //    printf("\n");
  //}


  // F is computed by matrix exponatniation, here it is truncated to 3st order
  A = A * dt;
  F = Matrix<double>().Identity(48) + A + A * A * 0.5 + A * A * A * (1 / 6.0);

  // determine H 
  /*
    * H:
    * -[1]^x 0_3 [1]^* I_3 
    * -[2]^x 0_3 ... 
    * -[3]^x 0_3 ...
    * -[4]^x 0_3 ...
    * -[5]^x 0_3 ...
    * -[6]^x 0_3 ...
    * -[7]^x 0_3 ... [7]^* I_3
    *  -> 7*3 x 3 + 3 + (6 + 3) * 7 = 21 x 69
    */

  // first populate the skews 

  H = Matrix<double>(21, 48);

  // mag is sensor 1 

  H.data.at(0).at(1) = mag_model.data.at(2);
  H.data.at(0).at(2) = - mag_model.data.at(1);

  H.data.at(1).at(0) = - mag_model.data.at(2);
  H.data.at(1).at(2) = mag_model.data.at(0);

  H.data.at(2).at(0) = mag_model.data.at(1);
  H.data.at(2).at(1) = - mag_model.data.at(0);

  // sun is sensor 1 + id 
  int ss_col = 3 + 3 * (sun_sensor_id - 1);
  int ss_row = 6 + 6 + (sun_sensor_id - 1) * 6;
  if (sun_sensor_id < 7) {
    H.data.at(ss_col).at(1) = sun_model.data.at(2);
    H.data.at(ss_col).at(2) = - sun_model.data.at(1);

    H.data.at(ss_col + 1).at(0) = - sun_model.data.at(2);
    H.data.at(ss_col + 1).at(2) = sun_model.data.at(0);

    H.data.at(ss_col + 2).at(0) = sun_model.data.at(1);
    H.data.at(ss_col + 2).at(1) = - sun_model.data.at(0);
  }

  // now get stars
  H.data.at(0).at(6) = mag_model.data.at(1);
  H.data.at(0).at(7) = mag_model.data.at(2);

  H.data.at(1).at(8) = mag_model.data.at(0);
  H.data.at(1).at(9) = mag_model.data.at(2);

  H.data.at(2).at(10) = mag_model.data.at(0);
  H.data.at(2).at(11) = mag_model.data.at(1);


  if (sun_sensor_id < 7) {
    H.data.at(ss_col).at(ss_row) = sun_model.data.at(1);
    H.data.at(ss_col).at(ss_row + 1) = sun_model.data.at(2);

    H.data.at(ss_col + 1).at(ss_row + 2) = sun_model.data.at(0);
    H.data.at(ss_col + 1).at(ss_row + 3) = sun_model.data.at(2);

    H.data.at(ss_col + 2).at(ss_row + 4) = sun_model.data.at(0);
    H.data.at(ss_col + 2).at(ss_row + 5) = sun_model.data.at(1);
  }


  //for (int i = 0; i < 3; ++i) {
  //    for (int j = 0; j < 15; ++j) {
  //        printf("%f ", H.data.at(i).at(j));
  //    }
  //    printf("\n");
  //}
  //    printf("\n");
  //if (sun_sensor_id > 5 && sun_sensor_id < 7) {
  //sun_model.print();
  //for (int i = 15; i < 21; ++i) {
  //    for (int j = 0; j < 48; ++j) {
  //        if (j > 5 && j < 36) {
  //        } else {
  //          printf("%f ", H.data.at(i).at(j));
  //        }
  //    }
  //    printf("\n");
  //}
  //}
  x_hat = F * x; 
  P_bar = F * P * F.transpose() + Q ;
}

void Vasco::update(Vector<double> mag_measured, Vector<double> sun_measured, int sun_sensor_id, Vector<double> mag_refrence, Vector<double> sun_refrence) {
  Vector<double> y = Vector<double>(21);
  mag_measured = mag_measured.normalize();
  mag_refrence = mag_refrence.normalize();
  mag_refrence = q_i_c.to_rotation_matrix().inverse() * mag_refrence;
  mag_measured = mag_measured.normalize();
  mag_measured = mag_measured - mag_refrence;

  sun_measured = sun_measured.normalize();
  sun_refrence = sun_refrence.normalize();
  sun_refrence = q_i_c.to_rotation_matrix().inverse() * sun_refrence; 
  sun_measured = sun_measured.normalize();
  sun_measured = sun_measured - sun_refrence;

  y = y.set_sub_vector(mag_measured, 0, 3);
  if (sun_sensor_id < 7) {
    y = y.set_sub_vector(sun_measured, 3 + (sun_sensor_id - 1) * 3, 6 + (sun_sensor_id - 1) * 3);
  }

  K = P_bar * H.transpose() * (H * P_bar * H.transpose() + R).inverse();
  residuals.push_back(y - H * x_hat);
  x = x_hat + K * (y - H * x_hat);

  P = (Matrix<double>().Identity(48) - K * H) * P_bar;
}

void Vasco::set_attitude_error_variance(double variance) {
  Vector<double> var(3);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  P.diag(P.diag().set_sub_vector(var, 0, 3));
}

void Vasco::set_gyro_error_variance(double variance) {
  Vector<double> var(3);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  P.diag(P.diag().set_sub_vector(var, 3, 6));
}

void Vasco::set_mag_error_variance(double variance) {
  Vector<double> var(6);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  P.diag(P.diag().set_sub_vector(var, 6, 12));
}

void Vasco::set_sun_error_variance(double variance) {
  Vector<double> var(36);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  P.diag(P.diag().set_sub_vector(var, 12, 48));
}

void Vasco::set_mag_measurement_variance(double variance) {
  Vector<double> var(3);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  R.diag(R.diag().set_sub_vector(var, 0, 3));
}

void Vasco::set_sun_measurement_variance(double variance) {
  Vector<double> var(18);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  R.diag(R.diag().set_sub_vector(var, 3, 21));
}

void Vasco::set_attitude_process_variance(double variance) {
  Vector<double> var(3);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  Q.diag(Q.diag().set_sub_vector(var, 0, 3));
}
void Vasco::set_gyro_process_variance(double variance) {
  Vector<double> var(3);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  Q.diag(Q.diag().set_sub_vector(var, 3, 6));
}
void Vasco::set_mag_process_variance(double variance) {
  Vector<double> var(6);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  Q.diag(Q.diag().set_sub_vector(var, 6, 12));
}

void Vasco::set_sun_process_variance(double variance) {
  Vector<double> var(36);
  for (uint i = 0; i < var.size; ++i) {
    var.data.at(i) = variance;
  }
  Q.diag(Q.diag().set_sub_vector(var, 12, 48));
}

bool Vasco::has_converged() {
  Vector<double> diagonal = P.diag();
  double sum = 0;
  for (uint i = 0; i < diagonal.size; ++i) {
    sum += diagonal.data.at(i);
    if (diagonal.data.at(i) > convergence_threshold) {
      return false;
    }
  }
  //printf("Sum %f\n", sum);
  return true;
}

void Vasco::feed_forward() {
  double angle = x.sub_vector(0, 3).norm();
  Vector<double> axis = x.sub_vector(0, 3).normalize();
  Quaternion<double> q_theta = Quaternion<double>(angle, axis);
  q_bk = q_i_c * q_theta.conjugate();
}

void Vasco::full_feedback(double elapsed_time) {
  if (!has_converged()) {
    return;
  }


  feed_forward();

  // gyro
  bias = x.sub_vector(3, 6);

  // mag
  mag_misalignment.data.at(0).at(1) += x.data.at(6);
  mag_misalignment.data.at(0).at(2) += x.data.at(7);
  mag_misalignment.data.at(1).at(0) += x.data.at(8);
  mag_misalignment.data.at(1).at(2) += x.data.at(9);
  mag_misalignment.data.at(2).at(0) += x.data.at(10);
  mag_misalignment.data.at(2).at(1) += x.data.at(11);

  // sun
  for (uint i = 0; i < 6; ++i) {
    sun_sensor_misalignment.at(i).data.at(0).at(1) += x.data.at(12 + i * 6);
    sun_sensor_misalignment.at(i).data.at(0).at(2) += x.data.at(13 + i * 6);
    sun_sensor_misalignment.at(i).data.at(1).at(0) += x.data.at(14 + i * 6);
    sun_sensor_misalignment.at(i).data.at(1).at(2) += x.data.at(15 + i * 6);
    sun_sensor_misalignment.at(i).data.at(2).at(0) += x.data.at(16 + i * 6);
    sun_sensor_misalignment.at(i).data.at(2).at(1) += x.data.at(17 + i * 6);
  }
  x = Vector<double>(48);

  q_i_c = q_bk.normalize();
}

void Vasco::partial_feedback() {
  if (!has_converged()) {
    return;
  }

  feed_forward();

  // restet the error state
  x.data.at(0) = 0;
  x.data.at(1) = 0;
  x.data.at(2) = 0;

  q_i_c = q_bk;
}
  

Vector<double> Vasco::correct_gyro(Vector<double> w) { return w - bias; }
Vector<double> Vasco::correct_mag(Vector<double> mag) {
  return mag_misalignment.inverse() * mag;
}
Vector<double> Vasco::correct_sun(Vector<double> sun, int id) {
  if (id > 6) {
    return sun;
  }
  return sun_sensor_misalignment.at(id - 1).inverse() * sun;
}

void Vasco::correct_matricies() {
  mag_misalignment.data[1][0] -= mag_misalignment.data[0][1];
  mag_misalignment.data[1][0] /= 2;
  mag_misalignment.data[0][1] = - mag_misalignment.data[1][0];

  mag_misalignment.data[2][0] -= mag_misalignment.data[0][2];
  mag_misalignment.data[2][0] /= 2;
  mag_misalignment.data[0][2] = - mag_misalignment.data[2][0];

  mag_misalignment.data[2][1] -= mag_misalignment.data[1][2];
  mag_misalignment.data[2][1] /= 2;
  mag_misalignment.data[1][2] = - mag_misalignment.data[2][1];

  for (int i = 0; i < 6; ++i) {
    sun_sensor_misalignment.at(i).data[1][0] -= sun_sensor_misalignment.at(i).data[0][1];
    sun_sensor_misalignment.at(i).data[1][0] /= 2;
    sun_sensor_misalignment.at(i).data[0][1] = - sun_sensor_misalignment.at(i).data[1][0];

    sun_sensor_misalignment.at(i).data[2][0] -= sun_sensor_misalignment.at(i).data[0][2];
    sun_sensor_misalignment.at(i).data[2][0] /= 2;
    sun_sensor_misalignment.at(i).data[0][2] = - sun_sensor_misalignment.at(i).data[2][0];

    sun_sensor_misalignment.at(i).data[2][1] -= sun_sensor_misalignment.at(i).data[1][2];
    sun_sensor_misalignment.at(i).data[2][1] /= 2;
    sun_sensor_misalignment.at(i).data[1][2] = - sun_sensor_misalignment.at(i).data[2][1];
  }
}

void Vasco::run(std::vector<double> dt, std::vector<Vector<double>> w, std::vector<Vector<double>> mag_model, std::vector<Vector<double>> mag_measured, std::vector<Vector<double>> sun_model, std::vector<Vector<double>> sun_measured, std::vector<int> sun_sensor_id, double threashold, std::vector<Quaternion<double>> attitude) {
  q_i_c = triad(sun_model[0], sun_measured[0], mag_model[0], mag_measured[0]);
  double elapsed_time = 0;
  double partial_time = 0;
  double full_time = 0;
  convergence_threshold = threashold;
  for (uint i = 0; i < w.size(); ++i) {
    propagate(correct_gyro(w.at(i)), dt.at(i));
    predict(dt.at(i), correct_gyro(w.at(i)), mag_model.at(i), sun_model.at(i), sun_sensor_id.at(i));
    update(correct_mag(mag_measured.at(i)), correct_sun(sun_measured.at(i), sun_sensor_id.at(i)), sun_sensor_id.at(i) , mag_model.at(i), sun_model.at(i));


    if (elapsed_time == 1500) {
      partial_feedback();
      partial_time = 0;
    }

    if (elapsed_time > 1500) {
      if (full_time > 500) {
        full_feedback(elapsed_time);
        full_time = 0;
        partial_time = 0;
      }

      if (partial_time > 50) {
        partial_feedback();
        partial_time = 0;
      }

    }
    if (elapsed_time > 21000) {
      break;
    }

    store_data(elapsed_time);
    print(elapsed_time);

    printf("Time %f [s]\n", elapsed_time);
    elapsed_time += dt.at(i);
    partial_time += dt.at(i);
    full_time += dt.at(i);
  } 
  full_feedback(elapsed_time);
}

void Vasco::print(double time) {
  printf("Time %f [s]\n", time);
  x.print();
  printf("Attitude \n");
  q_i_c.print();
  printf("Bias\n");
  bias.print();
  printf("Mag Misalignment\n");
  mag_misalignment.print();
  printf("Sun Misalignment\n");
  for (uint i = 0; i < 6; ++i) {
    sun_sensor_misalignment.at(i).print();
  }
}


Matrix<double> Vasco::b_to_c() {
  Matrix<double> b_to_c = Matrix<double>().Identity(3);
  b_to_c.data.at(0).at(1) = x.data.at(2);  
  b_to_c.data.at(0).at(2) = -x.data.at(1);
  b_to_c.data.at(1).at(0) = -x.data.at(2);
  b_to_c.data.at(1).at(2) = x.data.at(0);
  b_to_c.data.at(2).at(0) = x.data.at(1);
  b_to_c.data.at(2).at(1) = -x.data.at(0);
  return b_to_c;
}


void Vasco::store_data(double time) {
  elapsed_time.push_back(time);
  mag_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + mag_misalignment));
  sun_sensor_xp_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(0)));
  sun_sensor_xm_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(1)));
  sun_sensor_yp_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(2)));
  sun_sensor_ym_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(3)));
  sun_sensor_zp_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(4)));
  sun_sensor_zm_misalignment_estimates.push_back(Quaternion<double>(Matrix<double>().Identity(3) + sun_sensor_misalignment.at(5)));
  cov_matrix_diag.push_back(P.diag());
  attitude_estimates.push_back(q_i_c);
}
