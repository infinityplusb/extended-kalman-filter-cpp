#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    if(estimations.size() != ground_truth.size() || estimations.size() == 0)
    {
      return rmse;
    }

    for(unsigned int i = 0; i < estimations.size(); i++ )
    {
      // E
      VectorXd residuals = estimations[i] - ground_truth[i];
      // S
      residuals = residuals.array() * residuals.array() ;
      rmse += residuals;
    }

    // M
    rmse = rmse / estimations.size();

    // R
    rmse = rmse.array().sqrt();

    return rmse ;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
   // cout << "Calculating the Jacobian! " << endl;
   MatrixXd Hj(3,4) ;

   // extract the state being passed
   float p_x = x_state(0);
   float p_y = x_state(1);
   float v_x = x_state(2);
   float v_y = x_state(3);

   float a1 = p_x*p_x + p_y*p_y;
   float b1 = sqrt(a1);
   // cout << a1 << endl;
   // if we can't devide by zero
   if(fabs(a1) < 0.0001)
   {
     // cout << "Returning zero! " << endl;
     return MatrixXd::Zero(3,4);
   }

   Hj << (p_x/b1),      (p_y/b1),      0,      0,
        -(p_y/a1),      (p_x/a1),      0,      0,
        p_y*(v_x*p_y - v_y*p_x)/(a1*b1),
          p_x*(p_x*v_y - p_y*v_x)/(a1*b1),
          p_x/b1,
          p_y/b1;

  // cout << "Returning the Jacobian! " << endl;
   return Hj;
}
