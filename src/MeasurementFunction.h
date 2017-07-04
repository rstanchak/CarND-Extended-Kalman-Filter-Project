#ifndef MEASUREMENT_FUNCTION_H_
#define MEASUREMENT_FUNCTION_H_

#include <iostream>
#include "Eigen/Dense"

class MeasurementFunction {
 public:
    virtual Eigen::VectorXd residual(const Eigen::VectorXd & z, const Eigen::VectorXd & x) {
		return z-evaluate(x);
	}
    virtual Eigen::VectorXd evaluate(const Eigen::VectorXd & x) = 0; // => z
    virtual const Eigen::MatrixXd& jacobian(const Eigen::VectorXd& x) = 0;
};

class RadarMeasurementFunction : public MeasurementFunction {
	// rho, phi, rho-dot
 private:
    Eigen::MatrixXd Hj_;

 public:
    RadarMeasurementFunction(void): Hj_(3, 4) {
    }

	double normalize_angle(double theta) {
		return atan2( sin(theta), cos(theta) );
	}

    Eigen::VectorXd residual(const Eigen::VectorXd & z, const Eigen::VectorXd & x) {
		Eigen::VectorXd y = z-evaluate(x);
		y[1] = normalize_angle(y[1]);
		return y;
	}
    Eigen::VectorXd evaluate(const Eigen::VectorXd & x) {
        double px, py, vx, vy;
        double rho;
        px = x[0];
        py = x[1];
        vx = x[2];
        vy = x[3];
        rho = sqrt(px*px + py*py);
        assert(px != 0 && py != 0);
        return Eigen::Vector3d(
            rho,
            atan2(py, px), /* phi */
            (px*vx + py*vy)/rho /* rho-dot */);
    }
    const Eigen::MatrixXd& jacobian(const Eigen::VectorXd& x) {
        // recover state parameters
        double px = x(0);
        double py = x(1);
        double vx = x(2);
        double vy = x(3);


        // check division by zero
        if (px == 0 && py == 0) {
            assert(px != 0 && py != 0);
            std::cerr << "ERROR: divide by zero" << std::endl;
        }

        // compute the Jacobian matrix
        double rho2 = px*px+py*py;
        double rho = sqrt(rho2);
		double rho3_2 = rho2*rho;

        Hj_ << px/rho, py/rho, 0, 0,
            -py/rho2, px/rho2, 0, 0,
            py*(vx*py-vy*px)/rho3_2, px*(vy*px-vx*py)/rho3_2, px/rho, py/rho;
        return Hj_;
    }
};

#endif // MEASUREMENT_FUNCTION_H_
