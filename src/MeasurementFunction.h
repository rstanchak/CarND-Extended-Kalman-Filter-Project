#ifndef MEASUREMENT_FUNCTION_H
#define MEASUREMENT_FUNCTION_H

#include <iostream>
#include "Eigen/Dense"

class MeasurementFunction {
public:
    virtual Eigen::VectorXd evaluate( const Eigen::VectorXd & x ) = 0; // => z
    virtual const Eigen::MatrixXd& jacobian( const Eigen::VectorXd& x ) = 0;
};

class RadarMeasurementFunction : public MeasurementFunction {
private:
    Eigen::MatrixXd Hj_;

public:
    RadarMeasurementFunction(void): Hj_(3,4) {
    }
    Eigen::VectorXd evaluate( const Eigen::VectorXd & x ){
        double px, py, vx, vy;
        double rho;
        px = x[0];
        py = x[1];
        vx = x[2];
        vy = x[3];
        rho = sqrt( px*px + py*py );
        assert(px!=0 && py!=0);
        return Eigen::Vector3d( 
            rho, 
            atan2( py, px ), /* phi */
            (px*vx + py*vy)/rho /* rho-dot */);
    } 
    const Eigen::MatrixXd& jacobian( const Eigen::VectorXd& x) {
        //recover state parameters
        float px = x(0);
        float py = x(1);
        float vx = x(2);
        float vy = x(3);


        //check division by zero
        if(px==0 && py==0)
        {
            assert(px!=0 && py!=0);
            std::cerr<<"ERROR: divide by zero"<<std::endl;
        }

        //compute the Jacobian matrix
        float rho2 = px*px+py+py;
        float rho = sqrt(rho2);

        Hj_ << px/rho, py/rho, 0, 0,
            -py/rho2, px/rho2, 0, 0,
            py*(vx*py-vy*px)/pow(rho2,1.5), px*(vy*px-vx*py)/pow(rho2,1.5), px/rho, py/rho;
        return Hj_;
    }
};

#endif /* MEASUREMENT_FUNCTION_H */
