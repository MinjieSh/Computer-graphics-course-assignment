#ifndef CAMERA_H
#define CAMERA_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <math.h>
#include <Eigen/Dense>

class Camera {

    public:
        
        // Attributes
        double eye[3];
        double look[3];
        double up[3];
        double d = 0;
        double bounds[4];
        int res[2];
        Eigen::Vector3d E;
        Eigen::Vector3d W;
        Eigen::Vector3d U;
        Eigen::Vector3d V;

        // Methods
        void readDriver(std::ifstream& driver);
        void placement();

};

#endif
