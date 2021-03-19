#ifndef LIGHT_H
#define LIGHT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <math.h>
#include <limits>
#include <Eigen/Dense>
#include "Model.h"
#include "Sphere.h"

class Light {

    public:
        
        // Attributes
        double ambient[3];
        std::vector< std::array<double,7> > diffuse;
        Eigen::Vector3d color;
        int depth;
        double eta_outside = 1.0;
        
        bool flag;
        bool isSphere;
        double t = std::numeric_limits<double>::max();
        

        // Methods
        void readDriver(std::ifstream& driver);
       
        void illumination(Eigen::Vector3d &pixelPt, Eigen::Vector3d &shoot, Eigen::Vector3d &color, Model& model, Sphere& sphere, double (&reflectionAttenuation)[3], int level);
        void rayFind(Eigen::Vector3d &pixelPt, Eigen::Vector3d &shoot, Model& model, Sphere& sphere);
};

#endif
