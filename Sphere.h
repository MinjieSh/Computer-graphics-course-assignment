#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <math.h>
#include <Eigen/Dense>
#include <limits>

class Sphere {

    public:
    
        // Attributes
        std::vector< Eigen::Vector3d > Center;
        std::vector< double > radius;
        std::vector< std::array<double,3> > sphereAmbient;
        std::vector< std::array<double,3> > sphereDiffuse;
        std::vector< std::array<double,3> > sphereSpecular;
        std::vector< std::array<double,3> > sphereAttenuation;
        std::vector< std::array<double,3> > sphereOpacity;
        std::vector<double> spherePhong;
        std::vector<double> eta;
        Eigen::Vector3d exitPoint;
        Eigen::Vector3d exitDirection;
        Eigen::Vector3d exitRay;
        
        int size; // The number of spheres.
        double best = std::numeric_limits<double>::max();
        bool flag;
        int bestSphere = 0;
        Eigen::Vector3d bestHit; // The hit position of the closest shpere.
        
        // Methods
        void readDriver(std::ifstream& driver);
        void raySphereIntersection(Eigen::Vector3d pixelPt, Eigen::Vector3d shoot);
        Eigen::Vector3d refractionRay(Eigen::Vector3d reverseIncidentRay, Eigen::Vector3d hit, Eigen::Vector3d N, double eta_out, double eta_in);
        Eigen::Vector3d refractionExit(int tempBestSphere, Eigen::Vector3d reverseIncidentRay,Eigen::Vector3d Hit,double eta_outside);
};

#endif
