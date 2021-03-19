#ifndef MODEL_H
#define MODEL_H

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

class Model {

    public:
        
        // Attributes
        int size; // The number of models.
        int counter; // The number of vertices.
        std::vector< std::array<double,4> > axisAngleRotation;
        std::vector< double > scalingFactor;
        std::vector< std::array<double,3> > translation;
        Eigen::Matrix4d Matrix;
        std::vector<bool> smooth;
        std::vector< std::string > modelName;
        std::ifstream model;
        std::vector< std::string > materialFileName;
        std::ifstream material;
        std::vector< std::string > materialssss;
        std::vector< double > Ns;
        std::vector< std::array<double,3> > Ka;
        std::vector< std::array<double,3> > Kd;
        std::vector< std::array<double,3> > Ks;
        std::vector< std::array<double,3> > vertices;
        int materialId;
        std::vector< std::array<int,4> > triangles;
        std::vector< Eigen::Vector3d > sn;
        std::vector<int> record; // Helper for converting the indices for vertices.
        
        double best = std::numeric_limits<double>::max(); // The closest distance between the starting location of a ray and the closest triangle. 
        bool flag; // If a ray intersects all the triangles.
        Eigen::Vector3d bestSurfaceNormal;// The surface normal of the closest triangle.
        int bestModel = 0; // The index of the model which has the closest triangle. 
        int bestTriangle = 0;
        
        // Methods
        int readDriver(std::ifstream& driver);
        int readModel(std::ifstream& model);
        void readMaterial(std::ifstream& material);
        void transformationMatrix(std::array<double,4> axisAngleRotation, double scalingFactor, std::array<double,3> translation); 
        void transform(std::array<double,3> &arr); 
        void rayTriangleIntersection(Eigen::Vector3d pixelPt, Eigen::Vector3d shoot);
};

#endif
