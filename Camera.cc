#include "Camera.h"

using namespace std;
using namespace Eigen;
    
/*
 * This method is to read the camera specification from the driver file.
 * 
 * The general format for the camera specification is as follows: 
 * eye ex ey ez 
 * look x y z
 * up x y z
 * d distance
 * bounds left bottom right top 
 * res width height 
 * 
 * "eye" represents the location of the focal point (the Eye). 
 * "look" represents the look at point.
 * "up" represents the up vector.
 * "d" represents the focal length - the distance from the focal point to the image plane (near clipping plane).
 * "bounds" values indicate the minimum and maximum extend of the bounded image rectangle 
 *  on the infinite image plane in the camera horizontal and vertical directions respectively.
 * "res" values separately indicate the pixel sampling resolution across the horizontal 
 *  and vertical dimensions of the bounded rectangle.
 */ 
void Camera::readDriver(ifstream& driver){
    string line = "";
    while(getline(driver, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id[0] == '#') continue;
            else if(id == "eye"){
                
                /*
                 * Read the focal point.
                 */
                for(int i = 0; i < 3; i++){iss >> eye[i];}
            }
            else if(id == "look"){
                
                /*
                 * Read the look at point.
                 */
                for(int i = 0; i < 3; i++){iss >> look[i];}
            }
            else if(id == "up"){
                
                /*
                 * Read the up vector.
                 */
                for(int i = 0; i < 3; i++){iss >> up[i];}
            }
            else if(id == "d"){
                
                /*
                 * Read the focal length.
                 */
                iss >> d;
            }
            else if(id == "bounds"){
                
                /*
                 * Read the bounds of the image plane.
                 */
                for(int i = 0; i < 4; i++){iss >> bounds[i];}
            }
            else if(id == "res"){
                
                /*
                 * Read the resolution.
                 */
                for(int i = 0; i < 2; i++){iss >> res[i];}
            } else {continue;}
        } else {continue;}
    }
    driver.close(); 
}
    
    
/*
 * This method place the camera.
 */    
void Camera::placement(){

    E << eye[0], eye[1], eye[2];
    Vector3d L;
    L << look[0], look[1], look[2];
    Vector3d UP;
    UP << up[0], up[1], up[2];
    W = E - L; 
    W.normalize();
    U = UP.cross(W); 
    U.normalize();
    V = W.cross(U);
    
}    
