#include "Sphere.h"

using namespace std;
using namespace Eigen;

/*
 * This method is to read the material information of each sphere from the driver file.
 * 
 * The general format for each sphere is as follows: 
 * sphere x y z radius; Ka_r Ka_g Ka_b; Kd_r Kd_g Kd_b; Ks_r Ks_g Ks_b; Kr_r Kr_g Kr_b
 * 
 * The first three values are the x, y, and z coordinates of the center of the sphere in world coordinates. 
 * The fourth value is the radius of the sphere. 
 * First triplet after the radius represents (Ka_red, Ka_green, Ka_blue) ambient coefficients.
 * Next triplet represents diffuse (Kd_red, Kd_green, Kd_blue) coefficients.
 * Next we have (Ks_red, Ks_green, Ks_blue) specular coefficients. 
 * Finally, we have (Kr_red, Kr_green, Kr_blue) attenuation coefficients.
 * 
 * For spheres use phong constant = 16.
 * 
 * There could be zero or more spheres. 
 */ 
void Sphere::readDriver(ifstream& driver){
    size = 0;
    
    string line = "";
    while(getline(driver, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id == "sphere"){
                
                /*
                 * Sphere center and radius: x y z radius. 
                 */
                array<double,3> cen;
                for(int i = 0; i < 3; i++){iss >> cen[i];}

                Vector3d Cen;
                Cen << cen[0], cen[1], cen[2];
                Center.push_back(Cen);
                
                double rad;
                iss >> rad;
                radius.push_back(rad);
                
                /*
                 * Read the ambient coefficients.
                 */
                array<double,3> amb;
                for(int i = 0; i < 3; i++){iss >> amb[i];}
                sphereAmbient.push_back(amb);
                
                /*
                 * Read the diffuse coefficients.
                 */
                array<double,3> dif;
                for(int i = 0; i < 3; i++){iss >> dif[i];}
                sphereDiffuse.push_back(dif);
                
                /*
                 * Read the specular coefficients.
                 */
                array<double,3> spe;
                for(int i = 0; i < 3; i++){iss >> spe[i];}
                sphereSpecular.push_back(spe);
                
                /*
                 * Read the attenuation coefficients.
                 */
                array<double,3> att;
                for(int i = 0; i < 3; i++){iss >> att[i];}
                sphereAttenuation.push_back(att);
                
                /*
                 * Read the opacity coefficients.
                 */
                array<double,3> opa;
                for(int i = 0; i < 3; i++){iss >> opa[i];}
                sphereOpacity.push_back(opa);
                
                /*
                 * Read the exponent used to control the apparent size of specular highlights. 
                 */
                double spow;
                iss >> spow;
                spherePhong.push_back(spow);
                
                /*
                 * Read the index of refraction for the material.
                 * 1.0 for air and typically 1.5 for glass. 
                 */
                double refraction;
                iss >> refraction;
                eta.push_back(refraction);
                
                /*
                 * The number of spheres.
                 */
                size++;
                
            } else {continue;}
        } else {continue;}
    }
    
    driver.close(); 
}
    
    
/*
 * This method does Ray-Sphere-Intersection and get the minimum distance from
 * the starting location of the ray to all the spheres.
 * By the way, find the closest sphere.
 * 
 * The two parameters are the starting location and the direction of the ray. 
 */     
void Sphere::raySphereIntersection(Vector3d pixelPt, Vector3d shoot){
            flag = false;
            best = numeric_limits<double>::max();
            
            for(int k = 0; k < size; k++){
                // For each shpere
                //test code/////////////////////////////////////////////////////////////////////
                //cout << "Run the shpere test one time!\n";
                //cout << "so what is Ray L: " << pixelPt[0] << ", "
                //        << pixelPt[1] << ", " << pixelPt[2] << "\n";
                //cout << "so what is Ray D: " << shoot[0] << ", "
                //    << shoot[1] << ", " << shoot[2] << "\n";
                
                Vector3d LC;
                LC = Center[k] - pixelPt;
                double v = LC.dot(shoot);
                double cSquare = LC.dot(LC);
                double dSquare = radius[k] * radius[k] - (cSquare -v*v);
                if (dSquare > 0){
                    double tSphere = v - sqrt(dSquare);
                    if((tSphere > 0.00001) && (tSphere < best)){
                        flag = true;
                        best = tSphere;
                        bestHit = pixelPt + tSphere * shoot;
                        bestSphere = k;
                        //test code/////////////////////////////////////////////////////////////////////
                        //if(flag == true)
                        //cout << "HIT ONE TIME!\n";
                        //cout << "hit distance: " << best << "\n";
                    }
                }
            }
}    

/*
 * This method is to get the refraction ray at the first time. (from outside the shpere to inside the sphere)
 */    
Vector3d Sphere::refractionRay(Vector3d reverseIncidentRay,Vector3d hit,Vector3d N,double eta_out,double eta_in){
    //test code/////////////////////////////////////////////////////////////////////
    //cout << "Run first refraction one time!\n";
    //cout << "reverseIncidentRay: " << reverseIncidentRay[0] << ", "<< reverseIncidentRay[1] << ", " << reverseIncidentRay[2] << "\n";
            
    Vector3d T;
    T << 0.0,0.0,0.0;         
    double   etaRatio = eta_out / eta_in;
    double   alpha = -etaRatio;
    double   wDotN = reverseIncidentRay.dot(N);
    //test code/////////////////////////////////////////////////////////////////////
    //cout << "W: " << reverseIncidentRay[0] << ", "<< reverseIncidentRay[1] << ", " << reverseIncidentRay[2] << "\n";
    //cout << "N: " << N[0] << ", "<< N[1] << ", " << N[2] << "\n";
    //cout << "w dot n" << wDotN << "\n";
    
    double   square = pow(etaRatio, 2) * (pow(wDotN, 2) - 1) + 1;
    //cout << "square" << square << "\n";
    
    if(square < 0.0){
        return T;
    } else {
        double beta = (etaRatio * wDotN) - sqrt(square);
        T =  alpha * reverseIncidentRay + beta * N;
        //test code/////////////////////////////////////////////////////////////////////
        //cout << "alpha: " << alpha << "\n";
        //cout << "beta: " << beta << "\n";
        //cout << "T: " << T[0] << ", "<< T[1] << ", " << T[2] << "\n";
        return T;
    }
}

/*
 *  This method is to get the refraction ray at the second time.(from inside the sphere to outs the sphere)
 */    
Vector3d Sphere::refractionExit(int tempBestSphere, Vector3d reverseIncidentRay, Vector3d Hit, double eta_outside){
    //test code/////////////////////////////////////////////////////////////////////
    //cout << "Run second refraction one time!\n";
            
    exitRay << 0.0,0.0,0.0;  
    Vector3d N = Hit - Center[tempBestSphere]; 
    N.normalize();
    Vector3d refractionRayT = refractionRay(reverseIncidentRay, Hit, N ,eta_outside, eta[tempBestSphere]);
    bool cannot = ((refractionRayT[0] + refractionRayT[1] + refractionRayT[2]) == 0.0);
    if(!cannot){
        exitPoint = Hit + 2 * (Center[tempBestSphere] - Hit).dot(refractionRayT) * refractionRayT;
        Vector3d InsideSurfaceNormal = Center[tempBestSphere] - exitPoint;
        InsideSurfaceNormal.normalize();
        exitDirection = refractionRay(-1 * refractionRayT, exitPoint, InsideSurfaceNormal, eta[tempBestSphere], eta_outside);
        exitDirection.normalize();
        exitRay = exitPoint + exitDirection;
        //test code/////////////////////////////////////////////////////////////////////
        //cout << "exitPoint: " << exitPoint[0] << ", "<< exitPoint[1] << ", " << exitPoint[2] << "\n";
        //cout << "exitDirection: " << exitDirection[0] << ", "<< exitDirection[1] << ", " << exitDirection[2] << "\n";
    }
    return exitRay;
}
