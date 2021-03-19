#include "Light.h"

using namespace std;
using namespace Eigen;


/*
 * This method is to read the information of each light from the driver file.
 * 
 * The general format for each model is as follows: 
 * "recursionLevel depth"
 * "ambient r g b"
 * "light x y z w r g b"
 * 
 * "recursionLevel" specifies the recursion depth.
 * "ambient" specifies the ambient illumination in the scene.
 * The first four values given are the x, y, z, and w coordinates of the light source in world coordinates.
 * "w" is generally one, but a zero indicates a light source at infinity in the direction specified by x, y, and z.
 * The last three values indicate the red, green, and blue levels of the light source on a zero to one scale.
 * 
 * There would be multiple light sources in a driver file.
 */ 
void Light::readDriver(ifstream& driver){
    string line = "";
    while(getline(driver, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id == "ambient"){
                
                /*
                 * Read the ambient illumination in the scene.
                 */
                for(int i = 0; i < 3; i++){iss >> ambient[i];}
                
            } else if(id == "recursionLevel"){
                
                /*
                 * Read the recursion depth.
                 */
                iss >> depth;
                
            } else if(id == "light"){
                
                /*
                 * Read the light sources one by one.
                 */
                array<double,7> light;
                for(int i = 0; i < 7; i++){iss >> light[i];}
                diffuse.push_back(light);
                
            } else if(id == "eta_outside"){
                
                /*
                 * Read the eta_outside.
                 */
                iss >> eta_outside;
                
            } else {continue;}
        } else {continue;}
    }
    
    driver.close(); 
}
    
/*
 * This method is a recursive ray tracer which illuminates each pixel.
 * It implements ambient illumination, diffuse illumination and specular illumination for each pixel.
 * It also implements reflection in terms of the recursion depth.
 */    
void Light::illumination(Vector3d &pixelPt, Vector3d &shoot, Vector3d &accumulation, Model& model, Sphere& sphere, double (&reflectionAttenuation)[3], int levelParameter){
            
            //test code/////////////////////////////////////////////////////////////////////
            //cout << "level" << levelParameter << "\n";
            
            rayFind(pixelPt, shoot, model, sphere);
            
            Vector3d colorTemp;
            if(flag == true){
                //test code/////////////////////////////////////////////////////////////////////
                //cout << "into illumination\n";
                int bestSphere = 0;
                //int bestModel = 0;
                double Ns = 0;
                double Ka[3];
                double Kd[3];
                double Ks[3];
                double Kr[3];
                double Ko[3];
                Eigen::Vector3d sn_min;
                
                if(isSphere){
                    sn_min = sphere.bestHit - sphere.Center[sphere.bestSphere]; 
                    sn_min.normalize();
                    Ns = sphere.spherePhong[sphere.bestSphere];
                    for(int i = 0; i < 3; i++){Ka[i] = sphere.sphereAmbient[sphere.bestSphere][i];}
                    for(int i = 0; i < 3; i++){Kd[i] = sphere.sphereDiffuse[sphere.bestSphere][i];}
                    for(int i = 0; i < 3; i++){Ks[i] = sphere.sphereSpecular[sphere.bestSphere][i];}
                    for(int i = 0; i < 3; i++){Kr[i] = sphere.sphereAttenuation[sphere.bestSphere][i];}
                    for(int i = 0; i < 3; i++){Ko[i] = sphere.sphereOpacity[sphere.bestSphere][i];}
                    bestSphere = sphere.bestSphere;
                } else {
                    sn_min = model.bestSurfaceNormal;
                    int materialIndex = (int) model.triangles[model.bestTriangle][3];
                    Ns = model.Ns[materialIndex];
                    for(int i = 0; i < 3; i++){Ka[i] = model.Ka[materialIndex][i];}
                    for(int i = 0; i < 3; i++){Kd[i] = model.Kd[materialIndex][i];}
                    for(int i = 0; i < 3; i++){Ks[i] = model.Ks[materialIndex][i];}
                    for(int i = 0; i < 3; i++){Kr[i] = 1.0;}
                    for(int i = 0; i < 3; i++){Ko[i] = 1.0;}
                    //bestModel = model.bestModel;
                }
            
                Vector3d Hit = pixelPt + shoot * t;
                colorTemp << Ka[0] * ambient[0], Ka[1] * ambient[1], Ka[2] * ambient[2];
                
                for(unsigned int i = 0; i< diffuse.size(); i++){
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "Run trace insideloop one time!\n";
            
                    // For each light source
                    Vector3d lightPostion;
                    if(diffuse[i][3] == 0){
                        lightPostion << diffuse[i][0]*1000000000, diffuse[i][1]*1000000000, diffuse[i][2]*1000000000;  
                    } else {
                        lightPostion << diffuse[i][0], diffuse[i][1], diffuse[i][2];
                    }
                    // From the hit position on the sphere to the light
                    Vector3d L = lightPostion - Hit; 
                    double distance = L.norm();
                    L.normalize();
                    
                    ///////////////////////////////////////////////////////////////////////////////
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "into shadow!\n";
                    //cout << "distance: " << distance << "\n";
                    //cout << "so what is L: " << L[0] << ", "
                    //    << L[1] << ", " << L[2] << "\n";
                    //cout << "so where is Hit: " << Hit[0] << ", "
                    //    << Hit[1] << ", " << Hit[2] << "\n";
                
                    bool shadow = false;
                    sphere.raySphereIntersection(Hit, L);
                    model.rayTriangleIntersection(Hit, L);
                    if(((sphere.flag == true) && (sphere.best < distance)) || ((model.flag == true) && (model.best < distance))) 
                        shadow = true;
                    
                    ///////////////////////////////////////////////////////////////////////////////
                    
                    if(sn_min.dot(L) > 0.0 && shadow == false){
                        Vector3d KB;
                        KB << Kd[0] * diffuse[i][4], Kd[1] * diffuse[i][5], Kd[2] * diffuse[i][6];
                        colorTemp += KB * (sn_min.dot(L));
                        
                        Vector3d toCamera = pixelPt - Hit; 
                        toCamera.normalize();
                        Vector3d reflection = (2 * sn_min.dot(L) * sn_min) - L;
                        double CdR = toCamera.dot(reflection);
                        Vector3d KS;
                        KS << Ks[0] * diffuse[i][4], Ks[1] * diffuse[i][5], Ks[2] * diffuse[i][6];
                        if(CdR > 0.0001){
                            colorTemp += KS * pow(CdR, Ns);
                        }
                        
                    }
                }
                
                for(int y = 0; y < 3; y++){
                    accumulation[y] += reflectionAttenuation[y] * Ko[y] * colorTemp(y);
                }
                
                //test code/////////////////////////////////////////////////////////////////////
                //cout << "colorTemp: " << colorTemp[0] << ", "
                //     << colorTemp[1] << ", " << colorTemp[2] << "\n";
                //cout << "accum after illumination: " << accumulation[0] << ", "
                //     << accumulation[1] << ", " << accumulation[2] << "\n";
                //cout << "So now we test for Ko: " << Ko[0] << ", "
                //     << Ko[1] << ", " << Ko[2] << "\n";
                //cout << "So now we test for reflectionAttenuation: " << reflectionAttenuation[0] << ", "
                //     << reflectionAttenuation[1] << ", " << reflectionAttenuation[2] << "\n";
                /////////////////////////////////////////////////////////////////////////
                double newRefAtt[3];
                for(int y = 0; y < 3; y++){
                        newRefAtt[y] = Kr[y] * reflectionAttenuation[y]; 
                }
                    
                if(levelParameter > 0){
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "first recursion\n";
                    
                    Vector3d reflectColor;
                    reflectColor << 0.0, 0.0, 0.0;
                    pixelPt << Hit[0], Hit[1], Hit[2];
                    Vector3d inverse = -1 * shoot;
                    Vector3d reflectionDirection = (2 * sn_min.dot(inverse) * sn_min) - inverse;
                    illumination(pixelPt, reflectionDirection, reflectColor, model, sphere, newRefAtt, levelParameter-1);
                    for(int y = 0; y < 3; y++){
                        accumulation[y] += reflectionAttenuation[y] * Ko[y] * reflectColor[y];
                    }
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "accum in first recursion: " << accumulation[0] << ", "
                    //    << accumulation[1] << ", " << accumulation[2] << "\n";
                }
                
                // Double Recursion
                // Check the tansparency of this sphere
                //test code/////////////////////////////////////////////////////////////////////
                //cout << "isSphere? " << isSphere << "\n";
                //cout << "Again: isSphere? " << sphere.flag << "\n";
                //cout << "AgainAgain: isModel? Impossible!" << model.flag << "\n";
                
                bool transperant = ((Ko[0] + Ko[1] + Ko[2]) < 3.0);
                if(levelParameter > 0 && transperant){
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "second recursion\n";
                    
                    Vector3d thru;
                    thru << 0.0, 0.0, 0.0;
                    //test code/////////////////////////////////////////////////////////////////////
                    //cout << "HIT: " << Hit[0] << ", "
                    //    << Hit[1] << ", " << Hit[2] << "\n";
                        
                    Vector3d outRay = sphere.refractionExit(bestSphere, -1 * shoot, Hit, eta_outside);
                    bool cannot = ((outRay[0] + outRay[1] + outRay[2]) == 0.0);
                    //cout << "out ray: " << outRay[0] << ", "
                    //    << outRay[1] << ", " << outRay[2] << "\n";
                    //cout << "exitPoint: " << sphere.exitPoint[0] << ", "<< sphere.exitPoint[1] << ", " << sphere.exitPoint[2] << "\n";
                    //cout << "exitDirection: " << sphere.exitDirection[0] << ", "<< sphere.exitDirection[1] << ", " << sphere.exitDirection[2] << "\n";
                    
                    if(!cannot){
                        //test code/////////////////////////////////////////////////////////////////////
                        //cout << "ray L: " << sphere.exitPoint[0] << ", "
                        //     << sphere.exitPoint[1] << ", " << sphere.exitPoint[2] << "\n";
                        //cout << "ray D: " << sphere.exitDirection[0] << ", "
                        //     << sphere.exitDirection[1] << ", " << sphere.exitDirection[2] << "\n";
                        illumination(sphere.exitPoint, sphere.exitDirection, thru, model, sphere, newRefAtt, levelParameter-1);
                        //test code/////////////////////////////////////////////////////////////////////
                        //cout << "after thru: " << thru[0] << ", "
                        //     << thru[1] << ", " << thru[2] << "\n";
                        //cout << "refatt: " << reflectionAttenuation[0] << ", "
                        //     << reflectionAttenuation[1] << ", " << reflectionAttenuation[2] << "\n";
                        //cout << "Ko in second recursion: " << Ko[0] << ", "
                        //    << Ko[1] << ", " << Ko[2] << "\n";     
                        
                        for(int y = 0; y < 3; y++){
                            accumulation[y] += reflectionAttenuation[y] * (1.0 - Ko[y]) * thru[y];
                        }
                        //test code/////////////////////////////////////////////////////////////////////
                        //cout << "accum in second recursion: " << accumulation[0] << ", "
                        //      << accumulation[1] << ", " << accumulation[2] << "\n";
                        }
                    }
            }
}   

/*
 * This method checks if the ray intersects either a sphere or a model.
 * By the way, get the surface normal for that closest sphere or triangle.
 */    
void Light::rayFind(Vector3d &pixelPt, Vector3d &shoot, Model& model, Sphere& sphere){
    //test code/////////////////////////////////////////////////////////////////////
    //cout << "Why do I run test? in find!\n";
    
    model.rayTriangleIntersection(pixelPt, shoot);
    sphere.raySphereIntersection(pixelPt, shoot);
                    
    t = numeric_limits<double>::max(); // The closest distance.
    flag = false; // The ray either intersects a shpere or a model.
    isSphere = false; // The ray intersects a shpere.
            
    if((model.flag == true) && (sphere.flag == true)){
        //test code/////////////////////////////////////////////////////////////////////
        //cout << "First!\n";
        
        flag = true;
        if(sphere.best < model.best){
            t = sphere.best;
            isSphere = true;
        } else {
            t = model.best;
        }
    } else if(sphere.flag == true){
        //test code/////////////////////////////////////////////////////////////////////
        //cout << "Second!\n";
    
        flag = true;
        isSphere = true;
        t = sphere.best;
    } else if(model.flag == true){
        //test code/////////////////////////////////////////////////////////////////////
        //cout << "Third!\n";
    
        flag = true;
        t = model.best;
    }
    
    //test code/////////////////////////////////////////////////////////////////////
    //if(flag == false)
    //cout << "I didn't find!\n";
}

