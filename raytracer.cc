#include "Model.h"
#include "Light.h"
#include "Camera.h"
#include "Sphere.h"

using namespace std;
using namespace Eigen;

/*
 * This program takes two command line arguments.
 * The first argument is a driver file. 
 * The second argument is name of the image your program will write back to disk 
 * in the same folder where your executable is located.
 */ 
int main(int argc, char **argv){
    
    /*
     * This code segment checks the number of arguments this program takes.
     */
    if(argc != 3){
        cerr << "Usage: ./raytracer driver00.txt driver00.ppm.\n";
        return 1;
    }
    
    /*
     * This code segment checks if the driver file can be read successfully.
     */
    ifstream driver(argv[1]);
    
    if (!driver){
        cerr << "Fail to open the driver file: " << argv[1] << "\n"; 
        return 1;
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    
    /*
     * This code segment reads the .obj models from the driver.
     */
    Model model;
    if(model.readDriver(driver) == 1) return 1;
    
    //test code/////////////////////////////////////////////////////////////////////
    //for(unsigned int i = 0; i < model.triangles.size(); i++){
    //    cout << "test material indices: " << model.triangles[i][0] << " " << model.triangles[i][1] << " "
    //        << model.triangles[i][2] << " " << model.triangles[i][3] << "\n";
    //}
    //for(unsigned int i = 0; i < model.materialssss.size(); i++){
    //    cout << "test material contents: " << model.Kd[i][0] << " " << model.Kd[i][1] << " "
    //        << model.Kd[i][2] << "\n";
    //}
    
    
    /*
     * This code segment reads the spheres from the driver.
     */
    driver.open(argv[1]);
    Sphere sphere;
    sphere.readDriver(driver);
    
    /*
     * This code segment reads the camera specification from the driver.
     */
    driver.open(argv[1]);
    Camera camera;
    camera.readDriver(driver);
    camera.placement();
    
    /*
     * This code segment reads the light sources from the driver.
     */
    driver.open(argv[1]);
    Light light;
    light.readDriver(driver);

    ///////////////////////////////////////////////////////////////////////////////
    
    ofstream ppm(argv[2]);
    if (!ppm){
            cerr << "Fail to open the output file!\n"; 
            return 1;
        }
    ppm << "P3\n";
    ppm << camera.res[0] << " " << camera.res[1] << " " << 255 << "\n";
    
    ///////////////////////////////////////////////////////////////////////////////

    for(int j = 0; j < camera.res[1]; j++){
        for(int i = 0; i < camera.res[0]; i++){
            //test code/////////////////////////////////////////////////////////////////////
            //int i = 9; int j = 9;
            //cout << "(i,j): " << i << " " << j << "\n";
            
            double px = (double)i/(camera.res[0]-1) * (camera.bounds[2] - camera.bounds[0]) + camera.bounds[0];    
            double py = (double)j/(camera.res[1]-1) * (camera.bounds[1] - camera.bounds[3]) + camera.bounds[3];    
            Vector3d pixelPt,shoot;
            pixelPt = camera.E + ((-camera.d) * camera.W) + (px * camera.U) + (py * camera.V);
            shoot = pixelPt - camera.E; 
            shoot.normalize();
        
            // For each pixel: ray trace
            // Initialize accumulation color [accum] in SageMath
            light.color << 0,0,0; 
            
            // Initialize reflection Attenuation [refatt] in SageMath
            double reflectionAttenuation[3];
            for(int i = 0; i<3; i++){reflectionAttenuation[i] = 1;} 
            
            light.illumination(pixelPt, shoot, light.color, model, sphere, reflectionAttenuation, light.depth); 
            
            //test code/////////////////////////////////////////////////////////////////////
            //cout << "color: " << light.color[0] << ", "
            //    << light.color[1] << ", " << light.color[2] << "\n";
                
            // map the color from "0 to 1" to "0 to 255"
            for(int y = 0; y < 3; y++){
                light.color(y) = (int)floor(light.color(y) * 255);
                if(light.color(y) > 255){
                    light.color(y) = 255;
                } else if (light.color(y) < 0){
                    light.color(y) = 0;
                }
            } 
            
            ppm << (int)light.color(0) << " " << (int)light.color(1) << " " << (int)light.color(2) << " ";  
            
        }
    }
    ppm.close();
    
    return 0;
}
