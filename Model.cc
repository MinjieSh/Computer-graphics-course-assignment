#include "Model.h"

using namespace std;
using namespace Eigen;

/*
 * This method is to read the transformation and shading information of each .obj model from the driver file.
 * 
 * The general format for each model is as follows: 
 * "model wx wy wz theta scale tx ty tz option model.obj"
 * 
 * The triplet "wx wy wz" represents the axis about which to rotate 
 * followed by "theta" in degrees representing the angle by which to rotate (thus Axis-Angle format). 
 * The "scale" is a uniform scaling factor to shrink or grow a model.
 * The triplet "tx ty tz" defines a model-to-world translation.
 * The "option" represents two shading options: sharp and smooth (Gouraud).
 * "model.obj" is the model file on which you will be applying the defined 3D transformation.
 * 
 * The order of transformation is rotate, scale, and then translate.
 * 
 * There would be multiple model lines in a driver file.
 */ 
int Model::readDriver(ifstream& driver){
    
    record.push_back(0);
    counter = 0;
    size = 0;
    
    string line = "";
    while(getline(driver, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id == "model"){
                // For each model:
        
                /*
                 * Axis Angle Rotation: wx wy wz theta. 
                 */
                array<double,4> axis;
                for(int i = 0; i < 4; i++){iss >> axis[i];}
                axisAngleRotation.push_back(axis);
                
                /*
                 * scaling factor: scale. 
                 */
                double sca;
                iss >> sca;
                scalingFactor.push_back(sca);
                
                /*
                 * Translation: tx ty tz. 
                 */
                array<double,3> tra;
                for(int i = 0; i < 3; i++){iss >> tra[i];}
                translation.push_back(tra);
                
                /*
                 * Option: sharp or smooth. 
                 */
                string option;
                iss >> option;
                if(option == "smooth"){
                    smooth.push_back(true);
                } else if (option == "sharp"){
                    smooth.push_back(false);
                }
                
                /*
                 * Model file: model.obj. 
                 */
                string mod;
                iss >> mod;
                modelName.push_back(mod);
                
                /*
                 * This code segment checks if the model file can be read successfully.
                 */
                model.open(modelName[size]);
                if(!model) {cerr << "Fail to open the .obj file: " + modelName[size] + "\n"; return 1;}
                if(readModel(model) == 1) return 1;
                model.close();
                
                /*
                 * The number of models.
                 */
                size ++; 
                
                /*
                 * Record the number of vertices for each model.
                 */
                record.push_back(counter);
            } else { continue; }
        } else { continue; }
    }
    
    driver.close(); 
    return 0;
}


/*
 * This method is to read each .obj model from the model file.
 * 
 * The general format for each model is as follows: 
 * mtllib [external .mtl file name]
 * v v1 v2 v3 ...
 * f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
 * "mtllib" represents the model file.
 * "v" represents vertices.
 * "f" represents faces.
 * The format for the faces is: vertex index/vertex texture coordinate index/vertex normal index
 * Or without the vertex texture: vertex index//vertex normal index
 */ 
int Model::readModel(ifstream& model){

    /*
     * Get transformation Matrix.
     */
    transformationMatrix(axisAngleRotation[size], scalingFactor[size], translation[size]);
    
    string line = "";
    while(getline(model, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id == "v"){
            
                /*
                 * Read all the vertices.
                 */
                array<double,3> before;
                for(int i = 0; i < 3; i++){iss >> before[i];}
            
                /*
                 * Do transformation for each vertex.
                 */
                transform(before);
                vertices.push_back(before);
                
                /*
                 * The number of vertices.
                 */
                counter++; 
                
                } else if(id == "usemtl"){
                  
                    string materialName;
                    iss >> materialName;
                    vector<string>::iterator it = find(materialssss.begin(), materialssss.end(), materialName);
                    materialId = distance(materialssss.begin(), it);
                    
                } else if(id == "f"){
            
                /*
                 * Read all the faces.
                 */
                array<int, 4> faces; 
                for(int i = 0; i < 3; i++){
                    /*
                     * Read the vertex indices of a face.
                     */
                    string str;
                    iss >> str; 
                    istringstream input(str);
                    input >> faces[i];
                    /*
                     * Convert the indices.
                     * Combine all the vertices into a single model.
                     * For example, 0,1,2,3; 0,1,2,3,4,5; 0,1,2
                     * After the conversion: 0,1,2,3;
                     * 4 (0+4),5 (1+4),6 (2+4),7 (3+4),8 (4+4),9 (5+4);
                     * 10 (0+10),11 (1+10),12 (2+10) 
                     * Pay attention: the vertices starting from 1 in faces!
                     */
                    faces[i] += record[size]-1; 
                }
                faces[3] = materialId;
                triangles.push_back(faces);
                
            } else if(id == "mtllib"){
            
                /*
                 * Read the material file.
                 */
                string mat;
                iss >> mat;
                materialFileName.push_back(mat);
            
                /*
                 * This code segment checks if the material file can be read successfully.
                 */
                material.open(materialFileName[size]);
                if (!material){cerr << "Fail to open the material file: " + materialFileName[size] + "\n"; return 1;}
                readMaterial(material);  
                material.close();
            } else {continue;} 
        } else {continue;}
    }
    return 0;
}    
    
    
/*
 * This method is to read each model's material information from the material file.
 * 
 * The general format for the material information is as follows: 
 * Ns represents the phong constant.
 * Ka represents (Ka_red, Ka_green, Ka_blue) ambient coefficients.
 * Kd represents diffuse (Kd_red, Kd_green, Kd_blue) coefficients. 
 * Ks represents (Ks_red, Ks_green, Ks_blue) specular coefficients.
 * 
 */     
void Model::readMaterial(ifstream& model){
    string line = "";
    while(getline(material, line)){
        istringstream iss(line);
        string id = "";
        if (iss >> id){
            if(id == "newmtl"){
                
                /*
                 * Read the material name.
                 */
                string materialName;
                iss >> materialName;
                materialssss.push_back(materialName);
                
            }else if (id == "Ns"){
                
                /*
                 * Read the phong constant.
                 */
                double pho;
                iss >> pho;
                Ns.push_back(pho);
                
            } else if(id == "Ka"){
                
                /*
                 * Read the ambient coefficients.
                 */
                array<double,3> amb;
                for(int i = 0; i < 3; i++){iss >> amb[i];}
                Ka.push_back(amb);
                
            } else if(id == "Kd"){
                
                /*
                 * Read the diffuse coefficients.
                 */
                array<double,3> dif;
                for(int i = 0; i < 3; i++){iss >> dif[i];}
                Kd.push_back(dif);
                
            } else if(id == "Ks"){
                
                /*
                 * Read the specular coefficients.
                 */
                array<double,3> spe;
                for(int i = 0; i < 3; i++){iss >> spe[i];}
                Ks.push_back(spe);
                
            } else {continue;}
        } else {continue;}
    }
}    
    
    
    
/*
 * This method is to get the transformation matrix.
 * The order of transformation is rotate, scale, and then translate.
 */      
void Model::transformationMatrix(array<double,4> axisAngleRotation, double scalingFactor, array<double,3> translation){
       
       /*
        * Axis Angle Rotation
        */
        Vector3d w;
        w << axisAngleRotation[0], axisAngleRotation[1], axisAngleRotation[2];
        w.normalize();
        
        Vector3d m;
        m << abs(axisAngleRotation[0]), abs(axisAngleRotation[1]), abs(axisAngleRotation[2]);
        int j = -1;
        double minOfW = m(0);
        for (int i = 0; i < 3; i++){
            if(m(i) <= minOfW)
                j = i;
        }
        m(j) = 1.0;
        
        Vector3d u;
        u = w.cross(m); 
        u.normalize();
        
        Vector3d v;
        v = w.cross(u);
        
        Matrix4d r;
        r << u(0), u(1), u(2), 0,
             v(0), v(1), v(2), 0,
             w(0), w(1), w(2), 0,
             0, 0, 0, 1;
             
        Matrix4d rt;
        rt = r.transpose();
        
        Matrix4d rz;
        double theta = (axisAngleRotation[3] / 180) * 3.14159265359;
        double cosOfTheta = cos(theta);
        double sinOfTheta = sin(theta);
        rz << cosOfTheta, -sinOfTheta, 0, 0,
              sinOfTheta, cosOfTheta, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;
              
        Matrix4d rotationMatrix;
        rotationMatrix = rt * rz * r;
        
        /*
         * Scale
         */
        Matrix4d scaleMatrix;
        scaleMatrix << scalingFactor, 0, 0, 0,
                       0, scalingFactor, 0, 0,
                       0, 0, scalingFactor, 0,
                       0, 0, 0, 1;
        
        /*
         * Translation
         */
        Matrix4d translationMatrix;
        translationMatrix << 1, 0, 0, translation[0],
                             0, 1, 0, translation[1],
                             0, 0, 1, translation[2],
                             0, 0, 0, 1;
         
        /*
         * Get transformation matrix
         */                     
        Matrix = translationMatrix * scaleMatrix * rotationMatrix;
}    

/*
 * This method is to transform each vertex.
 */   
void Model::transform(array<double,3> &arr){
    
    MatrixXd pts(4,1), after(4,1);
            
    pts(0,0) = arr[0];
    pts(1,0) = arr[1];
    pts(2,0) = arr[2];
    pts(3,0) = 1;
            
    after = Matrix * pts;
    
    arr[0] = (double)after(0,0);
    arr[1] = (double)after(1,0);
    arr[2] = (double)after(2,0);
}
    
/*
 * This method does Ray-Triangle-Intersection and get the minimum distance from
 * the starting location of the ray to all the triangle.
 * By the way, get the index of the closest model, 
 * and the surface normal of the closest triangle.
 * 
 * The two parameters are the starting location and the direction of the ray. 
 */    
void Model::rayTriangleIntersection(Vector3d pixelPt, Vector3d shoot){
    
        flag = false;
        best = numeric_limits<double>::max();
        double beta[triangles.size()], gamma[triangles.size()], t[triangles.size()];
        Vector3d sn[triangles.size()];
        bool intersect[triangles.size()];
        
        for(unsigned int k = 0; k < triangles.size(); k++){
            // For each triangle
            Vector3d A, B, C;
            
            int indexA = triangles[k][0];
            array<double,3> vertex = vertices[indexA];
            A << vertex[0], vertex[1], vertex[2];
            
            int indexB = triangles[k][1];
            vertex = vertices[indexB];
            B << vertex[0], vertex[1], vertex[2];
            
            int indexC = triangles[k][2];
            vertex = vertices[indexC];
            C << vertex[0], vertex[1], vertex[2];
           
            Vector3d col1, col2, col3;
            col1 = A - B;
            col2 = A - C;
            col3 = A - pixelPt;
            
            sn[k] = col1.cross(col2); 
            sn[k].normalize(); // Calculate the true normal.
        
            intersect[k] = false;
            beta[k] = -1, gamma[k] = -1, t[k] = -1;
            
            Matrix3d M, M1, M2, M3;
            M << col1(0), col2(0), shoot(0),
                 col1(1), col2(1), shoot(1),
                 col1(2), col2(2), shoot(2);
                
            double deterOfM = M.determinant();
            if(abs(deterOfM) < 0.00001){
                continue;
            }     
                     
            M1 << col3(0), col2(0), shoot(0),
                  col3(1), col2(1), shoot(1),
                  col3(2), col2(2), shoot(2); 
            if(deterOfM != 0){
                beta[k] = M1.determinant() / deterOfM;
            }      
            if(!(beta[k] >= 0)) continue;
            
            M2 << col1(0), col3(0), shoot(0),
                  col1(1), col3(1), shoot(1),
                  col1(2), col3(2), shoot(2);
            if(deterOfM != 0){
                gamma[k] = M2.determinant() / deterOfM;
            }      
            if(!(gamma[k] >= 0)) continue;
            if(!(beta[k] + gamma[k] <= 1)) continue;
                
            M3 << col1(0), col2(0), col3(0),
                  col1(1), col2(1), col3(1),
                  col1(2), col2(2), col3(2);
            if(deterOfM != 0){
                t[k] = M3.determinant() / deterOfM;
            } 
            if(!(t[k] > 0.00001)) continue;
            
            if((t[k] > 0.00001) && (beta[k] >= 0) && (gamma[k] >= 0) && (beta[k] + gamma[k] <= 1)){
                intersect[k] = true;
                flag = true;
            }
            
        }
            
        /* 
         * Find the closest triangle.
         * By the way, get the index of the closest model, 
         * the minimum distance, 
         * the surface normal of the closest triangle.
         */
        int minIndex = -1;
        if(flag == true){
            minIndex = -1;
            best = numeric_limits<double>::max();
            for(unsigned int k = 0; k < triangles.size(); k++){
                if(intersect[k] == true){
                    if(t[k] < best){
                        best = t[k]; 
                        minIndex = k;
                    }
                }
            }
         
            /*
             * For example, the record is 0,4,10,13 (size=3)
             * the model is 0,1,2
             * 0,1,2,3; 0,1,2,3,4,5; 0,1,2
             * 0,1,2,3; 4,5,6,7,8,9; 10,11,12
             * if the minIndex = 3, index of model is 0
             * if the minIndex = 5, index of model is 1
             * if the minIndex = 10, index of model is 2
             */
            bestTriangle = minIndex;
            for(int i = 0; i < size; i ++){
                if(minIndex >= record[i] && minIndex < record[i+1])
                    bestModel = i;
            }
            
        
            if(smooth[bestModel]){
                Vector3d N[3];
                double threshold = cos((22.5 / 180) * 3.14159265359);
                for(int i = 0 ; i < 3; i++){
                    int count = 1;
                    N[i] << sn[minIndex][0], sn[minIndex][1], sn[minIndex][2];
                        for(int j = 0; j < (int)triangles.size(); j++){
                            //cout <<"face:" <<  minIndex << "go through" << j << "\n";
                            if(j == minIndex) continue;
                            if(triangles[j][0] == triangles[minIndex][i] 
                                || triangles[j][1] == triangles[minIndex][i] 
                                || triangles[j][2] == triangles[minIndex][i]){
                                //cout << "hit" << j << "\n";
                                if(sn[minIndex].dot(sn[j]) >= threshold){
                                    N[i] += sn[j];
                                    count ++;
                                }
                            }
                        }
                        N[i] = N[i] / (double)count;
                        N[i].normalize();
                    }
                    
                bestSurfaceNormal = (1- beta[minIndex] - gamma[minIndex]) * N[0] + beta[minIndex] * N[1] + gamma[minIndex] * N[2];
                bestSurfaceNormal.normalize();
                
            } else {
                if(shoot.dot(sn[minIndex]) > 0.0){
                    for(int x = 0; x < 3; x++){sn[minIndex](x) = -sn[minIndex](x);}    
                }
                bestSurfaceNormal = sn[minIndex];
            }
            
       }      
        
}    

