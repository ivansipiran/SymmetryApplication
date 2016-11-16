#include <iostream>
#include <cmath>
#include "Util/rotation.h"

using namespace std;

float Rotation::rotationMatrix[4][4];
float Rotation::inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
float Rotation::outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};

void Rotation::multiplyMatrix()
{
    for(int i = 0; i < 4; i++ ){
        for(int j = 0; j < 1; j++){
            outputMatrix[i][j] = 0;
            for(int k = 0; k < 4; k++){
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
}

void Rotation::setUpRotationMatrix(float angle, float u, float v, float w)
{
    float L = (u*u + v * v + w * w);
    angle = angle * 3.141592 / 180.0; //converting to radian value
    float u2 = u * u;
    float v2 = v * v;
    float w2 = w * w;

    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][3] = 0.0;

    rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][3] = 0.0;

    rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix[2][3] = 0.0;

    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;
}

float * Rotation::rotateVertex(float x, float y, float z, float angle, float u, float v, float w){
    //Setup input matrix
    inputMatrix[0][0] = x;
    inputMatrix[1][0] = y;
    inputMatrix[2][0] = z;
    inputMatrix[3][0] = 1.0;

    //Setup rotation matrix
    setUpRotationMatrix(angle, u, v, w);

    //Calculate
    multiplyMatrix();

    //Output result
    float * result = new float[3];
    result[0]=outputMatrix[0][0];
    result[1]=outputMatrix[1][0];
    result[2]=outputMatrix[2][0];
    return result;
}
