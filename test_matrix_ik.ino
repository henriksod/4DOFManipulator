/**********************************************************************************************
 * 3DOF Kinematics Solver
 * by Henrik Söderlund <henrik.a.soderlund@gmail.com>
 *
 * Copyright (c) 2018 Henrik Söderlund

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 **********************************************************************************************/

#include "MatrixMath.h"
#include <MemoryUsage.h>

int idx(int m, int j, int i=0) {
    return m * j + i;
}

template <typename T = float>
void Slice(T* A, int n, int a, int b, T* B) {
    int i, j;
    for (i = 0; i < a; i++)
        for(j = 0; j < b; j++) {
            B[b * i + j] = A[n * i + j];
        }
}

template <typename T = float>
void Eye4(T* M) { // TODO: Make generic for m x n   
    M[idx(4,0,0)] = 1; M[idx(4,0,1)] = 0; M[idx(4,0,2)] = 0; M[idx(4,0,3)] = 0;
    M[idx(4,1,0)] = 0; M[idx(4,1,1)] = 1; M[idx(4,1,2)] = 0; M[idx(4,1,3)] = 0;
    M[idx(4,2,0)] = 0; M[idx(4,2,1)] = 0; M[idx(4,2,2)] = 1; M[idx(4,2,3)] = 0;
    M[idx(4,3,0)] = 0; M[idx(4,3,1)] = 0; M[idx(4,3,2)] = 0; M[idx(4,3,3)] = 1;
}

template <typename T = float>
void ScrewZ(T theta, T d, T* M) {
    T c = cos(theta);
    T s = sin(theta);
    
    M[idx(4,0,0)] = c; M[idx(4,0,1)] = -s; M[idx(4,0,2)] = 0; M[idx(4,0,3)] = 0;
    M[idx(4,1,0)] = s; M[idx(4,1,1)] = c; M[idx(4,1,2)] = 0; M[idx(4,1,3)] = 0;
    M[idx(4,2,0)] = 0; M[idx(4,2,1)] = 0; M[idx(4,2,2)] = 1; M[idx(4,2,3)] = 0;
    M[idx(4,3,0)] = 0; M[idx(4,3,1)] = 0; M[idx(4,3,2)] = 0; M[idx(4,3,3)] = 1;
}

template <typename T = float>
void ScrewX(T alpha, T a, T* M) {
    T c = cos(alpha);
    T s = sin(alpha);

    M[idx(4,0,0)] = 1; M[idx(4,0,1)] = 0; M[idx(4,0,2)] = 0; M[idx(4,0,3)] = a;
    M[idx(4,1,0)] = 0; M[idx(4,1,1)] = c; M[idx(4,1,2)] = -s; M[idx(4,1,3)] = 0;
    M[idx(4,2,0)] = 0; M[idx(4,2,1)] = s; M[idx(4,2,2)] = c; M[idx(4,2,3)] = 0;
    M[idx(4,3,0)] = 0; M[idx(4,3,1)] = 0; M[idx(4,3,2)] = 0; M[idx(4,3,3)] = 1;
}

template <typename T = float>
T Norm(T* d, int n) {
    T diff_sq = 0;
    for(int j = 0; j < n; j++)
        diff_sq += d[j]*d[j];
    return sqrt(diff_sq);
}

template <typename T = float>
void forward_kinematics(T* angles, T* lengths, T* output, int joint_number = 3) {
    T* buffer4x4_1 = new T[idx(4,4)];
    T* buffer4x4_2 = new T[idx(4,4)];
    
    T* T00 = new T[idx(4,4)];
    Eye4(T00);
    
    T* T01 = new T[idx(4,4)];
    ScrewZ<T>(angles[0], 0, buffer4x4_1);
    Matrix.Multiply(T00, buffer4x4_1, 4, 4, 4, buffer4x4_2);
    ScrewX<T>(0, lengths[0], buffer4x4_1);
    Matrix.Multiply(buffer4x4_2, buffer4x4_1, 4, 4, 4, T01);

    if (joint_number == 1) {
      output[0] = T01[idx(4,0,3)];
      output[1] = T01[idx(4,1,3)];
      output[2] = T01[idx(4,2,3)];

      delete[] buffer4x4_1;
      delete[] buffer4x4_2;
      delete[] T00;
      delete[] T01;
      
      return;
    }
    
    T* T02 = T00; // Reuse pointer
    ScrewZ<T>(angles[1], 0, buffer4x4_1);
    Matrix.Multiply(T01, buffer4x4_1, 4, 4, 4, buffer4x4_2);
    ScrewX<T>(0, lengths[1], buffer4x4_1);
    Matrix.Multiply(buffer4x4_2, buffer4x4_1, 4, 4, 4, T02);

    if (joint_number == 2) {
      output[0] = T02[idx(4,0,3)];
      output[1] = T02[idx(4,1,3)];
      output[2] = T02[idx(4,2,3)];

      delete[] buffer4x4_1;
      delete[] buffer4x4_2;
      delete[] T00;
      delete[] T01;
      
      return;
    }
    
    T* T03 = T00; // Reuse pointer
    ScrewZ<T>(angles[2], 0, buffer4x4_1);
    Matrix.Multiply(T02, buffer4x4_1, 4, 4, 4, buffer4x4_2);
    ScrewX<T>(0, lengths[2], buffer4x4_1);
    Matrix.Multiply(buffer4x4_2, buffer4x4_1, 4, 4, 4, T03);

    output[0] = T03[idx(4,0,3)];
    output[1] = T03[idx(4,1,3)];
    output[2] = T03[idx(4,2,3)];
    
    delete[] buffer4x4_1;
    delete[] buffer4x4_2;
    delete[] T00;
    delete[] T01;
}

template <typename T = float>
uint8_t inverse_kinematics(
    T x,
    T y,
    T z,
    T tool_angle,
    T tool_radius,
    T* lengths,
    T* output
) {
    T* buffer4x4 = new T[idx(4,4)];
    T* buffer4x4_2 = new T[idx(4,4)];
    T* buffer3x3 = new T[idx(3,3)];
    T* buffer3 = new T[3];
    
    T* o = new T[3];
    o[0] = x; o[1] = y; o[2] = z;
    
    T* Tnorm = new T[idx(4,4)];
    T* Rnorm = new T[idx(3,3)];
    ScrewZ<T>(0 /* Surface normal angle */, 0, Tnorm);
    Slice<T>(Tnorm, 4, 3, 3, Rnorm);
    
    // Compute rotation matrix of tool angle towards target
    ScrewZ<T>(tool_angle, 0, buffer4x4);
    T* Rtool = new T[idx(3,3)];
    Matrix.Multiply(Tnorm, buffer4x4, 4, 4, 4, buffer4x4_2);
    Slice<T>(buffer4x4_2, 4, 3, 3, Rtool);
    
    // Calculate the wrist center position from the end effector tool
    // position, offset with an angle towards the target.
    T* oc = new T[3];
    oc[0] = -1; oc[1] = 0; oc[2] = 0;
    //Matrix.Scale(oc, 3, 1, tool_radius);
    Matrix.Scale(oc, 3, 1, lengths[2]);
    Matrix.Multiply(Rtool, oc, 3, 3, 1, buffer3);
    Matrix.Copy(buffer3, 3, 1, oc);
    Matrix.Add(oc, o, 3, 1, buffer3);
    Matrix.Copy(buffer3, 3, 1, oc);
    buffer3[0] = -1; buffer3[1] = 0; buffer3[2] = 0;
    Matrix.Scale(buffer3, 3, 1, tool_radius);
    Matrix.Multiply(Rnorm, buffer3, 3, 3, 1, o); // Reuse o pointer
    Matrix.Subtract(oc, o, 3, 1, buffer3);
    Matrix.Copy(buffer3, 3, 1, oc);
    
    // Compute distance from frame 0 to wrist center position.
    T dist = Norm<T>(oc, 3);
    
    // Compute the trigonometric equation that describes the angle of
    // the elbow joint based on the location of the wrist center using
    // law of cosines.
    T c_phi = (dist*dist - lengths[0]*lengths[0] - lengths[1]*lengths[1]) / (2*lengths[0]*lengths[1]);
    T s_phi = sqrt(1 - c_phi*c_phi);
    
    // Assign the angles for shoulder and elbow.
    T a1 = atan2(oc[1], oc[0]) + atan2(lengths[1]*s_phi, lengths[0]+lengths[1]*c_phi);
    T a2 = -atan2(s_phi, c_phi);
    
    // Extract the rotation matrix that describes how to rotate from
    // frame 2 to frame 3 based on the now known rotation matrix R02.
    ScrewZ<T>(a2, 0, buffer4x4);
    T* R02 = new T[idx(4,4)];
    Matrix.Copy(buffer4x4, 4, 4, R02);
    ScrewZ<T>(a1, 0, buffer4x4);
    Matrix.Multiply(buffer4x4, R02, 4, 4, 4, buffer4x4_2);
    Matrix.Copy(buffer4x4_2, 4, 4, R02);
    Slice<T>(R02, 4, 3, 3, buffer3x3);
    uint8_t solution = Matrix.Invert(buffer3x3, 3);
    T* R23 = Rtool; // Reuse pointer
    Matrix.Multiply(buffer3x3, Rnorm, 3, 3, 3, R23);
    
    if (solution == true) {
        output[0] = a1;
        output[1] = a2;
        output[2] = atan2(R23[idx(3,1,0)], R23[idx(3,0,0)]) + tool_angle;
    }
    
    delete[] buffer4x4;
    delete[] buffer4x4_2;
    delete[] buffer3x3;
    delete[] buffer3;
    delete[] Tnorm;
    delete[] Rnorm;
    delete[] o;
    delete[] Rtool;
    delete[] oc;
    delete[] R02;
    
    return solution;
}


#define NUM_JOINTS 4

#define STEP_ITERATIONS 10
#define MIN_STEP_SIZE 1
#define MAX_STEP_SIZE 50
#define MIN_X 150
#define MAX_X 200

void benchmark() {

  long n = 0;
  float mean_time = 0;
  float variance_time = 0;
  unsigned long min_time = 1000000;
  unsigned long max_time = 0;
  unsigned long start_time = 0;
  unsigned long end_time = 0;
  float delta_time = 0;
  float sum_xyz_error = 0;
  unsigned long sum_num_iterations = 0;
  
  mtx_type angles[] = {0, 0, 0};
  mtx_type lengths[] = {200, 200, 270};
  mtx_type end_position[] = {0, 0, 0};
  

  // Config
  Serial.print("Iteration");
  Serial.print("\t");
  Serial.print("Num joints");
  Serial.print("\t");
  Serial.print("Tolerance [mm]");
  Serial.print("\t");
  Serial.print("Move distance [mm]");
  Serial.print("\t");
  // Results
  Serial.print("Solvable");
  Serial.print("\t");
  Serial.print("Mean iterations");
  Serial.print("\t");
  Serial.print("Mean xyz error [mm]");
  Serial.print("\t");
  Serial.print("Mean exec time [us]");
  Serial.print("\t");
  Serial.print("Variance exec time [us]");
  Serial.print("\t");
  Serial.print("Min exec time [us]");
  Serial.print("\t");
  Serial.print("Max exec time [us]");
  Serial.print("\t");
  Serial.print("FPS");
  Serial.println();
  
  for (int step_i = STEP_ITERATIONS; step_i > 0; step_i--) {
      float step_size = ((float)MAX_STEP_SIZE-MIN_STEP_SIZE)/step_i;
  
      uint8_t result_status = true;
      float xyz = MIN_X;
      for (xyz = MIN_X; xyz < MAX_X; xyz+=step_size) {
        
        start_time = micros();
        result_status = result_status && inverse_kinematics<mtx_type>(xyz, xyz/10, 0, -HALF_PI, 0, lengths, angles);
        end_time = micros();

        forward_kinematics<mtx_type>(angles, lengths, end_position);

        float solved_x = end_position[0];
        float solved_y = end_position[1];
        float solved_z = end_position[2];
  
        n += 1;
        delta_time = end_time - start_time;
        float prev_mean_time = mean_time;
        mean_time += (delta_time - prev_mean_time) / n;
        variance_time += (delta_time - prev_mean_time)*(delta_time - mean_time);
        float diff_x = xyz-solved_x;
        float diff_y = xyz/10-solved_y;
        float diff_z = 0-solved_z;
        
        sum_xyz_error += sqrt(diff_x*diff_x +  diff_y*diff_y + diff_z*diff_z);
        sum_num_iterations += 1;
        
        if (min_time > delta_time) min_time = delta_time;
        if (max_time < delta_time) max_time = delta_time;
      }
      
      // Config
      Serial.print(STEP_ITERATIONS-step_i);
      Serial.print("\t");
      Serial.print(NUM_JOINTS);
      Serial.print("\t");
      Serial.print(0);
      Serial.print("\t");
      Serial.print(sqrt(3) * step_size, 2);
      Serial.print("\t");
      // Results
      Serial.print(result_status);
      Serial.print("\t");
      Serial.print((float)sum_num_iterations / n, 2);
      Serial.print("\t");
      Serial.print((float)sum_xyz_error / n, 2);
      Serial.print("\t");
      Serial.print(mean_time);
      Serial.print("\t");
      Serial.print(sqrt(variance_time / (n-1)));
      Serial.print("\t");
      Serial.print(min_time);
      Serial.print("\t");
      Serial.print(max_time);
      Serial.print("\t");
      Serial.print((1 / (mean_time/1000000)));
      Serial.println();
  
      n = 0;
      mean_time = 0;
      variance_time = 0;
      min_time = 100000;
      max_time = 0;
      sum_xyz_error = 0;
      sum_num_iterations = 0;
    }  
}

void setup() {
    Serial.begin(115200);

    MEMORY_PRINT_TOTALSIZE;
    STACKPAINT_PRINT;
    FREERAM_PRINT;
  
    mtx_type angles[] = {0, 0, 0};
    mtx_type end_position[] = {0, 0, 0};

    // Define arm
    mtx_type x = 150;
    mtx_type y = 50;
    mtx_type z = 0;
    mtx_type tool_angle = -HALF_PI;
    mtx_type tool_radius = 0;
    mtx_type lengths[] = {200, 200, 270};
    
    // Solve IK
    unsigned long prev_time = millis();
    inverse_kinematics<mtx_type>(x, y, z, tool_angle, tool_radius, lengths, angles);
    Serial.print("Solved IK in "); Serial.print(millis() - prev_time); Serial.println(" ms");
    
    Serial.println("Angles:");
    for (int i = 0; i < 3; i++) {
        Serial.print(angles[i] * RAD_TO_DEG);
        Serial.print("\t");
    }
    Serial.println();
    
    prev_time = millis();
    forward_kinematics<mtx_type>(angles, lengths, end_position);
    Serial.print("Solved FK in "); Serial.print(millis() - prev_time); Serial.println(" ms");
    
    Serial.println("End position:");
    for (int i = 0; i < 3; i++) {
        Serial.print(end_position[i]);
        Serial.print("\t");
    }
    Serial.println();

    STACKPAINT_PRINT;
    FREERAM_PRINT;
    
    benchmark();
    
    STACKPAINT_PRINT;
    FREERAM_PRINT;
}

void loop() {

}
