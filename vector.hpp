#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>

class vector3D {
private:
    double X, Y, Z;
public:
    // Initialize the vector
    void load(double x0, double y0, double z0);
    // Get the components
    double x(void) const { return X; }
    double y(void) const { return Y; }
    double z(void) const { return Z; }
    // Show the vector
    void show(void) const;
    //-------------------------
    // Vectorial operators
    //-------------------------
    // Equal
    void operator=(vector3D v2);
    // Sum
    vector3D operator+(vector3D v2) const;
    void operator+=(vector3D v2);
    // Subtraction
    vector3D operator-(vector3D v2) const;
    void operator-=(vector3D v2);
    // Scalar multiplication
    vector3D operator*(double a) const;
    void operator*=(double a);
    friend vector3D operator*(double a, vector3D v1);
    // Scalar division
    vector3D operator/(double a) const;
    void operator/=(double a);
    // Dot product
    double operator*(vector3D v2) const;
    // Cross product
    vector3D operator^(vector3D v2) const;
    // Norm operations
    double norm2(void) const;
    double norm(void) const;
    // Angle between two vectors
    friend double angle(vector3D v1, vector3D v2);
};

#endif // VECTOR_H

