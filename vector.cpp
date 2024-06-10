#include "vector.hpp"

// Initialize the vector
void vector3D::load(double x0, double y0, double z0) {
    X = x0;
    Y = y0;
    Z = z0;
}

// Show the vector
void vector3D::show(void) const {
    std::cout << "(" << X << ", " << Y << ", " << Z << ")\n";
}

// Equal
void vector3D::operator=(vector3D v2) {
    X = v2.X;
    Y = v2.Y;
    Z = v2.Z;
}

// Sum
vector3D vector3D::operator+(vector3D v2) const {
    vector3D total;
    total.X = X + v2.X;
    total.Y = Y + v2.Y;
    total.Z = Z + v2.Z;
    return total;
}

void vector3D::operator+=(vector3D v2) {
    X += v2.X;
    Y += v2.Y;
    Z += v2.Z;
}

// Subtraction
vector3D vector3D::operator-(vector3D v2) const {
    vector3D total;
    total.X = X - v2.X;
    total.Y = Y - v2.Y;
    total.Z = Z - v2.Z;
    return total;
}

void vector3D::operator-=(vector3D v2) {
    X -= v2.X;
    Y -= v2.Y;
    Z -= v2.Z;
}

// Scalar multiplication
vector3D vector3D::operator*(double a) const {
    vector3D total;
    total.X = X * a;
    total.Y = Y * a;
    total.Z = Z * a;
    return total;
}

vector3D operator*(double a, vector3D v1) {
    v1.X *= a;
    v1.Y *= a;
    v1.Z *= a;
    return v1;
}

void vector3D::operator*=(double a) {
    X *= a;
    Y *= a;
    Z *= a;
}

// Scalar division
vector3D vector3D::operator/(double a) const {
    vector3D total;
    total.X = X / a;
    total.Y = Y / a;
    total.Z = Z / a;
    return total;
}

void vector3D::operator/=(double a) {
    X /= a;
    Y /= a;
    Z /= a;
}

// Dot product
double vector3D::operator*(vector3D v2) const {
    return X * v2.X + Y * v2.Y + Z * v2.Z;
}

// Cross product
vector3D vector3D::operator^(vector3D v2) const {
    vector3D total;
    total.X = Y * v2.Z - Z * v2.Y;
    total.Y = Z * v2.X - X * v2.Z;
    total.Z = X * v2.Y - Y * v2.X;
    return total;
}

// Norm operations
double vector3D::norm2(void) const {
    return X * X + Y * Y + Z * Z;
}

double vector3D::norm(void) const {
    return std::sqrt(norm2());
}

// Angle between two vectors
double angle(vector3D v1, vector3D v2) {
    return std::acos((v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z) / (std::sqrt(v1.norm2()) * std::sqrt(v2.norm2())));
}
