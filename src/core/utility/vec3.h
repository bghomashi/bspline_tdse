#pragma once

#include <vector>
#include <cmath>

struct Vec3 {
    double x, y, z;
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline double length(const Vec3& a) {
    return std::sqrt(dot(a,a));
}
inline bool is_orthogonal(const Vec3& a, const Vec3& b) {
    return std::abs(dot(a,b)) < 1e-10;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    Vec3 out;
    out.x = a.y*b.z - a.z*b.y;
    out.y = a.z*b.x - a.x*b.z;
    out.z = a.x*b.y - a.y*b.x;
    return out;
}

inline Vec3 operator*(const Vec3& a, double b) {
    return Vec3{a.x*b, a.y*b, a.z*b};
}
inline Vec3 operator*(double b, const Vec3& a) {
    return a*b;
}
inline Vec3 operator/(const Vec3& a, double b) {
    return Vec3{a.x/b, a.y/b, a.z/b};
}
inline Vec3 operator+(const Vec3& a, const Vec3& b){
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 operator-(const Vec3& a, const Vec3& b){
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3& operator+=(Vec3& a, const Vec3& b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}
inline Vec3& operator-=(Vec3& a, const Vec3& b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}


inline void normalize(Vec3& a) {
    a = a/length(a);
}
inline Vec3 normal(const Vec3& a) {
    return a/length(a);
}
