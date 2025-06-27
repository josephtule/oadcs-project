#include "attitude.h"

mat3 skew(vec3 v) {
    mat3 s;
    s << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;

    return s;
}