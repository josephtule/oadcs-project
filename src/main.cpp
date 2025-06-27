#include <array>
#include <iostream>

#include "EulerAngleAtt.h"
#include "QuaternionAtt.h"
#include "units.h"

int main() {

    vec3 eas = vec3(.5, 0, .25);
    std::array<int, 3> seq = {3, 2, 1};
    EulerAngleAtt ea_att(eas, seq, UnitsAngle::RADIANS);
    std::cout << ea_att.to_string() << std::endl;

    vec4 q_init = vec4(0.1, 0.2, 0.3, 0.9);
    QuaternionAtt q_att(q_init);
    std::cout << q_att.to_string() << std::endl;
    q_att.from_dcm(mat3::Identity());
    std::cout << q_att.to_string() << std::endl;

    return 0;
}
