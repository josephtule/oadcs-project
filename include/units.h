#pragma once
#include "constants.h"
// #include "typedefs.h"
// #include <iostream>

enum class UnitsAngle {
    RADIANS,
    DEGREES,
    ARCSEC,
    MINUTES,
};

enum class UnitsLinear {
    MILLIMETER,
    METER,
    CENTIMETER,
    KILOMETER,
    FOOT,
    INCH,
    AU,
};

enum class UnitsTime {
    YEAR,
    MONTH,
    DAY,
    HOUR,
    MINUTE,
    SECOND,
};

enum class UnitsMass {
    KILOGRAMS,
    GRAMS,
    MILLIGRAMS,
    MICROGRAMS,
    METRIC_TON,
};

template <class T>
auto convertAngle(T val, UnitsAngle units_in, UnitsAngle units_out) -> T {
    using enum UnitsAngle;
    // convert to radians
    switch (units_in) {
    case RADIANS:
        break;
    case DEGREES:
        val *= deg2rad;
        break;
    case MINUTES:
        val *= deg2rad / 60.;
        break;
    case ARCSEC:
        val *= deg2rad / 3600.;
        break;
    }

    // convert to output
    switch (units_out) {
    case RADIANS:
        break;
    case DEGREES:
        val *= rad2deg;
        break;
    case MINUTES:
        val *= rad2deg * 60.;
        break;
    case ARCSEC:
        val *= rad2deg * 3600.;
        break;
    }

    return val;
}