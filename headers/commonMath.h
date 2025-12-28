#ifndef COMMONMATH_H
#define COMMONMATH_H

// Linear interpolation function
double lerp(double y0, double x0, double y1, double x1, double x) {
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

#endif