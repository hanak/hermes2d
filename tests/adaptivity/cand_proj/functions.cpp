#include <cmath>
#include <hermes2d.h>
#include "functions.h"

//x2
scalar func_x2_val(double x, double y) { return x*x; };
scalar func_x2_dx(double x, double y) { return 2*x; };
scalar func_x2_dy(double x, double y) { return 0; };

//x2y2
scalar func_x2y2_val(double x, double y) { return x*x * y*y; };
scalar func_x2y2_dx(double x, double y) { return 2*x * y*y; };
scalar func_x2y2_dy(double x, double y) { return x*x * 2*y; };

//x3y1
scalar func_x3y1_val(double x, double y) { return x*x*x * y + x*x; };
scalar func_x3y1_dx(double x, double y) { return 3*x*x * y + 2*x; };
scalar func_x3y1_dy(double x, double y) { return x*x*x; };

//x3y4
scalar func_x3y4y_val(double x, double y) { return x*x*x * y*y*y*y + y; };
scalar func_x3y4y_dx(double x, double y) { return 3*x*x * y*y*y*y; };
scalar func_x3y4y_dy(double x, double y) { return x*x*x * 4*y*y*y + 1; };

