#include <stdio.h>
#include <math.h>
#ifndef EPS
#define EPS 2.0
#endif

double reg_coul(double x_1, double y_1, double z_1,
				 double x_2, double y_2, double z_2) {
	double radius, coulval;
	radius = sqrt((x_1 - x_2) * (x_1 - x_2) + (y_1 - y_2) * (y_1 - y_2) +
	 				(z_1 - z_2) * (z_1 - z_2));
	coulval =  1 / sqrt(radius * radius + EPS * EPS);
	return coulval;
}
