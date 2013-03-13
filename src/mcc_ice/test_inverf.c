#include <stdio.h>

#include "inverf.h"

int main() {
	double x, y;
	long n=100;

	for (long i=1; i<n; i++) {
		x=1-2.*i/n;
		y=inverf(x);
		printf("%f %f\n", x, y);
	}
}

