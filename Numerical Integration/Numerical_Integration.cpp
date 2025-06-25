#include <iostream>
#include <cmath>

using namespace std;


// Function to integrate
double f(double x) {
    double k;
    k = sin(x); 
    return k ;
}

// Trapezoidal rule 
double trapezoidal(double a, double b, int n) {
   
    double h = (b - a) / n;
    
   
    double result = 0.5 * h * (f(a) + f(b));
    
    // loop over n = 1 to n - 1 and sum up the terms
    for (int i = 1; i < n; ++i) {
        result += h * f(a + i * h);
    }
    
    return result;
}


// Simpson's rule
double simpson(double a, double b, int n) {
    if (n % 2 != 0) {
        n++; // Ensure that n is even for Simpson's rule
    }
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i += 2) { // terms multiplied by 4
        sum += 4 * f(a + i * h);
    }
    for (int i = 2; i < n; i += 2) {    // terms multiplied by 2
        sum += 2 * f(a + i * h);
    }
    return sum * h / 3;
}


// Gaussian Quadrature (4-point method)
double gaussianQuadrature(double a, double b) {
    // Coefficients and abscissas for 4-point Gaussian Quadrature
    const int n = 4;
    double x[n] = {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
    double ci[n] = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};

    // Change of interval to map [-1,1] to [a,b]
    double r = 0.0;
    double c = 0.5 * (b + a);
    double m = 0.5 * (b - a);

    // loop to apply the Gaussian Quadrature formula
    for (int i = 0; i < n; ++i) {
        r += ci[i] * f(m * x[i] + c);
    }
    return r * m;
}
int main() {
double a, b;  // Declare variables for limits
    double n;
    
    //double a = 0.0;   // Lower limit
    //double b = M_PI;  // Upper limit
    //int n = 4;      // Number of steps
// Input for lower limit
    cout << "Lower limit: ";
    cin >> a;
    
    // Input for upper limit
    cout << "Upper limit: ";
    cin >> b;
    
    // Input for number of subintervals
    cout << "Subintervals: ";
    cin >> n;
    
    // Calculate the integral using the Trapezoidal rule
    double trap_result = trapezoidal(a, b, n);
    double simpson_result = simpson(a, b, n);
    double gauss_result = gaussianQuadrature(a, b);
    // Comparison to gaussian quadrature result
    double diff_trap = gauss_result - trap_result;
    double diff_simpson = gauss_result - simpson_result;
    // Display the result

    cout << "Trapezoidal rule result: " << trap_result <<endl;
    cout << "simpson rule result: " << simpson_result <<endl;
    cout << "Gaussian Quadrature result: " << gauss_result << endl;
    cout << "Exact result: 2.0" << endl;
    cout << " Difference for Trap.:"<< diff_trap << endl;
    cout << " Difference for simspon's.:"<< diff_simpson << endl;

    return 0;
}
