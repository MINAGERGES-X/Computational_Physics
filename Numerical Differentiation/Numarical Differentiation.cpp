#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Defining the function itself
double f(double x) {
    return x * x * cos(x);
}

// First derivative
double analyticFirstDerivative(double x) {
    return 2 * x * cos(x) - x * x * sin(x);
}

// Second derivative
double analyticSecondDerivative(double x) {
    return 2 * cos(x) - 4 * x * sin(x) - x * x * cos(x);
}

// Numerical methods for derivatives:

// First Derivative:
// Forward difference for first derivative
double FD1(double x, double h) {
    return (f(x + h) - f(x)) / h;
}

// Backward difference for first derivative
double BD1(double x, double h) {
    return (f(x) - f(x - h)) / h;
}

// Central difference for first derivative
double CD1(double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
}


// Second Derivative:
// Forward difference for second derivative
double FD2(double x, double h) {
    return (f(x + 2 * h) - 2 * f(x + h) + f(x)) / (h * h);
}

// Backward difference for second derivative
double BD2(double x, double h) {
    return (f(x) - 2 * f(x - h) + f(x - 2 * h)) / (h * h);
}

// Central difference for second derivative
double CD2(double x, double h) {
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
}

// Richardson extrapolation for first derivativeusing the first two values of the list
double richFirst(double x, double h1, double h2) {
    double R1 = CD1(x, h1);
    double R2 = CD1(x, h2);
    
    // Apply the modified Richardson extrapolation formula
    return R2 + (R2 - R1) / (pow((h1 / h2), 2) - 1);
}

// Richardson extrapolation for second derivative
double richSecond(double x, double h1, double h2) {
    double R3 = CD2(x, h1);
    double R4 = CD2(x, h2);
    
    // Apply the modified Richardson extrapolation formula
    return R4 + (R4 - R3) / (pow((h1 / h2), 2) - 1);
}


int main() {
    // Points to compute derivatives at
    double x1 = 1.0;
    double x2 = 2.0;
    
    // Step sizes
    double h[] = {0.5, 0.25, 0.125};

    // Output headers ( x_value, Numerical_method, h_vlaue, First_Derivative, Second_ Derivative
    cout << fixed << setprecision(10);
    cout << "x\tNumerical_Method\t\th\tFirst_Derivative\tSecond_Derivative\n";
    cout << "----------------------------------------------------------------------\n";
    
    // This loop is to assign the calculation for x1 and x2 values with the three values of h
    double points[] = {x1, x2};
    for (int i = 0; i < 2; ++i) {
        double point = points[i];
        cout << "At x = " << point << "\n";
        for (int i=0; i<3; ++i){
            double step = h[i];
            // Forward, Backward, Central, Richardson methods
            cout << point << "\tForward\t\t" << step << "\t"
                 << FD1(point, step) << "\t\t"
                 << FD2(point, step) << "\n";

            cout << point << "\tBackward\t" << step << "\t"
                 << BD1(point, step) << "\t\t"
                 << BD2(point, step) << "\n";

            cout << point << "\tCentral\t\t" << step << "\t"
                 << CD1(point, step) << "\t\t"
                 << CD2(point, step) << "\n";

            // Richardson extrapolation
            // "if" condition to select only th efirst two value of h as {h1 =0.5 & h2= 0.25)
            if (i < 2) {
                double next_step = h[i+1];
                cout << point << "\tRichardson\t" << step << "\t"
                     << richFirst(point, step, next_step) << "\t\t"
                     << richSecond(point, step, next_step) << "\n";
            }
        }

        // Priniting the analytic solution for the first and the second derivatives
        cout << "\nAnalytical first derivative at x = " << point << ": "
             << analyticFirstDerivative(point) << "\n";
        cout << "Analytical second derivative at x = " << point << ": "
             << analyticSecondDerivative(point) << "\n";
        cout << "----------------------------------------------------------------------\n";
    }

    return 0;
}
