#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
#define PI 3.14159265358979323846
const double TOL = 1e-6;  // Tolerance
const int MAX_ITER = 100; // iterations for Secant Method

// Define the first order differential equations of the system
double f1(double x, double y1, double y2) {
    return y2;  // dy1/dx = y2 = z
}

double f2(double x, double y1, double y2) {
    return -PI * PI / 4 * (y1 + 1);  // dy2/dx = -π^2/4 * (y1 + 1)
}

// Analytical solution
double analytical_solution(double x) {
    return cos(PI * x / 2) + 2 * sin(PI * x / 2) - 1;
}

// Runge-Kutta 4th order (RK4) for solving the IVP, with output logging to file
void RK4(double y1_0, double y2_0, double x_end, double dx, double& y1_end, ofstream& outfile) {
    double x = 0.0;
    double y1 = y1_0; // left Boundary
    double y2 = y2_0; // Guess of the derivative Z=y2

    // Write header for file with precision
    outfile << fixed << setprecision(4) << "x\t\t" << setprecision(5)
            << "numerical_u(x)\t" << setprecision(5) << "numerical_u'(x)\t"
            << setprecision(4) << "analytical_u(x)\n";

    while (x <= x_end) {
        outfile << setprecision(4) << x << "\t"
                << setprecision(5) << y1 << "\t"
                << setprecision(5) << y2 << "\t"
                << setprecision(4) << analytical_solution(x) << "\n";

        // Runge-Kutta calculations
        double k1_y1 = dx * f1(x, y1, y2);
        double k1_y2 = dx * f2(x, y1, y2);
        
        double k2_y1 = dx * f1(x + 0.5 * dx, y1 + 0.5 * k1_y1, y2 + 0.5 * k1_y2);
        double k2_y2 = dx * f2(x + 0.5 * dx, y1 + 0.5 * k1_y1, y2 + 0.5 * k1_y2);
        
        double k3_y1 = dx * f1(x + 0.5 * dx, y1 + 0.5 * k2_y1, y2 + 0.5 * k2_y2);
        double k3_y2 = dx * f2(x + 0.5 * dx, y1 + 0.5 * k2_y1, y2 + 0.5 * k2_y2);
        
        double k4_y1 = dx * f1(x + dx, y1 + k3_y1, y2 + k3_y2);
        double k4_y2 = dx * f2(x + dx, y1 + k3_y1, y2 + k3_y2);

        y1 += (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1) / 6.0;
        y2 += (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2) / 6.0;

        x += dx;
    }
    
    y1_end = y1;
}

// Secant Method to adjust initial slope y2_0
double shooting_method(double x_end, double dx, double target, ofstream& outfile) {
    double y1_0 = 0.0;     // Initial condition u(0) = 0
    double y2_0 = 0.5;     // First guess for u'(0)
    double y2_1 = 5.0;     // Second guess for u'(0)
    
    double y1_end_0, y1_end_1;
    RK4(y1_0, y2_0, x_end, dx, y1_end_0, outfile);
    RK4(y1_0, y2_1, x_end, dx, y1_end_1, outfile);

    cout << "Iter\t y2_0\t\t y2_1\t\t y1_end\t\t Error\n";
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        double error = y1_end_1 - target;
        cout << iter << "\t" << y2_0 << "\t" << y2_1 << "\t" << y1_end_1 << "\t" << error << endl;

        // Secant method iteration
        double y2_new = y2_1 - (y1_end_1 - target) * (y2_1 - y2_0) / (y1_end_1 - y1_end_0);
        
        // Update previous guesses
        y2_0 = y2_1;
        y1_end_0 = y1_end_1;
        y2_1 = y2_new;
        
        // Run RK4 with new guess
        RK4(y1_0, y2_1, x_end, dx, y1_end_1, outfile);
        
        // Check if we reached the target within tolerance
        if (fabs(y1_end_1 - target) < TOL) {
            return y2_1;  // Successful convergence
        }
    }
    
    cerr << "Shooting method did not converge within maximum iterations." << endl;
    return y2_1;
}

int main() {
    double x_end = 1.0;   // Boundary at x = 1
    double dx = 0.05;      // Step size for RK4
    double target = 1.0;  // Boundary condition u(1) = 1

    // save results
    ofstream outfile("results22.txt");
    if (!outfile) {
        cerr << "Error opening file for writing." << endl;
        return 1;
    }

    // Shooting method
    double y2_0 = shooting_method(x_end, dx, target, outfile);
    cout << "Initial slope u'(0) found: " << y2_0 << endl;
    outfile.close();

    return 0;
}
