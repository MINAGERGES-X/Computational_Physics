#include <iostream>
#include <cmath>

using namespace std;

// Function definition
double fx(double x) {
    return cos(2 * x) - 0.4 * x;
}

// Derivative of f(x) for Newton's method
double fpx(double x) {
    return -2 * sin(2 * x) - 0.4;
}

// This function calculates both f(x) and f'(x) and returns them by reference
void fx_and_fpx(double x, double& fx_out, double& fpx_out) {
    fx_out = fx(x);   // Corrected: Use '=' instead of '=='
    fpx_out = fpx(x); // Corrected: Use '=' instead of '=='
}

double newton(void(*f)(double, double&, double&), double x, double eps, int& flag)
{
    double fx, fpx, xc;
    int i, iter = 1000;
    i = 0;
    
    cout << "Iteration\tCurrent Root Estimate" << endl;
    
    do {
        i = i + 1;
        f(x, fx, fpx); // Call the function pointer to get f(x) and f'(x)
        xc = x - fx / fpx; // Newton's update formula
        x = xc;
        
        // Print the current root estimate
        cout << i << "\t\t" << xc << endl;
        
        if (i >= iter) break; // Stop if we hit the iteration limit
    } while (fabs(fx) >= eps); // Check if the function value is close enough to zero
    
    flag = i; // Set the flag to indicate how many iterations were used
    if (i == iter) flag = 0; // If max iterations were hit, flag it as 0 (failure)
    return xc; // Return the approximate root
}

int main() {
    double a, eps;
    int flag; // To store the iteration count or failure
    
    // Get user input
    cout << "Enter the initial guess x: ";
    cin >> a;
    cout << "Enter the tolerance epsilon: ";
    cin >> eps;

    // Call the Newton method
    double Newton_root = newton(fx_and_fpx, a, eps, flag);
    
    if (flag != 0) {
        cout << "The root found by Newton's method is: " << Newton_root << endl;
        cout << "f(root) = " << fx(Newton_root) << endl; // Should be close to 0
    } else {
        cout << "Newton's method did not converge within the iteration limit." << endl;
    }

    return 0;
}
