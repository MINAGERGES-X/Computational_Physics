#include <iostream>
#include <cmath>
#include <limits> // Required for std::numeric_limits

using namespace std;

// Function definition
double f(double x) {
    double k;
    k = cos(2 * x) - 0.4 * x;
    return k;
}

double bisect(double a, double b, double eps)
{
    double xl, x0, xr;
    if (f(a) * f(b) > 0.0) {
        cout << "The function has the same signs at a and b. No root guaranteed in the interval." << endl;
        return 999; // Return NaN
    }

    xl = a;
    xr = b;
    int iteration = 0;
    cout << "Iteration\tCurrent Root Estimate" << endl;
    while (fabs(xr - xl) >= eps) {
        x0 = (xr + xl) / 2.0;
        cout << iteration << "\t\t" << x0 << endl;
        if ((f(xl) * f(x0)) <= 0.0)
            xr = x0;
        else
            xl = x0;
        iteration++;
    }
    return x0;
}

int main() {
  
    double a, b, eps;
    // Get user input
    cout << "Enter the value of a: ";
    cin >> a;
    
    cout << "Enter the value of b: ";
    cin >> b;
    
    cout << "Enter the tolerance epsilon: ";
    cin >> eps;
  
    double Bisect_root = bisect(a, b, eps);
    if (isnan(Bisect_root)) {
        cout << "No root found in the given interval." << endl;
    } else {
        cout << "The root is: " << Bisect_root << endl;
    }
    return 0;
}
