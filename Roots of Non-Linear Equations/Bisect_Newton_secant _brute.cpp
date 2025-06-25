#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// Define the function f(x) = cos(2x) - 0.4x
double f(double x) {
    double k;
    k = cos(2 * x) - 0.4 * x;
    return k;
}
// Define the derivative of the function f'(x) = -2sin(2x) - 0.4
double df(double x) {
    double D;
    D = -2 * sin(2 * x) - 0.4;
    return D;
}

// Bisection Method
void bisection(double a, double b, double tol, int& iter) {
    if (f(a) * f(b) >= 0) {
        cout << "Bisection method: Incorrect initial guesses." << endl;
        return;
    }

    double c = a;
    iter = 0; // Initialize iter here
    cout << "Iteration\tCurrent Root Estimate" << endl;

    while ((b - a) >= tol) {
        // Find the midpoint
        c = (a + b) / 2;
        iter++; // Increment the iteration count

        // Check if the midpoint is a root
        if (f(c) == 0.0) {
            break;
        }
        // Decide the side to repeat the steps
        else if (f(c) * f(a) < 0) {
            b = c;
        } else {
            a = c;
        }
        
        cout << iter << "\t\t" << c << endl; // Output current estimate
    }

    cout << "Bisection method root: " << c << endl;
}

// Newton's Method
void newton(double x0, double tol, int max_iter) {
    double h = f(x0) / df(x0);
    int iter = 0;
    cout << "Iteration\tCurrent Root Estimate" << endl;
    
    while (abs(h) >= tol && iter < max_iter) {
        h = f(x0) / df(x0);
        x0 = x0 - h;
        cout << iter + 1 << "\t\t" << x0 << endl; // Corrected to show current iteration
        iter++;
    }
    
    if (iter == max_iter)
        cout << "Newton's method did not converge within the maximum iterations." << endl;
    else
        cout << "Newton's method root: " << x0 << endl;
}

// Secant Method
void secant(double x0, double x1, double tol, int max_iter) {
    double x2;
    int iter = 0;
    cout << "Iteration\tCurrent Root Estimate" << endl;

    while (iter < max_iter) {
        if (f(x0) == f(x1)) {
            cout << "Secant method: Division by zero error." << endl;
            return;
        }

        // Apply Secant formula
        x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));

        if (abs(x2 - x1) < tol) {
            cout << "Secant method root: " << x2 << endl;
            return;
        }

        // Update x0 and x1
        x0 = x1;
        x1 = x2;
        iter++;
        cout << iter << "\t\t" << x2 << endl; // Output current estimate
    }
    cout << "Secant method did not converge within the maximum iterations." << endl;
}

// Brute Force Method
void bruteForce(double start, double end, double step) {
    double x1 = start;
    double f1 = f(x1);

    for (double x2 = start + step; x2 <= end; x2 += step) {
        double f2 = f(x2);

        // Check for sign change indicating a root between x1 and x2
        if (f1 * f2 < 0) {
            // Use the midpoint to approximate the root
            double root = (x1 + x2) / 2.0;
            cout << "Root found: x = " << setprecision(6) << root << endl;
        }

        // Update x1 and f1 for next iteration
        x1 = x2;
        f1 = f2;
    }
}


// Main function
int main() {
    double a, b;
    double tol = 1e-6;
    int max_iter = 1000;
    double start = -10, end = 10, step = 0.01; // Interval for brute force method
    int iter = 0; // Variable for iteration count

    cout << fixed << setprecision(6);
    cout << "Enter the interval [a, b] for the Bisection Method: ";
    cin >> a >> b;

    // Apply Bisection Method
    bisection(a, b, tol, iter);

    // Apply Newton's Method
    double x0[3]; // Array for three initial guesses
    cout << "Enter three initial guesses for Newton's Method: ";
    for (int i = 0; i < 3; ++i) {
        cin >> x0[i];
    }
    for (int i = 0; i < 3; ++i) {
        newton(x0[i], tol, max_iter); // Call newton for each guess
    }

    // Apply Secant Method
    double x1;
    cout << "Enter two initial guesses for Secant Method (x0 and x1): ";
    cin >> x0[0] >> x1; // Corrected to use x0[0]
    secant(x0[0], x1, tol, max_iter);

    // Apply Brute Force Method
    cout << "Brute force method in interval [" << start << ", " << end << "] with step size " << step << ":" << endl;
    bruteForce(start, end, step);

    return 0;
}
