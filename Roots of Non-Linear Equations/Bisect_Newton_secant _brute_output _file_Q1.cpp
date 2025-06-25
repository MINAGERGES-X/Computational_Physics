#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;


// Define the function
double f(double x) {
    double k;
    k = cos(2 * x) - 0.4 * x;
    return k;
}
// Define the derivative
double df(double x) {
    double D;
    D = -2 * sin(2 * x) - 0.4;
    return D;
}
// Bisection Method
void bisection(double a, double b, double tol, int& iter, ofstream& outfile) {
    if (f(a) * f(b) >= 0) {
        outfile << "Bisection method: Incorrect initial guesses." << endl;
        return;
    }

    double c = a;
    iter = 0;
    outfile << "Iteration\tCurrent Root Estimate" << endl;

    while ((b - a) >= tol) {
       
        c = (a + b) / 2;
        iter++;

        if (f(c) == 0.0) {
            break;
        }
        else if (f(c) * f(a) < 0) {
            b = c;
        } else {
            a = c;
        }
        
        outfile << iter << "\t\t" << c << endl; // Output current estimate
    }

    outfile << "Bisection method root: " << c << endl;
}

// Newton's Method
void newton(double x0, double tol, int max_iter, ofstream& outfile) {
    double h = f(x0) / df(x0);
    int iter = 0;
    outfile << "Iteration\tCurrent Root Estimate" << endl;
    
    while (abs(h) >= tol && iter < max_iter) {
        h = f(x0) / df(x0);
        x0 = x0 - h;
        outfile << iter + 1 << "\t\t" << x0 << endl;
        iter++;
    }
    
    if (iter == max_iter)
        outfile << "Newton's method did not converge within the maximum iterations." << endl;
    else
        outfile << "Newton's method root: " << x0 << endl;
}

// Secant Method
void secant(double x0, double x1, double tol, int max_iter, ofstream& outfile) {
    double x2;
    int iter = 0;
    outfile << "Iteration\tCurrent Root Estimate" << endl;

    while (iter < max_iter) {
        if (f(x0) == f(x1)) {
            outfile << "Secant method: Division by zero error." << endl;
            return;
        }

        // Apply Secant formula
        x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));

        if (abs(x2 - x1) < tol) {
            outfile << "Secant method root: " << x2 << endl;
            return;
        }

        // Update x0 and x1
        x0 = x1;
        x1 = x2;
        iter++;
        outfile << iter << "\t\t" << x2 << endl;
    }
    outfile << "Secant method did not converge within the maximum iterations." << endl;
}

// Brute Force Method
void bruteForce(double start, double end, double step, ofstream& outfile) {
    double x1 = start;
    double f1 = f(x1);

    for (double x2 = start + step; x2 <= end; x2 += step) {
        double f2 = f(x2);

        // Check for sign change
        if (f1 * f2 < 0) {

            double root = (x1 + x2) / 2.0;
            cout << "Root found: x = " << setprecision(6) << root << endl;
            outfile << "Root found: x = "<< "\t\t" << root<< endl;
        }

        // Update x1 and f1
        x1 = x2;
        f1 = f2;
    }
}

// Main function
int main() {
    double a, b;
    double tol = 1e-6;
    int max_iter = 1000;
    double start = -2, end = 1, step = 0.01;
    int iter = 0;

    ofstream outfile("HW3_Q1_Results.txt");
    if (!outfile) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    outfile << fixed << setprecision(6);
    outfile << "Enter the interval [a, b] for the Bisection Method: ";
    cout << "Enter the interval [a, b] for the Bisection Method: ";
    cin >> a >> b;

    // Apply Bisection Method
    bisection(a, b, tol, iter, outfile);

    // Apply Newton's Method
    double x0[3];
    cout << "Enter three initial guesses for Newton's Method: ";
    outfile << "Enter three initial guesses for Newton's Method: ";
    for (int i = 0; i < 3; ++i) {
        cin >> x0[i];
        outfile << "Initial guess " << i + 1 << ": " << x0[i] << endl;
    }
    for (int i = 0; i < 3; ++i) {
        newton(x0[i], tol, max_iter, outfile);
    }

    // Apply Secant Method
    double x1;
    cout << "Enter two initial guesses for Secant Method (x0 and x1): ";
    outfile << "Enter two initial guesses for Secant Method (x0 and x1): ";
    cin >> x0[0] >> x1;
    secant(x0[0], x1, tol, max_iter, outfile);

    // Apply Brute Force Method
    outfile << "Brute force method in interval [" << start << ", " << end << "] with step size " << step << ":" << endl;
    bruteForce(start, end, step, outfile);

    outfile.close();
    return 0;
}
