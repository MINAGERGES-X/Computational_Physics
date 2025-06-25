#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// Define the function
double f(double x) {
    double k;
    k = x - cos(x);
    return k;
}

// Bisection Method
void bisection(double a, double b, double tol, int& iter, ofstream& outfile) {
    if (f(a) * f(b) >= 0) {
        outfile << "Bisection method: Incorrect initial guesses." << endl;
        return;
    }

    double c;
    iter = 0;
    outfile << "Iteration\tCurrent Root Estimate" << endl;

    while ((b - a) >= tol) {
        c = (a + b) / 2.0;
        iter++;

        if (f(c) == 0.0) {
            break;
        } else if (f(c) * f(a) < 0) {
            b = c;
        } else {
            a = c;
        }

        outfile << iter << "\t\t\t\t" << c << endl;
    }

    outfile << "Bisection method root: " << c << endl;
}

// False Position Method
void False_P(double a, double b, double tol, int max_iter, ofstream& outfile) {
    double x2;
    int iter = 0;
    outfile << "Iteration\tCurrent Root Estimate" << endl;

    while (iter < max_iter) {
        if (f(a) == f(b)) {
            outfile << "False_Position method: Division by zero error." << endl;
            return;
        }

        // Apply False position formula
        x2 = b - (f(b) * (b - a)) / (f(b) - f(a));

        if (fabs(x2 - b) < tol) {
            outfile << "False_Position method root: " << x2 << endl;
            return;
        }

        // Update a and b
        a = b;
        b = x2;
        iter++;
        outfile << iter << "\t\t\t\t" << x2 << endl;
    }

    outfile << "False_Position method did not converge within the maximum iterations." << endl;
}

// Main function
int main() {
    double a, b;
    double tol = 1e-6;
    int max_iter = 1000;
    int iter = 0;

    ofstream outfile("HW3_Q2_Results.txt");
    if (!outfile) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    outfile << fixed << setprecision(6);

    // Input interval [a, b]
    cout << "Enter the interval [a, b] for both methods: ";
    cin >> a >> b;

    // Apply Bisection Method
    outfile << "Calculating the root with bisection_Method.........." << endl;

    bisection(a, b, tol, iter, outfile);
  
    
   

    // Apply Fals Position Method
    outfile << "Calculating the root with false_position_Method.........." << endl;
    False_P(a, b, tol, max_iter, outfile);

    outfile.close();
    return 0;
}
