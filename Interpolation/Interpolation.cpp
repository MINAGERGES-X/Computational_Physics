#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

// Function to calculate y = sin(x^2)
double Function(double x) {
    double k;
    k = sin(x * x);
    return k;
}

// Linear Interpolation Function
double Lag_lin(double x, const vector<double>& xi, const vector<double>& yi, int imax) {
    double y;
    int j = 0;

    // If x is outside the xi[] interval, return the nearest y value
    if (x <= xi[0]) return yi[0];
    if (x >= xi[imax - 1]) return yi[imax - 1];

    // Loop to find j so that xi[j-1] < x < xi[j]
    while (j <= imax - 1) {
        if (xi[j] >= x) break;
        j++;
    }
    y = yi[j - 1] + (yi[j] - yi[j - 1]) * (x - xi[j - 1]) / (xi[j] - xi[j - 1]);
    return y;
}

//  Lagrange 2nd order Polynomial Interpolation Function
double lag_Pol(double x, const vector<double>& xi, const vector<double>& yi) {
    double y = 0.0;
    int n = xi.size();

    for (int i = 0; i < n; ++i) {
        double L = 1.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                L *= (x - xi[j]) / (xi[i] - xi[j]);
            }
        }
        y += L * yi[i];
    }

    return y;
}
// Generate points for Interpolation and the plotting section
vector<double> generatePoints(double start, double end, double step) {
    vector<double> points;
    for (double i = start; i <= end; i += step) {
        points.push_back(i);
    }
    return points;
}

int main() {
    // Generate data points for interpolation
    vector<double> xi;
    vector<double> yi;
    double start = 1.0; // Starting point
    double end = 5.0;   // Ending point
    double step = 0.1;  // Step size " here you change the step size to control the number of generated data
    xi = generatePoints(start, end, step);

    // Calculate y = sin(x^2) for the given x values
    for (size_t i = 0; i < xi.size(); ++i) {
        yi.push_back(Function(xi[i]));
    }

    // Interpolation range (x0 points for plotting)
    vector<double> x0;
    for (double x = 1; x <= 5; x += 0.01) {
        x0.push_back(x);
    }

    // Storing interpolated values
    vector<double> y_linear, y_polyNom, y_original;
    
    // Interpolate at each x0 point
    for (size_t i = 0; i < x0.size(); ++i) {
        y_linear.push_back(Lag_lin(x0[i], xi, yi, xi.size()));
        y_polyNom.push_back(lag_Pol(x0[i], xi, yi));
        y_original.push_back(Function(x0[i]));
    }

    // Write data to a CSV file for Python plotting
    ofstream outputFile("interpolation_data.csv");

    // Write the header for the CSV file
    outputFile << "x,Original,Linear,2nd_PolyNomial\n";

    // Write the data into the CSV file
    for (size_t i = 0; i < x0.size(); ++i) {
        outputFile << x0[i] << ","
                   << y_original[i] << ","
                   << y_linear[i] << ","
                   << y_polyNom[i] << "\n";
    }

    // Close the CSV file
    outputFile.close();

    // Notify the user that the CSV file has been created
    cout << "Data written to interpolation_data.csv" << endl;

    return 0;
}
