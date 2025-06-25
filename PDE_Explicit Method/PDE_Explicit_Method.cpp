#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream> // Include for file operations

using namespace std;
double h, k, T, L;
// Exact solution function
double exact_solution(double x, double t) {
    if (x == 0.0 || x == L) {
        return 0.0;
    }
    return exp(-M_PI * M_PI * t) * sin(M_PI * x);
}

int main() {
    // Get user inputs
    
    cout << "Enter the spatial step size (h): ";
    cin >> h;
    cout << "Enter the time step size (k): ";
    cin >> k;
    cout << "Enter the total simulation time (T): ";
    cin >> T;
    cout << "Enter the length of the rod (L): ";
    cin >> L;

    double alpha = k / (h * h);

    // Warn about stability condition but do not terminate
    if (alpha >= 0.5) {
        cerr << "Warning: Stability condition violated (alpha >= 0.5). Results may be inaccurate.\n";
    }

    int nx = static_cast<int>(L / h) + 1;
    int nk = static_cast<int>(T / k) + 1;

    // Grid and time step indices
    int i, j;

    // Define the grid
    vector<vector<double> > u(nk, vector<double>(nx, 0.0));

    // Set initial and boundary conditions
    for (i = 0; i < nx; ++i) {
        double x = i * h;
        u[0][i] = sin(M_PI * x / L); // Initial condition
    }

    for (j = 1; j < nk; ++j) {
        u[j][0] = 0.0;          // Boundary at x=0
        u[j][nx - 1] = 0.0;     // Boundary at x=L
    }

    // Time integration loop (explicit scheme)
    for (j = 0; j < nk - 1; ++j) {
        for (i = 1; i < nx - 1; ++i) {
            u[j + 1][i] = alpha * u[j][i - 1] + (1 - 2 * alpha) * u[j][i] + alpha * u[j][i + 1];
        }
    }

    // Open a file to save the data
    ofstream outfile("Heat_Diffusion_distribution.dat");
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file for writing!" << endl;
        return 1;
    }

    // Write data to file in a format suitable for plotting
    outfile << "x t u_exact u_numerical\n"; // Header
    for (j = 0; j < nk; ++j) {
        double t = j * k; // Current time
        for (i = 0; i < nx; ++i) {
            double x = i * h;
            double exact = exact_solution(x, t);
            outfile << x << " " << t << " " << exact << " " << u[j][i] << "\n";
        }
        outfile << "\n"; // Separate time steps with a blank line
    }

    outfile.close();
    cout << "Data saved to temperature_distribution.dat" << endl;

    // Output results for a few selected time steps (optional)
    cout << "Numerical and Exact Solutions for selected time steps:\n";
    cout << "------------------------------------------------------\n";
    for (j = 0; j < nk; j += nk / 10) { // Output 10 evenly spaced time steps
        double t = j * k; // Current time
        cout << "Time = " << t << "\n";
        cout << "x          Numerical      Exact         Error\n";
        cout << "-----------------------------------------------\n";

        for (i = 0; i < nx; ++i) {
            double x = i * h;
            double exact = exact_solution(x, t);
            double error = abs(u[j][i] - exact);

            cout << setw(10) << x
                 << setw(15) << u[j][i]
                 << setw(15) << exact
                 << setw(15) << error << "\n";
        }
        cout << "-----------------------------------------------\n";
    }

    return 0;
}
