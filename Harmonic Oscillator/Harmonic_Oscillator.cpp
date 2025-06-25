#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double euler2d(double(*d1x)(double, double, double),
               double(*d2x)(double, double, double),
               double ti, double xi, double vi, double tf,
               double& xf, double& vf, bool modified) {
    double dt = tf - ti;
    if (modified) {
        // Modified Euler method
        xf = xi + d1x(ti, xi, vi) * dt; //predicted
        vf = vi + d2x(ti, xi, vi) * dt; //predicted
        xf = xi + (d1x(ti, xi, vi) + d1x(ti, xf, vf)) * 0.5 * dt; // corrected
        vf = vi + (d2x(ti, xi, vi) + d2x(ti, xf, vf)) * 0.5 * dt; //corrected
    } else {
        // Simple Euler method
        xf = xi + d1x(ti, xi, vi) * dt;
        vf = vi + d2x(ti, xi, vi) * dt;
    }
    return 0.0;
}

double f1(double t, double x, double v) {
    return v; // dx/dt = v
}

double f2(double t, double x, double v) {
    return -x; // dv/dt = -x (harmonic oscillator)
}

int main() {
    // Initial conditions
    double A = 1.0; // amplitude
    double dt = 0.01; // time step
    double tf = 100.0; // final time
    int n_steps = static_cast<int>(tf / dt);
    
    // Vectors to store results
    vector<double> time(n_steps), simple_x(n_steps), simple_v(n_steps);
    vector<double> modified_x(n_steps), modified_v(n_steps), analytical_x(n_steps);
    vector<double> simple_energy(n_steps), modified_energy(n_steps);
    
    // Initial conditions
    simple_x[0] = A;
    simple_v[0] = 0.0;
    modified_x[0] = A;
    modified_v[0] = 0.0;
    
    // Initial energy calculation
    simple_energy[0] = 0.5 * simple_v[0] * simple_v[0] + 0.5 * simple_x[0] * simple_x[0];
    modified_energy[0] = 0.5 * modified_v[0] * modified_v[0] + 0.5 * modified_x[0] * modified_x[0];

    // Time integration using Simple Euler
    for (int i = 1; i < n_steps; i++) {
        time[i] = i * dt;
        euler2d(f1, f2, time[i - 1], simple_x[i - 1], simple_v[i - 1], time[i], simple_x[i], simple_v[i], false);
        
        // Calculate energy for Simple Euler
        simple_energy[i] = 0.5 * simple_v[i] * simple_v[i] + 0.5 * simple_x[i] * simple_x[i];
    }

    // Time integration using Modified Euler
    for (int i = 1; i < n_steps; i++) {
        euler2d(f1, f2, time[i - 1], modified_x[i - 1], modified_v[i - 1], time[i], modified_x[i], modified_v[i], true);
        
        // Calculate energy for Modified Euler
        modified_energy[i] = 0.5 * modified_v[i] * modified_v[i] + 0.5 * modified_x[i] * modified_x[i];
    }

    // Analytical solution
    for (int i = 0; i < n_steps; i++) {
        analytical_x[i] = A * cos(time[i]);
    }

    // Write results to files for plotting
    ofstream file("results.txt");
    file << "Time\tSimple_Euler_X\tModified_Euler_X\tAnalytical_X\tSimple_Energy\tModified_Energy\n";
    for (int i = 0; i < n_steps; i++) {
        file << time[i] << "\t" << simple_x[i] << "\t" << modified_x[i] << "\t" << analytical_x[i] << "\t"
             << simple_energy[i] << "\t" << modified_energy[i] << "\n";
    }
    file.close();

    cout << "Results written to results.txt" << endl;
    return 0;
}
