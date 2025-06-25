#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

// Constants
const double GM = 4 * M_PI * M_PI;  // gravitational constant (AU^3 / year^2)
const double r0 = 1.0;            // average Earth-Sun distance in AU
const double m = 1.0;             // mass of Jupiter in Earth masses
const double v0 = sqrt(GM / r0);    // ideal circular orbit velocity in AU/year

// Differential equations written in first-order form
double dx_dt(double vx) {
    return vx;  // dx/dt = vx
}

double dy_dt(double vy) {
    return vy;  // dy/dt = vy
}

double dvx_dt(double x, double y) {
    double r = sqrt(x * x + y * y);
    return -GM * x / pow(r, 3);  // dvx/dt = -GMx / r^3
}

double dvy_dt(double x, double y) {
    double r = sqrt(x * x + y * y);
    return -GM * y / pow(r, 3);  // dvy/dt = -GMy / r^3
}

// Runge-Kutta 4th-order step function
void rk4_step(double &x, double &y, double &vx, double &vy, double dt) {
    // k1 terms
    double k1x = dx_dt(vx) * dt;
    double k1y = dy_dt(vy) * dt;
    double k1vx = dvx_dt(x, y) * dt;
    double k1vy = dvy_dt(x, y) * dt;

    // k2 terms
    double k2x = dx_dt(vx + 0.5 * k1vx) * dt;
    double k2y = dy_dt(vy + 0.5 * k1vy) * dt;
    double k2vx = dvx_dt(x + 0.5 * k1x, y + 0.5 * k1y) * dt;
    double k2vy = dvy_dt(x + 0.5 * k1x, y + 0.5 * k1y) * dt;

    // k3 terms
    double k3x = dx_dt(vx + 0.5 * k2vx) * dt;
    double k3y = dy_dt(vy + 0.5 * k2vy) * dt;
    double k3vx = dvx_dt(x + 0.5 * k2x, y + 0.5 * k2y) * dt;
    double k3vy = dvy_dt(x + 0.5 * k2x, y + 0.5 * k2y) * dt;

    // k4 terms
    double k4x = dx_dt(vx + k3vx) * dt;
    double k4y = dy_dt(vy + k3vy) * dt;
    double k4vx = dvx_dt(x + k3x, y + k3y) * dt;
    double k4vy = dvy_dt(x + k3x, y + k3y) * dt;

    // Update position and velocity using weighted averages
    x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0;
    y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;
    vx += (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6.0;
    vy += (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6.0;
}

// Main function to simulate the orbit and save data to a file
void simulate_orbit(double velocity_ratio) {
    double v = velocity_ratio * v0;

    // Initial position and velocity
    double x = r0, y = 0;
    double vx = 0, vy = v ;
    
    // Time settings
    double dt = 1.0 / (365.25*100) ;  // time step
    double total_time = 3.0;         // 3 years
    int steps = total_time / dt;

    double r_min = r0, r_max = r0;
    double total_energy = 0, angular_momentum = 0;
    double a = 0, period = 0;

    // Open a file to store results
    string filename = "orbit_data_" + to_string(velocity_ratio) + ".txt";
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Header for the output file
    outfile << "# time(s) x(AU) y(AU) vx vy r(AU) kinetic_energy(J) potential_energy(J) total_energy(J) angular_momentum\n";

    // Run simulation and write results
    for (int i = 0; i < steps; ++i) {
        // Calculate the distance r
        double r = sqrt(x * x + y * y);

        // Track min/max radius for eccentricity calculation
        if (r < r_min) r_min = r;
        if (r > r_max) r_max = r;


        // Calculate energy and angular momentum
        double kinetic_energy = 0.5 * m * (vx * vx + vy * vy);
        double potential_energy = -(GM * m) / r;
        total_energy = kinetic_energy + potential_energy;
        angular_momentum = m * (x * vy - y * vx);

        // Save to file
        outfile << i * dt << " " << x << " " << y << " " << vx << " " << vy << " " << r << " "
                << kinetic_energy << " " << potential_energy << " "
                << total_energy << " " << angular_momentum << endl;

        // Perform one RK4 integration step
        rk4_step(x, y, vx, vy, dt);
    }

    // Eccentricity
    a = (r_min + r_max) / 2;  // Semi-major axis
    double eccentricity = (r_max - r_min) / (r_max + r_min);

    // Period
    period = 2 * M_PI * sqrt(a * a * a / GM);

    outfile.close();
    cout << fixed << setprecision(12) << endl;
    cout << "Final r_min: " << r_min << endl;
    cout << "Final r_max: " << r_max << endl;
    cout << "Eccentricity: " << eccentricity << endl;
    cout << "Period: " << period<< " years" << endl;
    cout << "Data for velocity ratio " << velocity_ratio << " saved to " << filename << endl;
    
    // Append period squared and semi-major axis cubed to a summary file
    string analysis_file = "orbit_analysis.txt";
    ofstream analysis_out(analysis_file, ios::app); // Open in append mode

    if (!analysis_out.is_open()) {
        cerr << "Error opening analysis file: " << analysis_file << endl;
        return;
    }

    // Write data for current velocity ratio
    analysis_out << velocity_ratio << " "
                 << pow(a, 3) << " "
                 << pow(period, 2) << endl;

    analysis_out.close();

}

int main() {
    vector<double> velocity_ratios;
    velocity_ratios.push_back(0.8);
    velocity_ratios.push_back(1.0);
    velocity_ratios.push_back(1.05);
    velocity_ratios.push_back(1.2);
    velocity_ratios.push_back(1.4);

    // Loop over velocity ratios
    for (size_t i = 0; i < velocity_ratios.size(); ++i) {
        double ratio = velocity_ratios[i];
        simulate_orbit(ratio);
    }
    return 0;
}
