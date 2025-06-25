
# MINAGERGES-X 🧮
### Computational Physics

A versatile suite of numerical methods implemented in C++. Each folder includes a `.cpp` implementation, sample results, and a detailed report.

## 📁 Directory Structure

* **BVP\_Shooting Method**
  Solves boundary value problems (BVPs) using the shooting method.

* **Harmonic Oscillator**
  Simulates simple & damped harmonic oscillators using numerical integration.

* **Interpolation**
  Implements linear, polynomial, and spline interpolation techniques.

* **Numerical Differentiation**
  Computes derivatives using forward, backward, and central difference formulas.

* **Numerical Integration**
  Supports trapezoidal, Simpson’s rules, and Gaussian quadrature methods.

* **PDE\_Explicit Method**
  Solves partial differential equations (PDEs) using explicit finite difference schemes.

* **Planetary Motion**
  Models orbital motion of planets using gravitational equations and numerical solvers.

* **Roots of Non‑Linear Equations**
  Finds roots using methods such as bisection, Newton–Raphson, and secant.

## 🛠️ Prerequisites

* C++ compiler (e.g., `g++` or `clang++`) with **C++11** or higher
* CMake (optional—see below)
* Standard C++ library support

## 🔧 Building & Running

### Using Make (for each folder)

1. `cd <folder_name>`
2. `g++ -std=c++11 *.cpp -o run`
3. `./run`

### Using CMake (optional)

```bash
mkdir build
cd build
cmake ..
make
./<executable_name>
```

Adjust `<executable_name>` to match your project folder.

## ✅ Contents of Each Folder

* `*.cpp` — Implementations of specific numerical algorithms.
* `results/` — Input/output sample files illustrating expected outcomes.
* `report.pdf` (or `.md`) — Documents methodology, analytical derivations, and benchmarking/results.

## 📚 Usage & Examples

Each program:

1. Reads predefined parameters (can be modified inside the `.cpp` or via config files).
2. Runs the numerical method.
3. Outputs results to terminal and/or files in `results/`.
4. Example results are included for comparison (see each folder’s report).

## 🎯 Purpose & Scope

This repository acts as a teaching and reference toolkit for:

* Learning numerical analysis methods
* Comparing algorithm efficiency and accuracy
* Visual demonstrations of method performance


## 🔗 References

Each folder’s report includes:

* Mathematical derivations
* Error analysis and convergence studies
* Bibliographical references for deeper reading

## 📝 License

[MIT License](LICENSE)
