# Hartmann Flow Simulation using Two-Fluid Theory

This project simulates a plasma using the two-fluid theory, specifically focusing on the Hartmann flow. The code allows you to simulate the velocity profile of the plasma and output the results for further analysis.

## How to Use

1. **Modify Constants**: To run the simulation with different parameters, modify the constants in the `constants.hpp` file. This file contains all the necessary physical and numerical constants for the simulation.

2. **Compile the Code**: After updating the constants (if needed), compile the code using `make`:

   ```bash
   make

3. **Run the Simulation**: Once compiled, run the program to generate the data output:

   ```bash
   ./Hartmann > datos.dat

4. **Analyze the Data**: The output is saved in `datos.dat`. To plot the velocity profile of the plasma, you can use the first two columns of the output data file.

## Requirements

- **C++ compiler**
- **Make utility**

## Contact

For any issues or questions regarding this project, feel free to contact the developer at:

- **Email**: jubeltranr@unal.edu.co

