# MetalPrinting
Aconity Metal 3D Printer Slicing Script

## Introduction
This MATLAB script is designed for slicing models to be used with the Aconity Metal 3D Printer. Developed by Dr. Connor Myant and Esperanza Shi at Imperial College London, Dyson School of Design Engineering in 2024, this script assists in preparing CLI files for printing.


## Contents

- **Main Script**: `SLICING.m`
This script is the main slicing tool for Aconity Metal 3D printers, also the main file where the parameters are set.\
@Line 46-51, Hatch_LaserPower/ScanSpeed will be overridden in the cli files where the LaserPowerMin/Max and ScanSpeedMin/Max are specified.\
@Line 203, you can change the number of divisions (how do you want the laser power&scan speed transit from the min to max).\

- **Function Script**: `BoundingBoxXYZ.m`
- **Function Script**: `Find_Raw_Slices_Vectorised.m`
- **Function Script**: `Sort_LayerLines.m`

- **Function Script**: `Add_Lines_Infill2Layer_turbo.m`
This is the file that generate a gradient laser power in the part. It\'92s also where you can be creative, just write your own script can call the new function in the Main Script.\

- **Function Script**: `Add_Support_Infill.m`
This script if for generating the infill. Parameters can be changed to adjust the density of the support.

- **Function Script**: `CLI_PLUS_ILT_Compiler_vk_2.m`
This script compiles coordinates of hatches/lines into Aconity format.

- **Function Script**: `CLI_PLUS_ILT_Compiler_vk_3.m`
This script only works with AconityStudio3.0 and above.

## Workflow

1. **Initial Setup**: Clear all existing variables, close all figures, and set up necessary parameters.

2. **Input Files**: Select the STL file you wish to slice using the MATLAB GUI. The STL file must be in millimeters (mm).

3. **Translations and Rotations**: This section ensures the part is placed correctly on the build plate and moves it to the bed center.

4. **Slice Generation**: Slice the model into layers and sort them into clockwise/counterclockwise order.

5. **Add Walls and Infill**: Add walls and infill patterns to each layer according to the specified print parameters.

6. **Support Structure Generation**: Optionally, generate support structures for overhangs and complex geometries.

7. **G-Code Compilation**: Compile the sliced layers into a single G-code file compatible with the Aconity printer.

8. **Visualization**: Optionally, visualize the sliced part and toolpaths for verification.

## Dependencies

- MATLAB software installed on your system.
- STL files in millimeters.
- Aconity Metal 3D Printer.

## Authors

- Dr. Connor Myant
- Esperanza (Qingyue) Shi
- Imperial College London, Dyson School of Design Engineering
