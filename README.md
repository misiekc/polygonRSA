# polygonRSA
A program to generate saturated random packings of arbitrary polygons and rounded polygons according to random sequential adsorption protocol. It is a stipped down version
of our full RSA simulation package https://github.com/misiekc/rsa3d

## Compilation
The project uses CMake build system. To prepare build files for the Release build in the `build` folder, go to project root folder and run:
```bash
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
```

GCC with C++17 support is required (version 7). If `gcc` and `g++` commands are aliased to a lower version, or to clang (the default behaviour on macOS), you have
to specify the path manually using `cmake` options `-DCMAKE_C_COMPILER=...`, `-DCMAKE_CXX_COMPILER=...`.

Now, in the `build` folder, you can use `make polygonRSA` to build the program. You can use `-j` option with a number of threads for parallel compilation.

## Running simulations

### Command line

The program has multiple modes, which has to be specified as the first positional argument. The brackets represent additional, mode-specific positional arguments.
* `simulate` - creates new packing(s)
* `test [*.bin packing file]` - tries to insert new particles to an existing packing. It can be used in debugging to check, if the packing is actually saturated
* `dat [directory with packings]` - generates `*.dat` file from packing `*.bin` in a given directory
* `povray [*.bin packing file]` - generates Povray input file to draw a packing 
* `wolfram [*.bin packing file]` - generates Wolfram Mathematice input file to draw a packing. For larger packings it is better to copy the contents and paste it
in Mathematica notebook manually for performance reasons rather that just open the file
* `bc_expand [*.bin packing input file] [*.bin packing output file]` - adds periodically reflected shapes near all 4 boundaries of the packing to the input packing
and stores in another, output file. It can be used to verify that periodic boundary conditions are properly applied or generate periodic packing image

Moreover, each mode accepts position-agnostic options to set parameters:
* `-f [input file]` - loads parameters from the input file specified. An example input file can be found in `polygonRSA/input.txt`
* `-i` - reads parameters from stdin. Note, that it has to receive EOF (`Ctrl+D`) to stop reading
* `-paramName paramValue` - can be used to change parameter values directly from the command line

Parameters can be read from multiple files and set using `-paramName paramValue` with the same keys - new values always override those specified earlier. All
unspecified parameters acquire their default values - see `polygonRSA/Parameters.h`.

### Output

Each simulation run produces `*.dat` file with 3 columns - packing number, number of shapes and dimensionless time of adding the last shape. Moreover, if requested
in the parameters, it can save each packing separately in `*.bin` binary file. The names of files encode most important parameters of the packings.

### Supported shapes

This version supports:

* **Polygon** - convex or concave polygon
* **Triangle** - convenient alias of **Polygon** for a triangle
* **RoundedPolygonAlgebraic** - convex or concave polygon with rounded corners using algebraic voxel elimination scheme
* **RegularRoundedPolygonAlgebraic** - **RoundedPolygonAlgebraic** optimised for regular polygons
* **RegularRoundedPolygonGeometric** - regular polygon with rounded corners using geometric voxel elimination scheme

The attibutes for each shape type can be found in code comments in shape classes from `polygonRSA/shape/shapes/polygon` directory.
