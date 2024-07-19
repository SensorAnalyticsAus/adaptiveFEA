## About
Finite Element Analysis (FEA) is a numerical modelling technique for predicting behaviour of solids, fluids, electromagnetic fields, and heat-flows. Adaptive FEA automatically adjusts the density of the mesh to provide more accurate results. It's done by calculating FEA error and using it to alter the mesh density until error is reduced to a tolerable limit. This archive contains four programs for adaptive FEA, **femgs1** and **adpgs1** perform FEA and error estimation respectively. Adaptive mesh generation is done by **faopgs2** . Output meshes with options to overlay error and stress data for visualisation is done by **plotp**.

## Prerequisites
Assuming `gcc` is already installed.
```
sudo apt update
sudo apt upgrade
sudo apt install gfortran
sudo apt install xpdf ghostscript
```

## Install
`cd` to this repo.
```
make distclean
make all
make clean
```

## Getting Started
Following steps will complete a cycle of adaptive FEA
# Plot the initial mesh
```
cd ./run
cp ../dot-dat/ex61.dat .
cp ex61.dat input.dat
```
Run the mesh plotting program to view the input mesh
```
../plotp
 output device is postscript file pltf*.plt
 Brightness controll factor [0.]...[1.]
0.7
 plot membrane elements ? (t/f)
t
 plot element errors or Domain dec? (t/f)
f
 plot stress ? (t/f)
f
 plot cable elements ?    (t/f)
f
 plot g strings ?   (t/f)
f
 number nodes ?           (t/f)
t
 number elements ? (t/f)
f
 element reduction factor ? (0.0 - 1.0)
0.
 auto scale ?             (t/f)
t
 plot single view ? t/f
t
 view no   1/2/3/4 ?
1
  isoyes =   F
    pltf1.plt     opened 
 enter the title
Initial Mesh ex61.dat
   4.6354998053636374E-310   0.0000000000000000        0.0000000000000000        0.0000000000000000
```
Convert output plot file to PDF format and view **ex61.dat**:
![ex61.png](/assets/png/ex61.png)
```
ps2pdf pltf1.plt
xpdf pltf1.plt.pdf
```
# FEA of the Initial Mesh
```
../femgs1 
Enter the input file (.dat assumed)
ex61
Example rect
0.001000 0.000000 
Guass [1]  JCC [0]
1

Semi-Band-Width before = 112

 Semi-Band-Width reduction Yes[1] No[0]
0

Semi-Band-Width after = 112

 F E analysis started with node = 96 and elements = 153
iln = 5
...
...
*** Time for analysis = 0.000000 min or 0 sec ***
*** no of iterations = 0 ***
*** The stresses are stored in ex61.str ***
*** The results of FE  are stored in ex61.fem ***
*** renumbered input data file ex6_.dat generated ***
```
# Adaptive error estimate
```
../adpgs1 
Enter the input file (.dat assumed)
ex61
Example rect
0.001000 0.000000 
opening ex61.str for reading the stresses
-5.27588e+00  2.28607e-01 -2.51780e-01
...
...
nita (percent) from the current mesh = 48.490119
 he_min = 0.871135  he_max = 106.717485
Enter your own hmin 
1.25
Magnification =   1.4 (hmin =   1.2)
*** Element errors are stored in ex61.err ***
*** Mesh parameters are in ex61.me 
*** Adaptivity results are in ex61.adp
```
# Re-mesh with adaptive refinement this time
```
../faopgs2 
Max of 50000 nodes and 70000 elements can be meshed

***Node para are going to be used*** 
Enter the project file name (without extension)
ex61
Example rect
0.001000 0.000000 
opening ex61.men for reading 
cannot open ex61.men  opening ex61.me for read
Enter [1] if coarse background mesh is to be refined
1
Enter [1] if mesh post-processing required
1
...
...
*** Mesh Compiled with 1600 Nodes and 3008 Element *** 
*** File ex61o.dat with 1600 nodes & 3008 elements has been generated ***
*** Neural net training data saved in ex61.inf ***
```
# Perform FEA on newly created mesh in ex61o.dat
```
../femgs1
Enter the input file (.dat assumed)
ex61o
Example rect
0.001000 0.000000 
Guass [1]  JCC [0]
1

Semi-Band-Width before = 2620

 Semi-Band-Width reduction Yes[1] No[0]
1

...
...
*** Time for analysis = 0.000000 min or 0 sec ***
*** no of iterations = 0 ***
*** The stresses are stored in ex61o.str ***
*** The results of FE  are stored in ex61o.fem ***
*** renumbered input data file ex61_.dat generated ***
```
# Calculate errors for ex61o.dat mesh
```
../adpgs1 
Enter the input file (.dat assumed)
ex61o
Example rect
0.001000 0.000000 
opening ex61o.str for reading the stresses
 1.61877e+00 -3.04040e-02 -4.70029e-01
...
...
nita (percent) from the current mesh = 32.172155
 he_min = 0.017925  he_max = 1371.344995
Enter your own hmin 
.9
Magnification =  50.2 (hmin =   0.9)
*** Element errors are stored in ex61o.err ***
*** Mesh parameters are in ex61o.me 
*** Adaptivity results are in ex61o.adp
```
# Copy output mesh, stress, and error files to plotp filenames
```
cp ex62o.dat input.dat
cp ex62o.str stress.dat
cp ex61o.err error.dat
```
# Plot the new mesh with say the Sxy (shear stress) overlay
```
../plotp 
 output device is postscript file pltf*.plt
 Brightness controll factor [0.]...[1.]
0.7
 plot membrane elements ? (t/f)
t
 plot element errors or Domain dec? (t/f)
f
 plot stress ? (t/f)
t
 plot cable elements ?    (t/f)
f
 plot g strings ?   (t/f)
f
 number nodes ?           (t/f)
f
 number elements ? (t/f)
f
 element reduction factor ? (0.0 - 1.0)
0.
 auto scale ?             (t/f)
t
 enter 1 2 3 4 5 for sx sy sxy smax smin
3
 plot single view ? t/f
t
 view no   1/2/3/4 ?
1
  isoyes =   F
    pltf1.plt     opened 
 enter the title
Refined Mesh ex61o.dat Sxy
   4.6354998053636374E-310   0.0000000000000000       -22.795900000000000        33.677500000000002   
   
```
Convert to PDF and view **ex61o.dat**:
![ex61o.png](/assets/png/ex61o.png)
```
ps2pdf pltf1.plt
xpdf pltf1.plt.pdf
```

The last step can be repeated to obtain a plot with the error overlay by setting plot element errors or Domain dec? option as 't' and the following option as 'f'.

## dot-dat\ Contains Data Files for Adaptive Finite Element Analysis
A *dot-dat* file specifies the starting mesh topology, material properties of the structure, boundary conditions, and applied force(s) i.e. load points on the structure. In this regard `try.pdf` contains information for creating custom *dot-dat* files. In such files the columns have to be maintained exactly in the same format as of `try.dat`. For instance the 3rd line of the try *dot-dat* file comprises first two fields of 6 columns i.e. 2 x 6 columns, followed by 9 fields of 5 columns i.e. 9 x 5 columns (a total of 57 vim columns). Each value must be fully contained within its own column. The FORTRAN graphics program `plotp` expect input values strictly being within the specified columns. 

Other *dot-dat* files may contain some data fields outside of `try.dat` populated spaces, such data fields are not read thus can be ignored.
## Further Information
[1]: PARALLEL FINITE ELEMENT COMPUTATIONS B.H.V. Topping and A.I. Khan [Chapter 5] (https://www.saxe-coburg.co.uk/pubs/descrip/btak.htm)

