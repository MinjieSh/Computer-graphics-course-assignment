README
======

This package includes the following files.

|-- raytracer.cc [Read the driver file, apply transformation, ray-triangle/sphere intersection, 
	illumination, shading, and reflection then write the ppm file.]
|-- Camera.cc [Camera class file]
|-- Camera.h [Camera header file]
|-- Light.cc [Light class file]
|-- Light.h [Light header file]
|-- Model.cc [Model class file]
|-- Model.h [Model header file]
|-- Sphere.cc [Sphere class file]
|-- Sphere.h [Sphere header file]

|-- driver01.txt [My driver file]
|-- driver01.ppm [My result for my driver file]
|-- driver02.txt [My driver file]
|-- driver02.ppm [My result for my driver file]
|-- driver03.ppm [My result for my driver file]
|-- Makefile [A Makefile that performs both a make clean as well as a make.]
|-- README.txt [This file]
|-- eigen-eigen-b3f3d4950030 [A folder contains the Eigen library.]

To compile:
    make

To run:
    ./raytracer driver##.txt driver##.ppm

For example;
    ./raytracer driver00.txt driver00.ppm
    

To clean:
    make clean
