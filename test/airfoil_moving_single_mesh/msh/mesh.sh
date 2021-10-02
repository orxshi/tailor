#!/bin/bash

np=$1

gmsh NACA0012_O.geo -3 -oneFilePerPart -part $(($np)) -format msh2
gmsh bg.geo -3 -oneFilePerPart -part $(($np)) -format msh2
