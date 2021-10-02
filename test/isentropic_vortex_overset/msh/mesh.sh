#!/bin/bash

np=$1

gmsh bg.geo -3 -oneFilePerPart -part $(($np)) -format msh2
gmsh overset_square_a.geo -3 -oneFilePerPart -part $(($np)) -format msh2
