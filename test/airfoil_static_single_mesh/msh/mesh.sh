#!/bin/bash

np=$1

gmsh NACA0012_C.geo -3 -oneFilePerPart -part $(($np)) -format msh2
