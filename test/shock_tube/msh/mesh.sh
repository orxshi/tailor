#!/bin/bash

np=$1

gmsh sod.geo -3 -oneFilePerPart -part $(($np)) -format msh2
