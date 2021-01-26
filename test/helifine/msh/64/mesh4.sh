#!/bin/bash

np=$1

gmsh fuspyl.geo -3 -oneFilePerPart -part $(($np)) -format msh2
