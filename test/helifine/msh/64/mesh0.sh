#!/bin/bash

np=$1

gmsh wing0.geo -3 -oneFilePerPart -part $(($np)) -format msh2
