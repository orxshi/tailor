#!/bin/bash

np=$1

gmsh wing3.geo -3 -oneFilePerPart -part $(($np)) -format msh2
