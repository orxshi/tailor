#!/bin/bash

np=$1

gmsh hubshaft.geo -3 -oneFilePerPart -part $(($np)) -format msh2
