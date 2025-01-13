#!/bin/bash

g++ 1d_example.cpp linear.cpp cubic.cpp utils.cpp -lm -o ex1d.exe
./ex1d.exe; rm ex1d.exe

python3 1d_analysis.py