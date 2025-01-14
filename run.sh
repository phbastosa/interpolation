#!/bin/bash

g++ 1d_example.cpp linear.cpp cubic.cpp utils.cpp -lm -o ex1d.exe
./ex1d.exe

python3 1d_analysis.py
rm *.bin *.exe

g++ 2d_example.cpp linear.cpp cubic.cpp utils.cpp -lm -o ex2d.exe
./ex2d.exe

python3 2d_analysis.py
rm *.bin *.exe