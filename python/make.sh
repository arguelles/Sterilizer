#!/bin/sh
rm build/temp.linux-x86_64-2.7/SterileSearchPy.o 
rm build/lib.linux-x86_64-2.7/SterileSearchPy.so 
python setup.py build_ext