#!/bin/bash

#python main.py -a 5,6.0,1.0,100 -p 0.01 --multiprocessing=3
python main.py -a 5,6.0,-1.0,100 -p 0.01 --multiprocessing=3
python main.py -a 5,6.0,0.0,100 -p 0.01 --multiprocessing=3
python main.py -a 5,-6.0,1.0,100 -p 0.01 --multiprocessing=3
python main.py -a 5,-6.0,-1.0,100 -p 0.01 --multiprocessing=3
python main.py -a 5,-6.0,0.0,100 -p 0.01 --multiprocessing=3
