To perform a test, execute the "OptrodeVersion5.py" file.
To do this in the linux OS:
- Open a terminal
- Change to the directory that contains the program -- "cd Documents/Olivier/LFF"
- Execute the program -- "python OptrodeVersion5.py"


Potential Complications:
- At times the program needs permission to access the devices -- To resolve this, enter the sudo password
- Errors when reading the DAQ card -- To resolve, simple exit the program and re-run it


Output:
This program gives three plots of output, these are:
1. 
2. 
3. 


In the main program, there are multiple other programs, as well as libraries that have been imported for use. These are:

Programs:

"DAQT7_Objective.py"
-- This read the information from the DAQ card.

"SeaBreeze_Objective.py"
-- This read the information from the Spectrometer.

"ThorlabsPM100_Objective.py"
-- This read the information from the Power Meter.


Libraries/Modules used in this program:

"import sys"
-- Used for interacting with the interpreter and command line

import socket
import subprocess
import struct

"import os"
"import os.path"
"import glob"
-- Used for creating the path for output data to the Records folder

import tempfile

"import datetime"
-- Used for timestamping the output data with the date and time

"import time"
-- Used for to keep track of time while performing the experiment

"import bisect"
-- Used for finding the index of an element that is not currently in a sorted list, if it were to be added to the list.
bisect.bisect(a, x, lo=0, hi=len(a)) -- returns an insertion index, i, of a list a, such that each element to the left of i is <= x, and each element to the right of i is > x. Lo and hi can be used to search shorter sublists, their default values search the whole list.

from multiprocessing import Process, Value, Array

"from Tkinter import *"
"from ttk import Button, Style, Label, Entry, Notebook, Scale"
"from tkFileDialog import askopenfilename"
"from PIL import Image, ImageTk"
-- Used for creating the GUI

"import h5py"
-- Used for reading HDF5 binary data

import DAQT7_Objective as DAQ
import SeaBreeze_Objective as SBO
import ThorlabsPM100_Objective as P100

"import numpy as np"
-- Useful for creating arrays/matrices of 0's
np.zeros(shape=(i, j), dtype="float") -- creates an ixj matrix of zeroes, dtype can be set to integers or floating point values.

"import matplotlib.pyplot as plt"
-- Used for creating the plots

