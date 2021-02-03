'''
Shorter program than "Compare_PM-PD.py"
This program plots avg pm/pd for each sample on one plot to compare (more than 6 samples).
Author: Benjamin McIntosh 
Date: 03/02/2021
'''
import os
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
from Tkinter import *
import tkFileDialog as filedialog

def read_data(filename):
    '''
    Reads the sample files and returns the data as a list
    '''
    data = h5py.File(filename, "r")        
    file_data = []
    #Reads groups and arrays in order -- Photodiode (Values, times), Powermeter (Values, times), Spectrometer (Values, times, wavelengths)
    for group in data.keys():
        for dataset in data[group].keys():
            file_data.append(data[group][dataset][:]) # adding [:] returns a numpy array
    return file_data

def browse(n):
    '''
    Allows user to browse directories for input data files
    '''
    filenames = filedialog.askopenfilenames(initialdir="/home/frederique/PhysicsLabPythonCodes/Optrode/Records/", filetypes=(("HDF5 files", "*.hdf5"),))
    sample_files.set(filenames)
    global file_list, file_titles
    file_list = sample_files.get().split(", ")
    file_titles = []
    trim(file_list, file_titles, num_samples.get())

def trim(f_list, f_titles, n):
    '''
    Adjust the names of the input parameters so that the program can read the correct files
    Adds the filenames to a list of titles for the plot labels
    
    For example, will trim ('/home/frederique/PhysicsLabPythonCodes/Optrode/Records/Green Cont - Saline test 1.hdf5',) to:
    filename: "/home/frederique/PhysicsLabPythonCodes/Optrode/Records/Green Cont - Saline test 1.hdf5"
    filetitle: "Green Cont - Saline test 1"
    '''
    for i in range(n):
        start = f_list[i].find("/")
        end = f_list[i].rfind("5")
        f_list[i] = f_list[i][start:end+1]
        start2 = f_list[i].rfind("/")
        end2 = f_list[i].rfind(".")
        f_titles.append(f_list[i][start2+1:end2])

def close_program():
    root.destroy()
    exit()

def generate_output(files, titles, img_title):
    '''
    Takes the input data and generates:
    - A plot for each sample comparing the Powermeter and Photodiode (Photodiode should be in blue and Powermeter in green)
    - A scatter plot where each point is the the average Photodiode reading (x) and average Powermeter reading (y)
    The first point on the plot will be the same sample as the first plot of the PD/PM comparison plots
    '''
    data = []
    for i in range(num_samples.get()):
        data.append(read_data(files[i]))

    Path_to_Records = "/home/frederique/PhysicsLabPythonCodes/Optrode/Records/"

    pdr, pdt, pmr, pmt = [], [], [], []
    for k in range(num_samples.get()):
        #Subtract first value to get scaled time
        pd_readings = np.asarray(data[k][0])
        pd_times = np.asarray(data[k][1])
        pd_times -= pd_times[0]
        pm_readings = np.asarray(data[k][2])
        pm_times = np.asarray(data[k][3])
        pm_times -= pm_times[0]

        #Cut off trailing zeroes
        for i in range(len(pd_readings)):
            if pd_readings[i] == float(0):
                pd_index = i
                break
        for j in range(len(pm_readings)):
            if pm_readings[j] == float(0):
                pm_index = j
                break

        pdr.append(pd_readings[:i])
        pdt.append(pd_times[:i])
        pmr.append(pm_readings[:j])
        pmt.append(pm_times[:j])

        #adjust scale
        pd_readings = pd_readings / 1000

    pd_avg, pm_avg = [], []
    for i in range(num_samples.get()):
        pd_readings = pdr[i]
        pm_readings = pmr[i]
        pd_avg.append(sum(pd_readings)/len(pd_readings))
        pm_avg.append(sum(pm_readings)/len(pm_readings))

    fig1 = plt.figure()
    plt.title("Average Photodiode and Powermeter")
    for i in range(num_samples.get()):
        plt.plot(pd_avg[i], pm_avg[i], "o", alpha=0.8, label=titles[i])
    plt.ylabel("Powermeter (mW)")
    plt.xlabel("Photodiode (V)")
    fig1.savefig('{}/{} - Graph.png'.format(Path_to_Records, img_title))

    fig2 = plt.figure()
    for i in range(num_samples.get()):
        plt.plot(0, 0, "o", label=titles[i])
    plt.legend(loc="center")
    plt.axis("off")
    fig2.savefig('{}/{} - Legend.png'.format(Path_to_Records, img_title))
    plt.show()

if __name__ == "__main__":
    '''
    This is the code for the GUI
    '''
    root = Tk()
    root.title("Average PM/PD")
    root.geometry("500x175")
    root.grid_columnconfigure(0, weight=1)

    sample_files = StringVar(value="")
    image_title = StringVar(value="")

    # Title frame
    frame1 = Frame(root)
    frame1.grid(row=0, column=0)

    title = Label(frame1, text="Average PM/PD", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    # Parameter frame
    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=1, column=0, padx=8, pady=4)

    frame3 = Frame(root)
    frame3.grid(row=2, column=0)

    global num_samples
    num_samples = IntVar(value=1)   
    numb = Entry(frame2, textvariable=num_samples, width=3)
    numb.grid(row=2, column=2, pady=5, padx=5, sticky=W)

    ttl = Label(frame2, text="Enter image title: ", font=(None, 11))
    ttl.grid(row=0, column=0, pady=5)

    input_box0 = Entry(frame2, textvariable=image_title, width=25)
    input_box0.grid(row=0, column=1, padx=5)

    input_box1 = Entry(frame2, textvariable=sample_files, width=25)
    input_box1.grid(row=2, column=0, padx=5)

    but1 = Button(frame2, text="Browse", command=lambda: browse(1))
    but1.grid(row=2, column=1)

    but3 = Button(frame3, text="Start", command=lambda: generate_output(file_list, file_titles, image_title.get()))
    but3.grid(row=5, column=0, padx=10, pady=10)

    but4 = Button(frame3, text="Close", command=lambda: close_program())
    but4.grid(row=5, column=1)

    root.mainloop()


