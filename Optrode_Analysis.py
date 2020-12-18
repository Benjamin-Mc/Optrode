"""
Date: 18/12/2020
Author: Benjamin McIntosh
This program generates a PDF with analysis plots for the Linear Flow Optrode system
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py

from matplotlib.backends.backend_pdf import PdfPages
from Tkinter import *
import tkFileDialog as filedialog


def integral_plot(intensities, wavelengths, export_pdf, times, subtitle):
    '''
    Plots the integral over time
    '''
    integrals = []
    for j in range(len(intensities)):
        intensities_2d = np.ndarray.reshape(intensities[j], (len(intensities[j])/len(wavelengths[j]), len(wavelengths[j])))
        intensities_2d.transpose()
        i = 0
        min_index = 0
        while wavelengths[j][i] <= upper_bound:
            if wavelengths[j][i] >= lower_bound:
                if min_index == 0:
                    min_index = i
            i += 1
        max_index = i
        area_arr = []
        subtotal = 0
        for i in range(len(intensities[j])/len(wavelengths[j])):
            for k in range(min_index, max_index):
                subtotal += intensities_2d[i, k]
            area_arr.append(subtotal)
            subtotal = 0
        integrals.append(area_arr)

    for f in range(len(integrals)):
        plt.plot(times[f], integrals[f], label=file_titles[f])
        plt.title("{} Integral over time".format(subtitle))
        plt.xlabel('time (s)')
        plt.ylabel('Sum of intensities of wavelengths between {}-{} nm (a.u.)'.format(lower_bound, upper_bound))

def dualplot(intensities, wavelengths, export_pdf, page_heading, background):
    '''
    Plots the sum of intensities and average intensity against the wavelength range
    '''
    #There may be multiple input files, background file always has one sample
    if background == True:
        num_samples = 1
    else:
        num_samples = num_files.get()
    

    sum_intensity, avg_intensity = [], []
    #Gives a list of lists for each sample
    for i in range(num_samples):
        intensity = intensities[i]
        intensities_2d = np.ndarray.reshape(intensity, (len(intensity)/len(wavelengths[i]), len(wavelengths[i])))
        intensities_2d.transpose()
        #2D Matrix, each column is a wavelength and the rows are the intesities for that wavelength
        intensities_sum = np.ndarray.sum(intensities_2d, axis=0)
        intensities_avg = np.ndarray.mean(intensities_2d, axis=0)
        sum_intensity.append(intensities_sum)
        avg_intensity.append(intensities_avg)

    #Plots sum of intensites against wavelengths
    plt.figure(figsize=(12, 12), dpi=1200)
    plt.suptitle(page_heading, fontsize=24) 
    plt.axis("off")
    plt.subplot(3, 1, 1) 
    for i in range(num_samples):
        plt.plot(wavelengths[i], sum_intensity[i], label=file_titles[i])
    plt.title("Sum of Intensities")
    #plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend()

    #Plots average intensity against wavelengths
    plt.subplot(3, 1, 2)
    for i in range(num_samples):    
        plt.plot(wavelengths[i], avg_intensity[i], label=file_titles[i])
    plt.title("Average Intensity")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend()

    #Writes some output statistics for the plots
    statistics = ""
    for i in range(num_samples):
        peak_sum = max(sum_intensity[i])
        peak_avg = max(avg_intensity[i])
        avg_across_wavelengths = sum(avg_intensity[i])/float(len(avg_intensity[i]))
        statistics += file_titles[i]
        statistics += ("Peak Value Sum of Intensities: {} a.u. \n".format(peak_sum))
        statistics += ("Peak Value Average Intensity: {} a.u. \n".format(peak_avg))
        statistics += ("Average Intensity across all Wavelengths: {} a.u. \n".format(avg_across_wavelengths))
        statistics += ("\n")
    plt.figtext(0.1, 0.2, statistics)

    #Saves the plots and data to the PDF
    export_pdf.savefig()
    plt.close()

def read_data(filename):
    '''
    Reads the sample files and returns the data as a list
    '''
    data = h5py.File(filename, "r")        
    file_data = []
    #Reads groups and arrays in order -- Photodiode (Values, times), Powermeter (Values, times), Spectrometer (Values, times, wavelengths)
    for group in data.keys():
        for dset in data[group].keys():
            file_data.append(data[group][dset][:]) # adding [:] returns a numpy array
    return file_data

def generate_output(sample_file, background_file, document_title):
    '''
    Main function for generating the PDF
    '''
    with PdfPages(r"/home/frederique/PhysicsLabPythonCodes/Optrode/Records/{}.pdf".format(document_title)) as export_pdf:

        #Title page
        plt.figure()
        plt.axis("off")
        plt.text(0.5, 0.5, "Linear Flow Fluorescence", fontsize=28, ha="center", va="center")
        export_pdf.savefig()
        plt.close()        
        
        background_data = read_data(background_file)
        intensities_b, times_b, wavelengths_b = [background_data[2]], [background_data[3]], [background_data[4]]

        intensities, times, wavelengths, net_intensities = [], [], [], []
        for i in range(num_files.get()):
            sample_data = read_data(sample_file[i]) 
            #List of lists with the data in order: Photodiode readings and times (0, 1), Power meter readings and times (2, 3), Spectrometer intensities, times and wavelengths (4, 5, 6)
            adjust = 2
            if len(sample_data) == 7:
                adjust = 0
            intensities.append(sample_data[4-adjust])
            times.append(sample_data[5-adjust])
            wavelengths.append(sample_data[6-adjust])
            #net_intensity = np.subtract(sample_data[4-adjust], intensities_b[0])
            #net_intensities.append(net_intensity)
            #photodiodes, times_d = sample_data[0], sample_data[1]
            #pm, times_p = sample_data[2], sample_data[3]

            '''
            plt.figure()
            plt.plot(times_d, photodiodes)
            plt.xlabel("Time (s)")
            plt.ylabel("PD")
            export_pdf.savefig()
            plt.close()

            plt.figure()
            plt.plot(times_p, pm)
            plt.xlabel("Time (s)")
            plt.ylabel("PM")
            export_pdf.savefig()
            plt.close()
            '''

            #Completes plots
        dualplot(intensities, wavelengths, export_pdf, "Sample Reading", False)
        dualplot(intensities_b, wavelengths_b, export_pdf, "Backgroud Reading", True)
        #dualplot(net_intensities, wavelengths, export_pdf, "Net Reading", False)

        plt.figure(figsize=(12, 12))
        plt.suptitle("Integrals", fontsize=24)
        plt.subplot(2, 1, 1)
        integral_plot(intensities, wavelengths, export_pdf, times, "Sample")
        plt.legend()
        #plt.subplot(2, 1, 2)
        #integral_plot(net_intensities, wavelengths, export_pdf, times, "Net")
        #plt.legend()            
        export_pdf.savefig()
        plt.close()

    print(" \nOutput generated")

def trim(file_list):
    '''
    Adjust the names of the input parameters so that the program can read the correct files
    Adds the filenames to a list of titles to label the plots
    '''
    for i in range(num_files.get()):
        start = file_list[i].find("/")
        end = file_list[i].rfind("5")
        file_list[i] = file_list[i][start:end+1]
        start2 = file_list[i].rfind("/")
        end2 = file_list[i].rfind(".")
        file_titles.append(file_list[i][start2+1:end2])

def browse(n):
    '''
    Allows user to browse directories for input data files
    '''
    filenames = filedialog.askopenfilenames(initialdir="/home/frederique/PhysicsLabPythonCodes/Optrode/Records/", filetypes=(("HDF5 files", "*.hdf5"),))
    if n==1:
        sample_files.set(filenames)
        global file_list, file_titles
        file_list = sample_files.get().split(", ")
        file_titles = []
        trim(file_list)
    else:
        background_file.set(filename)

def close_program():
    root.destroy()

if __name__ == '__main__':

    #Desired bounds for the integrals of the spectra
    global lower_bound, upper_bound
    lower_bound = 505
    upper_bound = 550

    root = Tk()
    root.title("Linear Flow Fluorescence Analysis")
    root.geometry("750x225")
    root.grid_columnconfigure(0, weight=1)

    document_title = StringVar(value="Test Optrode Analysis2")
    sample_files = StringVar(value="")
    background_file = StringVar(value="/home/frederique/PhysicsLabPythonCodes/Optrode/Records/bg.hdf5")
   

    #The following details all of the instructions for the Tkinter GUI

    frame1 = Frame(root)
    frame1.grid(row=0, column=0)

    title = Label(frame1, text="Linear Flow Fluorescence Analysis", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=2, column=0, padx=8, pady=6, sticky=W+E)

    global num_files
    num_files = IntVar(value=1)
    nums = [1, 2, 3, 4, 5]
    drop = OptionMenu(frame2, num_files, *nums)
    drop.grid(row=1, column=2)

    ynlbl = Label(frame2, text="Select number of Samples: ", font=(None, 11))
    ynlbl.grid(row=1, column=1, padx=28, pady=6, stick=W)

    ttl = Label(frame2, text="Enter document title: ", font=(None, 11))
    ttl.grid(row=0, column=1, padx=28, pady=6, sticky=W)

    fnlbl = Label(frame2, text="Select Sample Measurement(s): ", font=(None, 11))
    fnlbl.grid(row=2, column=1, padx=28, pady=6, sticky=W)

    xnlbl = Label(frame2, text="Select Background Measurement: ", font=(None, 11))
    xnlbl.grid(row=4, column=1, padx=28, pady=6, sticky=W)

    input_box0 = Entry(frame2, textvariable=document_title, width=40)
    input_box0.grid(row=0, column=2, padx=6)

    tag = Label(frame2, text=".pdf", font=(None, 10))
    tag.grid(row=0, column=3, sticky=W)

    input_box1 = Entry(frame2, textvariable=sample_files, width=40)
    input_box1.grid(row=2, column=2, padx=12)

    input_box2 = Entry(frame2, textvariable=background_file, width=40)
    input_box2.grid(row=4, column=2, padx=12)

    but1 = Button(frame2, text="Browse", command=lambda: browse(1))
    but1.grid(row=2, column=3)
    but2 = Button(frame2, text="Browse", command=lambda: browse(2))
    but2.grid(row=4, column=3)

    frame3 = Frame(root)
    frame3.grid(row=5, column=0)

    but3 = Button(frame3, text="Start", command=lambda: generate_output(file_list, background_file.get(), document_title.get()))
    but3.grid(row=5, column=0, padx=10, pady=10)

    but4 = Button(frame3, text="Close", command=close_program)
    but4.grid(row=5, column=1)

    root.mainloop()
