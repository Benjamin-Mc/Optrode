import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py

from matplotlib.backends.backend_pdf import PdfPages
from Tkinter import *
import tkFileDialog as filedialog

def integral_plot(intensities, wavelengths, export_pdf, times):

    intensities_2d = np.ndarray.reshape(intensities, (len(intensities)/len(wavelengths), len(wavelengths)))
    intensities_2d.transpose()
    i = 0
    min_index = 0
    while wavelengths[i] <= 550:
        if wavelengths[i] >= 505:
            if min_index == 0:
                min_index = i
        i += 1
    max_index = i
    area_arr = []
    subtotal = 0
    for i in range(len(intensities)/len(wavelengths)):
        for j in range(min_index, max_index):
            subtotal += intensities_2d[i, j]
        area_arr.append(subtotal)
        subtotal = 0
    plt.figure()
    plt.plot(times, area_arr)
    plt.title("Integral over time")
    plt.xlabel('time (s)')
    plt.ylabel('Sum of intensities of wavelengths between 505-550 nm')
    export_pdf.savefig()
    plt.close()


def dualplot(intensities, wavelengths, export_pdf):
        
    intensities_2d = np.ndarray.reshape(intensities, (len(intensities)/len(wavelengths), len(wavelengths)))
    intensities_2d.transpose()
    #2D Matrix, each column is a wavelengths and the rows are the intesities for that wavelength
    intensities_sum = np.ndarray.sum(intensities_2d, axis=0)
    intensities_avg = np.ndarray.mean(intensities_2d, axis=0)

    plt.figure(figsize=(12, 12))
    plt.subplot(2, 1, 1)
    plt.plot(wavelengths, intensities_sum)
    plt.title("Sum of Intensities")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')

    plt.subplot(2, 1, 2)       
    plt.plot(wavelengths, intensities_avg)
    plt.title("Average Intensity")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    export_pdf.savefig()
    plt.tight_layout(h_pad=10.0)
    plt.close()

def read_data(filename):
    data = h5py.File(filename, "r")        
    file_data = []
    #Reads groups and arrays in order -- Photodiode (Values, times), Powermeter (Values, times), Spectrometer (Values, times, wavelengths)
    for group in data.keys():
        for dset in data[group].keys() :
            file_data.append(data[group][dset][:]) # adding [:] returns a numpy array
    return file_data

def generate_output(sample_file, background_file, document_title):
    with PdfPages(r"/home/frederique/PhysicsLabPythonCodes/Optrode/Records/{}.pdf".format(document_title)) as export_pdf:
        sample_data = read_data(sample_file)
        background_data = read_data(background_file)
        
        wavelengths, intensities, times = sample_data[6], sample_data[4], sample_data[5]
        wavelengths_b, intensities_b, times_b = background_data[4], background_data[2], background_data[3]
        net_intensities = np.subtract(intensities, intensities_b)
        
        dualplot(intensities, wavelengths, export_pdf)
        dualplot(intensities_b, wavelengths, export_pdf)
        dualplot(net_intensities, wavelengths, export_pdf)

        integral_plot(intensities, wavelengths, export_pdf, times)
        integral_plot(net_intensities, wavelengths, export_pdf, times)            

    print("\\n Output generated")

def browse(n):
    filename = filedialog.askopenfilename(initialdir="/home/frederique/PhysicsLabPythonCodes/Optrode/Records/", filetypes=(("HDF5 files", "*.hdf5"),))
    if n==1:
        sample_file.set(filename)
    else:
        background_file.set(filename)

def close_program():
    root.destroy()

if __name__ == '__main__':

    global lower_bound, upper_bound
    lower_bound = 505
    upper_bound = 550

    root = Tk()
    root.title("Linear Flow Fluorescence Analysis")
    root.geometry("750x200")
    root.grid_columnconfigure(0, weight=1)

    document_title = StringVar()
    sample_file = StringVar()
    background_file = StringVar()

    frame1 = Frame(root)
    frame1.grid(row=0, column=0)

    title = Label(frame1, text="Linear Flow Fluorescence Analysis", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=2, column=0, padx=8, pady=6, sticky=W+E)

    ttl = Label(frame2, text="Enter document title", font=(None, 11))
    ttl.grid(row=1, column=1, padx=28, pady=6, sticky=W)

    fnlbl = Label(frame2, text="Select Sample Measurement", font=(None, 11))
    fnlbl.grid(row=2, column=1, padx=28, pady=6, sticky=W)

    xnlbl = Label(frame2, text="Select Background Measurement", font=(None, 11))
    xnlbl.grid(row=4, column=1, padx=28, pady=6, sticky=W)

    input_box0 = Entry(frame2, textvariable=document_title, width=40)
    input_box0.grid(row=1, column=2, padx=6)

    tag = Label(frame2, text=".pdf", font=(None, 10))
    tag.grid(row=1, column=3, sticky=W)

    input_box1 = Entry(frame2, textvariable=sample_file, width=40)
    input_box1.grid(row=2, column=2, padx=12)

    input_box2 = Entry(frame2, textvariable=background_file, width=40)
    input_box2.grid(row=4, column=2, padx=12)

    but1 = Button(frame2, text="Browse", command=lambda: browse(1))
    but1.grid(row=2, column=3)
    but2 = Button(frame2, text="Browse", command=lambda: browse(2))
    but2.grid(row=4, column=3)

    frame3 = Frame(root)
    frame3.grid(row=5, column=0)

    but3 = Button(frame3, text="Start", command=lambda: generate_output(sample_file.get(), background_file.get(), document_title.get()))
    but3.grid(row=5, column=0, padx=10, pady=10)

    but4 = Button(frame3, text="Close", command=close_program)
    but4.grid(row=5, column=1)

root.mainloop()
