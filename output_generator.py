#Output Script

from reportlab.lib.units import inch, cm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py

def hello(c):
    c.drawString(100, 100, "Hello World")

if __name__ == "__main__":

    c = canvas.Canvas("Optrode Measurement data.pdf")
    hello(c)
    c.showPage()
    c.save()

    path = "./Records/'test2'-2020-12-07-14-04-11.hdf5"
    data = h5py.File(path, 'r')

    data_list = []
    #Reads groups and arrays in order -- Photodiode (Values, times), Powermeter (Values, times), Spectrometer (Values, times, wavelengths)
    for group in data.keys():
        for dset in data[group].keys() :
            data_list.append(data[group][dset][:]) # adding [:] returns a numpy array
    
    wavelengths = data_list[6]
    intensities = data_list[4]
    times = data_list[5]
    intensities_2d = np.ndarray.reshape(intensities, (len(intensities)/len(wavelengths), len(wavelengths)))
    intensities_2d.transpose()
    #2D Matrix, each column is a wavelengths and the rows are the intesities for that wavelength
    intensities_sum = np.ndarray.sum(intensities_2d, axis=0)
    intensities_avg = np.ndarray.mean(intensities_2d, axis=0)

    plt.figure()
    for i  in range(len(intensities)/len(wavelengths)):
        plt.plot(wavelengths, intensities_2d[i, :])
    plt.title("Spectrometer Recordings")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.show()

    plt.figure()  
    plt.plot(wavelengths, intensities_sum)
    plt.title("Sum of Intensities")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.show()

    plt.figure()
    plt.plot(wavelengths, intensities_avg)
    plt.title("Average Intensity")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.show()

    integral_array = []
    wavelength_array = []
    time_array = []
    area_under_curve = 0 
    for i in range(len(wavelengths)):
        if wavelengths[i] >= 505 and wavelengths[i] <= 550:
            area_under_curve += intensities_avg[i]
            wavelength_array.append(wavelengths[i])
            integral_array.append(area_under_curve)

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
    for i in range(1250):
        for j in range(min_index, max_index):
            subtotal += intensities_2d[i, j]
        area_arr.append(subtotal)
        subtotal = 0

    plt.figure()
    plt.plot(wavelength_array, integral_array)
    plt.title("Integral area")
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Area under the curve')
    plt.show()

    plt.figure()
    plt.plot(times, area_arr)
    plt.title("Integral over time")
    plt.xlabel('time (s)')
    plt.ylabel('Sum of intensities of wavelengths between 505-550 nm')
    plt.show()

