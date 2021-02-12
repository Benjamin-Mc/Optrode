"""
Date: 18/12/2020
Author: Benjamin McIntosh
This program generates a PDF with analysis plots for the Linear Flow Optrode system
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
from matplotlib.backends.backend_pdf import PdfPages
from Tkinter import *
import tkFileDialog as filedialog

def integral_plot(intensities, wavelengths, export_pdf, times, subtitle):
    '''
    Plots the integral over time for a given dataset
    '''
    integrals = []
    #Finds the intensities that are within the specified wavelength range
    for j in range(len(intensities)):
        #Find index range
        intensities_2d = intensities[j]
        i, min_index = 0, 0
        while wavelengths[j][i] <= integral_upper_bound.get():
            if wavelengths[j][i] >= integral_lower_bound.get():
                if min_index == 0:
                    min_index = i
            i += 1
        max_index = i

        #Sum the intensities in the range
        area_arr = []
        for i in range(len(intensities_2d[0][:])):
            subtotal = 0
            for k in range(min_index, max_index):
                subtotal += intensities_2d[k, i]
            area_arr.append(subtotal)
            subtotal = 0
        integrals.append(area_arr)

    #Plots the integral for each sample
    for i in range(len(integrals)):
        plt.plot(times[i], integrals[i], label=file_titles[i])
        plt.title("{} Integral over time between {}-{} nm".format(subtitle, integral_lower_bound.get(), integral_upper_bound.get()))

def photodiode_plots(photodiodes, times, export_pdf, titles):
    '''
    Plots the photodiode readings for each sample
    '''
    fig = plt.figure()
    for i in range(len(photodiodes)):
        ls = ["-","--","-."][i%3]
        plt.plot(times[i], photodiodes[i], label=titles[i], linestyle=ls)
    plt.title("Photodiode Readings")
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")
    export_pdf.savefig()
    plt.close()

    #Plots the legend on a separate page, otherwise it usually cuts off
    fig = plt.figure()
    for i in range(num_files.get()):
        plt.plot(0, 0, "o", label=titles[i])
    plt.legend(loc="center")
    plt.axis("off")
    export_pdf.savefig()
    plt.close()

def powermeter_plots(power_readings, pm_times, export_pdf, title):
    '''
    Plots the powermeter readings for the background file
    '''
    plt.figure()
    plt.plot(pm_times[0], power_readings[0], label=title)
    plt.title("Powermeter Readings")
    plt.xlabel("Time (s)")
    plt.ylabel("Power (W)")
    export_pdf.savefig()
    plt.close()
    
def mainplots(intensities, wavelengths, export_pdf, page_heading, background, titles):
    '''
    Plots the sum of intensities and average intensity against the wavelength range
    If the paradigm type is continuous - plots every 100th spectrum for each sample
    If the paradigm type is multi integration - plots each spectrum normalised for the integration time (might look strange due to saturated reading)
    Calculates some output statistics for the plots, saves in a table
    '''
    #There may be multiple input files, background file always has one sample
    if background == True:
        num_samples = 1
    else:
        num_samples = num_files.get()
    
    sum_intensities, avg_intensities, selected_matrix, intensities_matrix = [], [], [], []
    #Gives a list of lists for each sample
    #Gives intensities as a 2D Matrix, each column is a wavelength and the rows are the intesities for that wavelength
    for i in range(num_samples):
        intensities_2d = np.asarray(intensities[i])
        #intensities_2d = intensities_2d.transpose()

        #Adds each 100th spectrum to a list
        j = 0
        selected_values = []
        while j < (len(intensities_2d[0, :])):
            selected_values.append(intensities_2d[:, j])
            j += 100
        selected_matrix.append(selected_values)
        count = len(selected_values)        

        intensities_matrix.append(intensities_2d)

        intensities_sum = np.ndarray.sum(intensities_2d, axis=1)
        intensities_avg = np.ndarray.mean(intensities_2d, axis=1)
        sum_intensities.append(intensities_sum)
        avg_intensities.append(intensities_avg)

    #If the samples are of a multi-integration paradigm
    if paradigm_mode.get() == "m":
        num_multi_integration = len(intensities_matrix[0][0, :])
        integration_times = [] #Integration time for the spectrometer in ms
        for i in range(num_multi_integration):
            time = 2**(i+3) #Adds integration times 8, 16, 32, 64, 128... to the list of integration times
            integration_times.append(time)
        if num_samples == 1:
            plt.figure(dpi=1200)
            plt.suptitle("{}\nNormalised Spectra".format(page_heading), fontsize=9)
            for j in range(len(integration_times)):
                intensities = intensities_matrix[0][:, j]
                normal_intensities = intensities / float(integration_times[j])
                ls = ["-", "--", "-."][j%3]
                plt.plot(wavelengths[0], normal_intensities, linestyle=ls)
                plt.title(titles[0], fontsize=9)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Intensity (a.u.)")

        else: #More that one sample
            #Creates a grid of plots, one for each sample
            num_cols = 2
            num_rows = num_samples // num_cols
            num_rows += num_samples % num_cols
            pos = range(1, num_samples+1) 
            fig = plt.figure(dpi=1200)
            fig.suptitle("{}\nNormalised Spectra".format(page_heading), fontsize=9)
            for k in range(num_samples):
                ax = fig.add_subplot(num_rows, num_cols, pos[k])
                for j in range(len(integration_times)):
                    intensities = intensities_matrix[k][:, j]
                    normal_intensities = intensities / float(integration_times[j])
                    ls = ["-","--","-."][j%3]
                    ax.plot(wavelengths[k], normal_intensities, linestyle=ls)
                ax.title.set_text(titles[k])
            for ax in fig.get_axes():
                ax.label_outer()
            fig.text(0.5, 0.04, "Wavelength (nm)", ha="center", fontsize=14)
            fig.text(0.04, 0.5, "Intensity (a.u.)", ha="center", rotation="vertical", fontsize=14)
        export_pdf.savefig()
        plt.close()

    #If the samples are of a continuous paradigm
    if paradigm_mode.get() == "c":
        if num_samples == 1:
            plt.figure(dpi=1200)
            plt.suptitle("{}\nSelected Spectra".format(page_heading), fontsize=9)
            for j in range(count):
                ls = ["-", "--", "-."][j%3]
                plt.plot(wavelengths[0], selected_matrix[0][j], linestyle=ls)
                plt.title(titles[0], fontsize=9)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Intensity (a.u.)")
        else: #More than one sample
            #Creates a grid of plots, one for each sample
            num_cols = 2
            num_rows = num_samples // num_cols
            num_rows += num_samples % num_cols
            pos = range(1, num_samples+1) 
            fig = plt.figure(dpi=1200)
            for k in range(num_samples):
                ax = fig.add_subplot(num_rows, num_cols, pos[k])
                for j in range(count):
                    ls = ["-","--","-."][j%3]
                    ax.plot(wavelengths[k], selected_matrix[k][j], linestyle=ls)
                ax.title.set_text(titles[k])
            fig.suptitle("{}\nSelected Spectra".format(page_heading), fontsize=9)
            for ax in fig.get_axes():
                ax.label_outer()
            fig.text(0.5, 0.04, "Wavelength (nm)", ha="center", fontsize=14)
            fig.text(0.04, 0.5, "Intensity (a.u.)", ha="center", rotation="vertical", fontsize=14)
        export_pdf.savefig()
        plt.close()

    #Plots sum of intensites against wavelengths
    fig = plt.figure(figsize=(12, 12), dpi=1200)
    fig.text(0.5, 0.04, "Wavelength (nm)", ha="center", fontsize=14)
    fig.text(0.04, 0.5, "Intensity (a.u.)", ha="center", rotation="vertical", fontsize=14)
    plt.suptitle(page_heading, fontsize=24) 
    plt.axis("off")
    plt.subplot(2, 1, 1) 
    for i in range(num_samples):
        ls = ["-","--","-."][i%3]
        plt.plot(wavelengths[i], sum_intensities[i], label=titles[i], linestyle=ls)
    plt.title("Sum of Intensities")
    plt.legend()

    #Plots average intensity against wavelengths
    plt.subplot(2, 1, 2)
    for i in range(num_samples):
        ls = ["-","--","-."][i%3]
        plt.plot(wavelengths[i], avg_intensities[i], label=titles[i], linestyle=ls)
    plt.title("Average Intensity")
    plt.legend()

    #Saves the plots and data to the PDF
    export_pdf.savefig(dpi=1200)
    plt.close()

    #Writes some output statistics for the plots
    col_titles = ["Sample", "Peak Value\n (Sum)", "Peak Value\n (Average)", "Average Intensity across\n all Wavelengths", 
                    "Standard deviation\n (Sum)", "Standard deviation\n (Average)"]
    row_titles = titles
    stats = [row_titles, [], [], [], [], []]
    for i in range(num_samples):
        stats[1].append(max(sum_intensities[i]))
        stats[2].append(max(avg_intensities[i]))
        stats[3].append(sum(avg_intensities[i])/float(len(avg_intensities[i])))
        stats[4].append(np.std(sum_intensities[i]))
        stats[5].append(np.std(avg_intensities[i]))

    df = pd.DataFrame({col_titles[0]:stats[0], col_titles[1]:stats[1], col_titles[2]:stats[2], col_titles[3]:stats[3], col_titles[4]:stats[4], col_titles[5]:stats[5]}, columns=col_titles)

    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis("off")
    ax.axis("tight")
    l = 0.05
    b = 0.05
    w = 0.9
    h = 0.9

    table = ax.table(cellText=df.values, colLabels=col_titles, loc="center", cellLoc="center", bbox=[l, b, w, h])
    table.auto_set_font_size(False)
    table.set_fontsize(5)
    table.auto_set_column_width(col=col_titles)
    export_pdf.savefig()

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

def generate_output(sample_file, background_file, document_title):
    '''
    Main function for generating the PDF
    '''
    print("\nGenerating Output...")
    with PdfPages(r"/home/frederique/PhysicsLabPythonCodes/Optrode/Records/{}.pdf".format(document_title)) as export_pdf:

        #Title page
        plt.figure()
        plt.axis("off")
        plt.text(0.5, 0.5, "Linear Flow Fluorescence", fontsize=28, ha="center", va="center")
        export_pdf.savefig()
        plt.close()        
        
        #The background file should be a Powermeter reading
        background_data = read_data(background_file[0])
        intensities_b, times_b, wavelengths_b = [background_data[4]], [background_data[5]], [background_data[7]]
        power_readings, pm_times = [background_data[2]], [background_data[3]]

        intensities, times, wavelengths, net_intensities, photodiodes, pd_times = [], [], [], [], [], []
        for i in range(num_files.get()):
            sample_data = read_data(sample_file[i]) 
            adjust = 2
            if len(sample_file[i]) == 8: #If the powermeter readings are included, adjusts to read the correct data
                adjust = 0
            #List of lists with the data in order: Photodiode readings and times (0, 1), Powermeter readings and times (2, 3), Spectrometer intensities, times and wavelengths (4, 5, 6)
            #If there is no Powermeter readings the data will be in the order: Photodiode readings and times (0, 1), Spectrometer intensities, times and wavelengths (2, 3, 4)
            #It may be necessary to adjust these indices, depending on how the Optrode saves the data
            intensities.append(sample_data[4-adjust])
            times.append(sample_data[5-adjust])
            wavelengths.append(sample_data[7-adjust])
            net_intensity = np.subtract(sample_data[4-adjust], intensities_b[0])
            net_intensities.append(net_intensity)

            photodiodes.append(sample_data[0])
            pd_times.append(sample_data[1])

        #Completes plots
        photodiode_plots(photodiodes, pd_times, export_pdf, file_titles)
        powermeter_plots(power_readings, pm_times, export_pdf, background_title)
        mainplots(intensities, wavelengths, export_pdf, "Sample Readings", False, file_titles)
        mainplots(intensities_b, wavelengths_b, export_pdf, "Backgroud Reading", True, background_title)
        mainplots(net_intensities, wavelengths, export_pdf, "Net Readings", False, file_titles)
        #Integrals
        fig = plt.figure(figsize=(12, 12))
        fig.text(0.5, 0.04, "Time (s)", ha="center", fontsize=14)
        fig.text(0.04, 0.5, "Intensity (a.u.)", ha="center", rotation="vertical", fontsize=14)
        plt.suptitle("Integrals", fontsize=24)
        plt.subplot(2, 1, 1)
        integral_plot(intensities, wavelengths, export_pdf, times, "Sample")
        plt.legend()
        plt.subplot(2, 1, 2)
        integral_plot(net_intensities, wavelengths, export_pdf, times, "Net")
        plt.legend()            
        export_pdf.savefig()
        plt.close()

    print("\n{}.pdf Generated".format(document_title))

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

def browse(n):
    '''
    Allows user to browse directories for input data files
    '''
    filenames = filedialog.askopenfilenames(initialdir="/home/frederique/PhysicsLabPythonCodes/Optrode/Records/", filetypes=(("HDF5 files", "*.hdf5"),))
    if n == 1: #sample files
        sample_files.set(filenames)
        global file_list, file_titles
        file_list = sample_files.get().split(", ")
        file_titles = []
        trim(file_list, file_titles, num_files.get())
    else: #background
        background_file.set(filenames)
        global background_list, background_title
        background_list = background_file.get().split(", ")
        background_title = []
        trim(background_list, background_title, 1)

def close_program():
    root.destroy()
    exit()

if __name__ == '__main__':

    root = Tk()
    root.title("Linear Flow Fluorescence Analysis")
    root.geometry("750x350")
    root.grid_columnconfigure(0, weight=1)

    document_title = StringVar(value="Optrode Analysis") #Default value
    sample_files = StringVar(value="")
    background_file = StringVar(value="")
    paradigm_mode = StringVar(value="c") # Paradigm mode ('c' or 'm')
    integral_lower_bound = IntVar(value=505) #Bounds for integral plot
    integral_upper_bound = IntVar(value=550)

    #The following details all of the instructions for the Tkinter GUI

    frame1 = Frame(root)
    frame1.grid(row=0, column=0)

    title = Label(frame1, text="Linear Flow Fluorescence Analysis", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=2, column=0, padx=8, pady=6, sticky=W+E)

    global num_files
    num_files = IntVar(value=1)
    nums = [1, 2, 3, 4, 5, 6]
    drop = OptionMenu(frame2, num_files, *nums)
    drop.grid(row=1, column=2, padx=12, pady=6, sticky=W)

    ttl = Label(frame2, text="Enter document title: ", font=(None, 11))
    ttl.grid(row=0, column=1, padx=28, pady=6, sticky=W)

    nmlbl = Label(frame2, text="Select number of Samples: ", font=(None, 11))
    nmlbl.grid(row=1, column=1, padx=28, pady=6, stick=W)

    ptlbl = Label(frame2, text="Select Paradigm type: ", font=(None, 11))
    ptlbl.grid(row=2, column=1, padx=28, pady=6, sticky=W)

    modes = [("Continuous", "c"), ("Multi Integration", "m"),]
    i = 2
    for text, mode in modes:
        box = Radiobutton(frame2, text=text, variable=paradigm_mode, value=mode)
        box.grid(row=i, column=2, padx=4, pady=4, sticky=W)
        i = i+1

    splbl = Label(frame2, text="Select Sample Measurement(s): ", font=(None, 11))
    splbl.grid(row=4, column=1, padx=28, pady=6, sticky=W)

    bglbl = Label(frame2, text="Select Background Measurement: ", font=(None, 11))
    bglbl.grid(row=5, column=1, padx=28, pady=6, sticky=W)

    lblbl = Label(frame2, text="Integral Lower Bound: ", font=(None, 11))
    lblbl.grid(row=6, column=1, padx=28, pady=6, sticky=W)

    ublbl = Label(frame2, text="Integral Upper Bound: ", font=(None, 11))
    ublbl.grid(row=7, column=1, padx=28, pady=6, sticky=W)

    input_box0 = Entry(frame2, textvariable=document_title, width=40)
    input_box0.grid(row=0, column=2, padx=6)

    tag = Label(frame2, text=".pdf", font=(None, 10))
    tag.grid(row=0, column=3, sticky=W)

    input_box1 = Entry(frame2, textvariable=sample_files, width=40)
    input_box1.grid(row=4, column=2, padx=12)

    input_box2 = Entry(frame2, textvariable=background_file, width=40)
    input_box2.grid(row=5, column=2, padx=12)

    input_box3 = Entry(frame2, textvariable=integral_lower_bound, width=15)
    input_box3.grid(row=6, column=2, padx=12, sticky=W)

    input_box4 = Entry(frame2, textvariable=integral_upper_bound, width=15)
    input_box4.grid(row=7, column=2, padx=12, sticky=W)

    but1 = Button(frame2, text="Browse", command=lambda: browse(1))
    but1.grid(row=4, column=3)
    but2 = Button(frame2, text="Browse", command=lambda: browse(2))
    but2.grid(row=5, column=3)

    frame3 = Frame(root)
    frame3.grid(row=5, column=0)

    but3 = Button(frame3, text="Start", command=lambda: generate_output(file_list, background_list, document_title.get()))
    but3.grid(row=5, column=0, padx=10, pady=10)

    but4 = Button(frame3, text="Close", command=lambda: close_program())
    but4.grid(row=5, column=1)

    root.mainloop()
