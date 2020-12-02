# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:00:10 2016
This code runs the Optrode version 2
@author: Yaqub Jonmohamadi
"""
#%%
import socket
import subprocess
import struct
import tempfile

import os
import os.path
import sys
import glob
import datetime
import time
import bisect
from multiprocessing import Process, Value, Array
from Tkinter import * #Tk, Text, BOTH, W, N, E, S, RAISED, Frame, Message, LEFT, TOP, BOTTOM, DISABLED, NORMAL, PhotoImage, StringVar, Toplevel
from ttk import Button, Style, Label, Entry, Notebook, Scale
from tkFileDialog import askopenfilename
from PIL import Image, ImageTk

import h5py
import DAQT7_Objective as DAQ
import SeaBreeze_Objective as SBO
import ThorlabsPM100_Objective as P100
import numpy as np
import matplotlib.pyplot as plt
#%%
time_start =  time.time()

#Naming the DAQ ports
Blue_Laser = "FIO1"
Blue_Shutter = "DAC1"
Green_Laser = "FIO0"
Green_Shutter = "DAC0"

#Initialising values
Open_Shutter = 5
Close_Shutter = 0

PhotoDiode_Port = "AIN1"
#Spectrometer_Trigger_Port = "DAC1"

Paradigm = 'm' # Paradigm refers to the continuous ('c') or multi-integration ('m'), for the performance of the Optrode
               # Default for the GUI is set to 'm'

def DAQ_Read_Process(num_DAC_Samples):
    '''
    A function for reading the DAQ analogue inpute on AINX
    Reading the Photodiode
    '''
    for i in range(num_DAC_Samples):
        DAQ_Signal[DAQ_Index[0]], DAQ_Time[DAQ_Index[0]] = DAQ1.readPort(PhotoDiode_Port)
        DAQ_Index[0] = i
    DAQ_Is_Read.value = 1

def DAQ_Speed_Test(num_DAQ_Tests):
    '''
    A function for testing the speed of the DAQ analogue inpute on AINX
    '''
    Scratch_Signal = 0
    Start_Time = time.time()
    for i in range(num_DAQ_Tests):
        Scratch_Signal, Scratch_Time = DAQ1.readPort(PhotoDiode_Port)
    Duration = time.time() - Start_Time
    Mean_Time = Duration/float(num_DAQ_Tests)
    print ('DAQ Analogue reading requires %s seconds per sample \n' %Mean_Time)
    return Mean_Time

def Power_Speed_Test(num_Power_Tests):
    '''
    A function for testing the speed of the Power meter
    '''
    Scratch_Signal = 0
    Start_Time = time.time()
    for i in range(num_Power_Tests):
        Scratch_Signal, Scratch_Time = Power_meter.readPower()
    Duration = time.time() - Start_Time
    Mean_Time = Duration/float(num_Power_Tests)
    print ('Power meter reading requires %s seconds per sample \n' %Mean_Time)
    return Mean_Time

def Power_Read_Process(num_power_samples):
    '''
    A function for reading the Power meter -- readPower method returns the Power meter value and the time at which the value was taken.
    Hence, the Power_Signal Array stores the values and Power_Time Array stores the times.
    '''
    for i in range(num_power_samples):
        Power_Signal[Power_Index[0]], Power_Time[Power_Index[0]] = Power_meter.readPower()
        Power_Index[0] = i
    Power_Is_Read.value = 1 #Has been read

def Timer_Multi_Process(Time_In_Seconds):
    '''
    Interrupt like delays (s)
    Usage Ex: Px = Process(target=Timer_Multi_Process, args=(Timer_length,))
    Px.start() and in your code constantly check for "Timer_Is_Done"
    '''
    if Timer_Is_Over.value is 1:
        print ('ERROR: The previous timer is still running. This timer can be only called once at a time.')
    time.sleep(Time_In_Seconds)
    Timer_Is_Over.value = 1

def Spectrometer_Init_Process(Integration_Time, Trigger_mode):
    '''
    A function for initializing the spectrometer (integration time and triggering mode)
    Integration_Time is given in milliseconds
    '''
    print ("Spectrometer is initialised")
    Spec1.setTriggerMode(Trigger_mode)
    time.sleep(0.01)
    Spec1.setIntegrationTime(Integration_Time*1000) #Converted to Microseconds
    time.sleep(0.01)
    Spectrometer_Init_Done.value = 1 #Completed

def Spectrometer_Speed_Test(num_Spec_Tests):
    '''
    A function for testing the spectrometer speed in free running mode
    '''
    Test_Integration_Time = Spec1.Handle.minimum_integration_time_micros/float(1000)
    Spectrometer_Init_Process(Test_Integration_Time, 0)
    Mean_Time = 0
    #I = 0
    Threshold = 1.2
    while True:
        Start_Time = time.time()
        Spectrometer_Read_Process(num_Spec_Tests)
        Spectrometer_Index[0] = 0
        '''
        while I < num_Spec_Tests:
            Scratch_Signal, Mean_Time = Spec1.readIntensity(True, True)
            I = I + 1
        '''
        Duration = (time.time() - Start_Time)*1000
        Mean_Time = Duration/float(num_Spec_Tests)
        if (Mean_Time < Threshold*float(Test_Integration_Time)):
            print ('Finished. Duration: %s' %Duration)
            print ('Mean time %s' %Mean_Time)
            print ('Test_Integration_Time %s' %Test_Integration_Time)
            break
        else: #The integration time was not enough for the spectromer so it is increased by 0.5 portion of the minimum integration time
            print ('Mean time %s' %Mean_Time)
            print ('Test_Integration_Time %s' %Test_Integration_Time)
            print ('Duration: %s' %Duration)
            Test_Integration_Time = Test_Integration_Time + Spec1.Handle.minimum_integration_time_micros/float(2000)
            Spectrometer_Init_Process(Test_Integration_Time, 0)
            #I = 0

    print ('Spectrometer minimum integration time is %output_file ms and the practical minimum integration time is %output_file ms \n' %(Spec1.Handle.minimum_integration_time_micros/float(1000), Test_Integration_Time))
    return Test_Integration_Time

def Spectrometer_Read_Process(num_Spectrometer_Samples):
    '''
    A function for reading the spectrometer intensities
    '''
    for i in range(num_Spectrometer_Samples):
        Intensity, Time = Spec1.readIntensity(True, True) #Returns the intensity value and time at which it was taken
        Current_Spec_Record = Spec1.readWavelength
        Intensity_Array = np.asanyarray(Intensity)
        Intensity_Matrix = np.ndarray.reshape(Intensity_Array, (len(Intensity_Array), 1))
        Full_Spec_Records[:, Spectrometer_Index[0]] = Intensity_Matrix[Min_Wave_Index:Max_Wave_Index]
        #Spec_Time[Spectrometer_Index[0]] = Time
        Spectrometer_Index[0] = i
        Spectrometer_is_read.value = 1
        #print ("spectrometer Index is %i" %Spectrometer_Index[0])

def Check_Output_Folder():
    '''
    A function that looks for the 'Records' output folder and creates it if it does not exist
    '''
    global Path_to_Records
    Path_to_Records = os.path.abspath(os.path.join(os.getcwd())) + "/Records"
    if not os.path.exists(Path_to_Records):
        os.makedirs(Path_to_Records)

def Multi_Integration_Paradigm(Integration_list_MilSec, Integration_Buffer_Time, Shutter_Delay, num_power_samples):
    '''
    Below paradigm is based on free running of the spectrometer.
    '''
    if (Power_meter.Error == 0): #No Error
        Pros_Power = Process(target=Power_Read_Process, args=(num_power_samples,))
        Pros_Power.start()

    Integration_Base = Integration_list_MilSec[-1] + 2*Integration_Buffer_Time

    Trigger_mode = 0    # Free running
    Process_Order = 2   # Order of Process: 0 = setting the timer for the duration of the laser exposure
               # 1 = turning off the laser (the exposure time is over) and waiting for the spectromer to finish the current integration cycle
               # 2 = integration time is over and spectrometer must be read and setting the timer for the buffer time of the integration cycle

    Spectrometer_Init_Done.value = 0 #Not Complete
    #Pros_Spectrometer_Init = Process(target = Spectrometer_Init_Process, args=(Integration_list_MilSec[Spectrometer_Index[0]], Trigger_mode))
    #Pros_Spectrometer_Init = Process(target = Spectrometer_Init_Process, args=(Integration_Base, Trigger_mode))
    #Pros_Spectrometer_Init.start()
    time.sleep(0.1)

    Spectrometer_is_read.value = 0 #Has not been read
    #Pros_Spec = Process(target=Spectrometer_Read_Process, args=(len(Integration_list_MilSec),))
    #Pros_Spec.start() #Begins reading

    print ('Step1, First integration does not have laser exposure', time.time())
    while (Spectrometer_Index[0] < len(Integration_list_MilSec)):
        if  (Timer_Is_Over.value == 1) and (Process_Order == 0): #Sets up the timer for the next integration cycle
            Timer_Is_Over.value = 0 #Resets timer
            P_Timer = Process(target=Timer_Multi_Process, args=(Integration_list_MilSec[Spectrometer_Index[0]-1]/float(1000),))
            P_Timer.start()
            DAQ1.writePort(Chosen_Shutter, Open_Shutter)
            Ref_Time[DAQ_Index[0]] = time.time()
            print ('Step2, Raising edge has started',  time.time())
            Process_Order = 1
        elif(Timer_Is_Over.value == 1) and (Process_Order == 1): #Closes the shutter
            DAQ1.writePort(Chosen_Shutter, Close_Shutter)
            print ('Step3, Falling edge has started',  time.time())
            Ref_Time[DAQ_Index[0]] = time.time()
            Timer_Is_Over.value = 0
            Process_Order = 2
        elif(Spectrometer_is_read.value == 1) and (Process_Order == 2): #Once an integration cycle is complete
            while Spectrometer_Index[0] < len(Integration_list_Milsec): #Repeats until all integration cycles are complete
                print ('Step4, Spectrometer has been read',  time.time())
                Spectrometer_is_read.value = 0
                #Full_Spec_Records[:, Spectrometer_Index[0]] = Current_Spec_Record #Updates the matrix with readings
                Spectrometer_Index[0] += 1

            Timer_Is_Over.value = 0
            P_Timer = Process(target=Timer_Multi_Process, args=(Integration_Buffer_Time/float(1000),))
            P_Timer.start()
            Pros_Spectrometer_Init = Process(target = Spectrometer_Init_Process, args=(Integration_Base, Trigger_mode))
            Pros_Spectrometer_Init.start()
            Process_Order = 0 #Reset Process
        u = (DAQ1.readPort(PhotoDiode_Port))
        DAQ_Signal[DAQ_Index[0]], DAQ_Time[DAQ_Index[0]] = DAQ1.readPort(PhotoDiode_Port)
        Ref_Time[DAQ_Index[0]] = time.time()
    Pros_Spec.terminate()

    Timer_Is_Over.value = 0
    P_Timer = Process(target=Timer_Multi_Process, args=(Integration_Buffer_Time/float(1000),))
    P_Timer.start()
    print ('time of Timer start %s output_file:' %time.time())
    while  Timer_Is_Over.value == 0: #Not over
        DAQ_Signal[DAQ_Index[0]], DAQ_Time[DAQ_Index[0]] = DAQ1.readPort(PhotoDiode_Port)
        #print (DAQ_Signal[DAQ_Index[0]])
        DAQ_Index[0] = DAQ_Index[0] + 1
        Ref_Time[DAQ_Index[0]] = time.time()
    P_Timer.terminate()
    if (Power_meter.Error == 0):
        Pros_Power.terminate()

def Continuous_Paradigm(Integration_Continuous, num_Spectrometer_Samples, num_DAC_Samples, num_power_samples, num_BakGro_Spec):
    '''
    A function used to read data for the Continuous Paradigm
    '''
    #First take the samples of the readings of the Power Meter
    if (Power_meter.Error == 0): #No Error
        Pros_Power = Process(target=Power_Read_Process, args=(num_power_samples,))
        Pros_Power.start()

    #Initialises Spectrometer with Integration time and Trigger mode
    Spectrometer_Init_Done.value = 0 #Not Completeime
    Pros_Spectrometer_Init = Process(target = Spectrometer_Init_Process, args=(Integration_Continuous, 0))
    Pros_Spectrometer_Init.start()

    Spectrometer_is_read.value = 0 #Not Read
    Pros_Spec = Process(target=Spectrometer_Read_Process, args=(num_Spectrometer_Samples,))
    Timer_Is_Over.value = 0 #Not Over
    P_Timer = Process(target=Timer_Multi_Process, args=(0.1,))
    P_Timer.start()
    while  Timer_Is_Over.value == 0:
        DAQ_Signal[DAQ_Index[0]], DAQ_Time[DAQ_Index[0]] = DAQ1.readPort(PhotoDiode_Port)
        DAQ_Index[0] += 1

    Spectrometer_Read_Process(num_Spectrometer_Samples)
    DAQ1.writePort(Chosen_Shutter, Open_Shutter)

    Pros_DAQ = Process(target=DAQ_Read_Process, args=(num_DAC_Samples,))
    Pros_DAQ.start()
    Pros_DAQ.join()
    Last_DAQ_Signals = 0 # This is a flag bit used to detect the last part of the continuous recording
    while (int(Spectrometer_Index[0]) < num_Spectrometer_Samples):
        if int(Spectrometer_Index[0]) == (num_Spectrometer_Samples-num_BakGro_Spec):
            DAQ1.writePort(Chosen_Shutter, Close_Shutter) #Closes shutter for final few readings
        DAQ_Signal[DAQ_Index[0]], DAQ_Time[DAQ_Index[0]] = DAQ1.readPort(PhotoDiode_Port)
        DAQ_Index[0] += 1
        Spectrometer_Index[0] += 1
        if  Spectrometer_is_read.value == 1:
            Spectrometer_is_read.value = 0
            #Full_Spec_Records[:, Spectrometer_Index[0]] = Current_Spec_Record #Update columns with spectrum

    if (Power_meter.Error == 0):
        Pros_Power.terminate()

def Begin_Test():
    '''
    Test if parameters are appropriate, update UI and start the test.
    '''
    MIN_INTEGRATION_TIME, MAX_INTEGRATION_TIME = Spec1.Handle.minimum_integration_time_micros/1000.0, 2000
    MIN_RECORD_TIME, MAX_RECORD_TIME = 1, 3600
    MIN_WAVELENGTH, MAX_WAVELENGTH = 100, 10000

    #Device detection errors
    if Spec1.Error == 1:
        debug("ERROR: Session failed, could not detect spectrometer.")
    elif DAQ1.Error == 1:
        debug("ERROR: Session failed, could not detect DAQ, try unplugging and plugging back in.")
    #Potential errors when the user navigates the GUI
    elif not is_number(integration_time.get()):
        debug("ERROR: Integration duration is not a number.")
    elif float(integration_time.get()) < MIN_INTEGRATION_TIME:
        debug("ERROR: Integration duration is smaller than " + str(MIN_INTEGRATION_TIME) + ".")
    elif float(integration_time.get()) > MAX_INTEGRATION_TIME:
        debug("ERROR: Integration duration is greater than " + str(MAX_INTEGRATION_TIME) + ".")
    elif not is_number(record_time.get()):
        debug("ERROR: Recording duration is not a number.")
    elif float(integration_time.get()) < MIN_RECORD_TIME:
        debug("ERROR: Recording duration is smaller than " + str(MIN_RECORD_TIME) + ".")
    elif float(integration_time.get()) > MAX_RECORD_TIME:
        debug("ERROR: Recording duration is greater than " + str(MAX_RECORD_TIME) + ".")
    elif not is_number(min_length.get()):
        debug("ERROR: Minimum wavelength is not a number.")
    elif float(min_length.get()) < MIN_WAVELENGTH:
        debug("ERROR: Minimum wavelength is smaller than " + str(MIN_WAVELENGTH) + ".")
    elif float(min_length.get()) > MAX_WAVELENGTH:
        debug("ERROR: Minimum wavelength is greater than " + str(MAX_WAVELENGTH) + ".")
    elif not is_number(max_length.get()):
        debug("ERROR: Maximum wavelength is not a number.")
    elif float(min_length.get()) < MIN_WAVELENGTH:
        debug("ERROR: Maximum wavelength is smaller than " + str(MIN_WAVELENGTH) + ".")
    elif float(min_length.get()) > MAX_WAVELENGTH:
        debug("ERROR: Maximum wavelength is greater than " + str(MAX_WAVELENGTH) + ".")
    elif float(min_length.get()) >= float(max_length.get()):
        debug("ERROR: Minimum wavelength is smaller than maximum wavelength.")
    elif paradigm_mode.get() != "c" and paradigm_mode.get() != "m":
        debug("ERROR: Invalid paradigm mode selected.")
    elif shutter_mode.get() != Green_Shutter and shutter_mode.get() != Blue_Shutter:
        debug("ERROR: Invalid shutter selected.")
    else:
        # If no errors, then start the test after 10ms, to allow UI time to update
        debug("SETTING UP...")
        but1.config(state=DISABLED)
        Disable_UI(root)
        root.after(10, Perform_Test)

def Perform_Test():
    '''
    Main function that performs entire spectrometer test.
    '''
    #Initialises values as integers equal to 0
    global Spectrometer_is_read, Spectrometer_Init_Done, DAQ_Is_Read, Power_Is_Read, Timer_Is_Over
    Spectrometer_is_read = Spectrometer_Init_Done = DAQ_Is_Read = Power_Is_Read = Timer_Is_Over = Value('i', 0)

    Integration_Time = 100                                        # Integration time in ms
    Spec1.setTriggerMode(3)                                       # It is set for free running mode
    #Spec1.setIntegrationTime(Integration_Time*1000)              # Integration time is in microseconds when using the library

    # Check to see if the output folder named Records exists, creates this folder if it does not already exist
    Check_Output_Folder()

    #Ensures shutters are closed before starting
    DAQ1.writePort(Green_Shutter, Close_Shutter)
    DAQ1.writePort(Blue_Shutter, Close_Shutter)

    # Initialising the variables
    Integration_list_MilSec = [8, 16, 32, 64, 128, 256, 512, 1024] #Integration time for the spectrometer in ms
    Shutter_Delay = 4  #ms

    num_DAQ_Tests = 20000
    DAQ_SamplingRate = DAQ_Speed_Test(num_DAQ_Tests)*1000  #Shows the sampling speed in ms

    #Takes data from devices to determine the range of appropriate wavelengths
    global Wavelengths, Min_Wave_Index, Max_Wave_Index, Spec_Time, Spectrometer_Index, Current_Spec_Record
    Wavelengths = Spec1.Handle.wavelengths()
    Min_Wave_Index = max(bisect.bisect(Wavelengths, float(min_length.get())-1), 0)
    Max_Wave_Index = bisect.bisect(Wavelengths, float(max_length.get()))
    Wavelengths = Wavelengths[Min_Wave_Index:Max_Wave_Index]
    Current_Spec_Record = np.zeros(shape=len(Wavelengths), dtype=float)
    num_Spec_Tests = 500
    Spec_Time = np.array(np.zeros(shape=(num_Spec_Tests, 1), dtype = float ))
    Spectrometer_Index = np.array(np.zeros(shape=(1, 1), dtype = int))
    
    #Commented out for speed Spec_SamplingRate = Spectrometer_Speed_Test(num_Spec_Tests)
    Integration_Buffer_Time = 100       #ms               # This is for the spectrometer. This is the time from the integration started till shutter opens
    #DurationOfReading = np.sum(Integration_list_MilSec)  + len(Integration_list_MilSec)*Delay_Between_Integrations   # Duration of reading in seconds.
    DurationOfReading = (Integration_list_MilSec[-1] + Integration_Buffer_Time + Shutter_Delay*3)*len(Integration_list_MilSec)     # Duration of reading in seconds.
    num_BakGro_Spec = 10       # This is for continuous reading and refers to the last few spectrom reading wich are background and the laser is off

    # Reads and stores parameter values the user inputs into the GUI
    global Chosen_Shutter
    Chosen_Shutter = shutter_mode.get() #Green Shutter or Blue Shutter
    Integration_Continuous = float(integration_time.get()) #Integration time in ms
    DurationOfReading = float(record_time.get())*1000.0 + num_BakGro_Spec*float(integration_time.get()) #Recording duration in s
    Filename_Prefix = filename.get() #Filename, default is OptrodeData
    if Filename_Prefix == "":
        Filename_Prefix = "OptrodeData"

    if (Power_meter.Error == 0):
        #Powermeter_SamplingRate = 5.1     #ms
        num_Power_Tests = 200
        Power_SamplingRate = Power_Speed_Test(num_Power_Tests)*1000            #Shows the sampling speed in ms

    #Defining the size of the arrays and matrices for recording the signals beased on the duration of the recording #
    #num_DAC_Samples = int(round((DurationOfReading + DurationOfReading/4) /DAQ_SamplingRate))        # Number of samples for DAQ analogue to digital converter (AINx).
    num_DAC_Samples = int((2 * DurationOfReading)/DAQ_SamplingRate)

    if (Power_meter.Error == 0):
        num_power_samples = int((1.5 * DurationOfReading) /Power_SamplingRate)
    else:
        num_power_samples = 0
    #No_Power_Sample = int(round(DurationOfReading/Powermeter_SamplingRate))
    # Number of samples for P100D Power meter to read.
    # Roughly P100 can read the power every 2.7 ms.

    if (paradigm_mode.get() == 'c'): #Continuous paradigm
        num_Spectrometer_Samples =  int(round(float(DurationOfReading)/float(float(Integration_Continuous))))  # Number of samples for spectrometer to read.
    else: #Multi Integration
        num_Spectrometer_Samples =  len(Integration_list_MilSec) # Number of samples for spectrometer to read.

    rerun = "First"
    while True:
        #Initialising Variables
        global Full_Spec_Records
        Full_Spec_Records = np.array(np.zeros(shape=(len(Wavelengths), num_Spectrometer_Samples), dtype = float)) #Values in 2D Array
        #Each column is a spectrum with each row being the intensity value for each wavelength
        Current_Spec_Record = np.zeros(shape=(len(Wavelengths), 1), dtype = float)
        global DAQ_Signal, DAQ_Time, DAQ_Index, DAQ_Index_Total, Ref_Signal, Ref_Time
        DAQ_Signal = DAQ_Time = np.array(np.zeros(shape=(num_DAC_Samples, 1), dtype = float))
        DAQ_Index = DAQ_Index_Total = np.array(np.zeros(shape=(1, 1), dtype = int))
        Ref_Signal = Ref_Time = np.array(np.zeros(shape=(num_DAC_Samples, 1), dtype = float))

        if Power_meter.Error == 0:
            global Power_Signal, Power_Time, Power_Index
            Power_Signal = np.array(np.zeros(shape=(num_power_samples, 1), dtype=float)) #Keeps track of values read from Power Meter
            Power_Time = np.array(np.zeros(shape=(num_power_samples, 1), dtype=float)) #Time at which values were taken
            Power_Index = np.array(np.zeros(shape=(1, 1), dtype=int)) #Counter

        # Wait for user to press Start
        but2.config(state=NORMAL)
        debug("Ready to start paradigm. Press start to begin.")
        but2.wait_variable(wait_var)
        but2.config(state=DISABLED)

        # Starting the chosen paradigm -- The main process
        if (paradigm_mode.get() == 'm'):
            Multi_Integration_Paradigm(Integration_list_MilSec, Integration_Buffer_Time, Shutter_Delay, num_power_samples)
        else:
            Continuous_Paradigm(float(Integration_Continuous), num_Spectrometer_Samples, num_DAC_Samples, num_power_samples, num_BakGro_Spec)

        #Loading the Spectrometer Array to a matrix before saving and plotting
        num_wavelengths = len(Wavelengths)
        # Closing the devices
        Spec_Details = Spec1.readDetails()
        DAQ_Details = DAQ1.getDetails()
        DAQ1.writePort(Chosen_Shutter, Close_Shutter)

        #The file containing the records (HDF5 format)
        os.chdir(Path_to_Records)
        if is_suff.get() == 1: #If a suffix is desired, a timestamp is added to the name of the output data
            File_name_Suffix = str('%s' %datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H-%M-%S'))+ ".hdf5"
            File_name = Filename_Prefix + '-' + File_name_Suffix
        else:
            File_name = Filename_Prefix

        #Opens a file to output data
        output_file = h5py.File(File_name, "w")

        # Saving the recorded signals in HDF5 format
        Optrode_DAQ = output_file.create_group('DAQT7')
        output_file.create_dataset('DAQT7/PhotoDiode', data = np.asanyarray(DAQ_Signal[:]))
        output_file.create_dataset('DAQT7/TimeIndex', data = np.asanyarray(DAQ_Time[:]))
        Optrode_DAQ.attrs['DAQT7 Details'] = np.string_(DAQ_Details)

        Spec_Records_1D = np.ndarray.reshape(Full_Spec_Records, len(Wavelengths)*num_Spectrometer_Samples)
        Optrode_Spectrometer = output_file.create_group('Spectrometer')
        output_file.create_dataset('Spectrometer/Intensities', data = np.asanyarray(Spec_Records_1D))
        output_file.create_dataset('Spectrometer/Time_Index', data = np.asanyarray(Spec_Time))
        output_file.create_dataset('Spectrometer/WaveLength', data = np.asanyarray(Wavelengths))
        Optrode_Spectrometer.attrs['Spectrometer Details'] = np.string_(Spec_Details)

        if (Power_meter.Error == 0):
            Optrode_Power = output_file.create_group('PM100_PowerMeter')
            output_file.create_dataset('PM100_PowerMeter/Power', data = np.asanyarray(Power_Signal[:]))
            output_file.create_dataset('PM100_PowerMeter/TimeIndex', data = np.asanyarray(Power_Time[:]))
            Optrode_DAQ.attrs['PowerMeter Details'] = np.string_(DAQ_Details)
        output_file.close()

        # Plotting the spectrometer and the photodiode recordings
        if plots[0].get() == 1:
            plt.figure()
            plt.plot(np.asarray(DAQ_Time[0:DAQ_Index[0]])-DAQ_Time[0], np.asanyarray(DAQ_Signal[0:DAQ_Index[0]]))
            plt.title('Photo diode')
            plt.xlabel('Elapsed time (s)')
            plt.ylabel('Voltage (v)')

        plt.figure()
        for i in range(num_Spectrometer_Samples):
            plt.plot(np.asarray(Spec1.readWavelength()[Min_Wave_Index:Max_Wave_Index]), Full_Spec_Records[:, i]) 
        #Plot the wavelengths vs each column/spectrum -- should be all on the same plot but separate spectrums denoted by colour
        plt.title('Spectrometer recordings')
        plt.xlabel('Wavelength (nano meter)')
        plt.ylabel('Intensity')
        plt.show()
        #plt.plot(np.asarray(Ref_Time[0:DAQ_Index[0]]) - DAQ_Time[0], Ref_Signal[0:DAQ_Index[0]])

        # Estimate the latencies of the devices
        if plots[1].get() == 1:
            plt.figure()
            plt.subplot(1,2,1)
            DAQ_Latency = np.asanyarray(DAQ_Time[0:DAQ_Index[0]])
            DAQ_Latency[0] = 0
            for i in range(1, DAQ_Index[0]):
                DAQ_Latency[i] = DAQ_Time[i] - DAQ_Time[i-1]
            plt.subplot(1,3,1)
            plt.plot(DAQ_Latency)
            plt.ylabel("Time (s)")
            plt.title("DAQ latencies")
            plt.show()

            plt.subplot(1,2,2)
            Spec_Latency = np.asarray(Spec_Time[0:np.int(Spectrometer_Index[0])])
            Spec_Latency[0] = 0
            for i in range(1,Spectrometer_Index[0]):
                Spec_Latency[i] = np.float(Spec_Time[i] - Spec_Time[i-1])
            plt.plot(Spec_Latency[1:])

            plt.ylabel("Time (s)")
            plt.title("Spectrometer integration durations")
            plt.show()

        #PLots the readings for the powermeter
        if plots[2].get() == 1 and Power_meter.Error == 0:
            plt.figure()
            plt.subplot(1,2,1)
            Power_Latency = np.asanyarray(Power_Time[0:Power_Index[0]])
            Power_Latency[0] = 0
            for I in range(1,int(Power_Index[0])):
                Power_Latency[I] = Power_Time[I] - Power_Time[I-1]
            plt.subplot(1,3,1)
            plt.plot(Power_Latency)
            plt.ylabel("Time (s)")
            plt.title("Power latencies")
            plt.show()
            plt.subplot(1,2,2)
            plt.plot(np.asarray(Power_Time[0:Power_Index[0]]) - Power_Time[0], np.asanyarray(Power_Signal[0:Power_Index[0]]))
            plt.title('Power Meter')
            plt.xlabel('Elapsed time (s)')
            plt.ylabel('Power (w)')
            plt.show()

        print('\n')
        print('Data is saved to the Records folder. \n')

        # Wait for user to decide what to do next -- Re-run, change or quit
        but3.config(state=NORMAL)
        but4.config(state=NORMAL)
        debug("Re-run test with same parameters, change parameters, or quit.")
        but3.wait_variable(wait_var)
        but3.config(state=DISABLED)
        but4.config(state=DISABLED)

        # If we clicked the 'change' but, quit loop, otherwise keep going.
        if wait_var.get() == 2:
            but1.config(state=NORMAL)
            Disable_UI(root, False)
            break

def Close_GUI():
    '''
    Closes GUI.
    '''
    time.sleep(0.1)
    DAQ1.close()
    Spec1.close()
    root.destroy()

def Disable_UI(parent, disable=True):
    '''
    Disables sections of UI that are not buttons. Functions recursively.
    '''
    for w in parent.winfo_children():
        if w.winfo_class() == "TEntry" or w.winfo_class() == "Radiobutton":
            if disable == True:
                w.config(state=DISABLED)
            else:
                w.config(state=NORMAL)
        else:
            Disable_UI(w, disable)

def debug(msg):
    '''
    Used for printing Error messages
    '''
    print(msg)
    error_msg.set(msg)

def is_number(s):
    '''
    Checks if a string is either a float or an integer
    '''
    try:
        float(s)
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    '''
    Main function that opens the UI and allows the user to setup/start the processes
    '''
    # Checks all devices are connected and operational
    print("")
    Spec1 = SBO.DetectSpectrometer()
    print("")
    DAQ1 = DAQ.DetectDAQT7()
    print("")
    Power_meter = P100.DetectPM100D()

    # Creating GUI
    root = Tk()
    root.geometry("380x450+150+150")
    root.grid_columnconfigure(0, weight=1)

    # Setting variables -- These are the default values that can be changed
    filename = StringVar(value="")      # Filename of output data
    is_suff = IntVar(value=0)           # 0 or 1 depending on whether or not to add a suffix to filename
    integration_time = StringVar(value="15")    # Integration time of spectrometer (ms)
    record_time = StringVar(value="10")    # Recording duration of spectrometer (s)
    min_length = StringVar(value="500")    # Minimum wavelength to record (nm)
    max_length = StringVar(value="600")    # Maximum wavelength to record (nm)
    paradigm_mode = StringVar(value="c")     # Paradigm mode ('c' or 'm')
    shutter_mode = StringVar(value=Blue_Shutter)   # Shutter to use (Blue_Shutter or Green_Shutter)

    error_msg = StringVar(value=" ")    # Variable that stores the error message displayed on GUI
    wait_var = IntVar(value=0)          # Variable used to keep GUI waiting for user input during tests

    small_entry = 6
    large_entry = 19

    '''
    The following code sets all of the parameters for the GUI
    '''

    # Title frame
    frame1 = Frame(root)
    frame1.grid(row=0, column=0)

    title = Label(frame1, text="Linear Flow Fluorescense", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    # Parameter frame
    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=1, column=0, padx=8, pady=4, sticky=W+E)

    row_filename, row_suffix, row_recording_duration, row_integration_time, row_wavelength_range = 1, 2, 3, 4, 5

    fnlbl = Label(frame2, text="Select Filename", font=(None, 11))
    fnlbl.grid(row=row_filename, column=1, padx=6, pady=(6,0), sticky=E)
    fntxt = Entry(frame2, textvariable=filename, width=large_entry)
    fntxt.grid(row=row_filename, column=2, padx=6, pady=(6,0), columnspan=3)

    sulbl = Label(frame2, text="Add Suffix", font=(None, 11))
    sulbl.grid(row=row_suffix, column=1, padx=9, pady=6, sticky=E)
    fntxt = Checkbutton(frame2, variable=is_suff)
    fntxt.grid(row=row_suffix, column=2, padx=0, pady=6, sticky=W)

    rdlbl = Label(frame2, text="Recording duration (s)", font=(None, 11))
    rdlbl.grid(row=row_recording_duration, column=1, padx=6, pady=6, sticky=E)
    rdtxt = Entry(frame2, textvariable=record_time, width=large_entry)
    rdtxt.grid(row=row_recording_duration, column=2, padx=6, pady=6, columnspan=3)

    itlbl = Label(frame2, text="Integration time (ms)", font=(None, 11))
    itlbl.grid(row=row_integration_time, column=1, padx=6, pady=6, sticky=E)
    ittxt = Entry(frame2, textvariable=integration_time, width=large_entry)
    ittxt.grid(row=row_integration_time, column=2, padx=6, pady=6, columnspan=3)

    wrlbl = Label(frame2, text="Wavelength range (nm)", font=(None, 11))
    wrlbl.grid(row=row_wavelength_range, column=1, padx=6, pady=6, sticky=E)
    wrtxt = Entry(frame2, textvariable=min_length, width=small_entry)
    wrtxt.grid(row=row_wavelength_range, column=2, padx=6, pady=6, sticky=W)
    wrlbl = Label(frame2, text="to", font=(None, 11))
    wrlbl.grid(row=row_wavelength_range, column=3, padx=0, pady=6)
    wrtxt = Entry(frame2, textvariable=max_length, width=small_entry)
    wrtxt.grid(row=row_wavelength_range, column=4, padx=6, pady=6, sticky=E)

    # Frame for checkbox frames
    framech = Frame(root, relief=RAISED, borderwidth=1)
    framech.grid(row=2, column=0, padx=8, pady=4, sticky=W+E)
    framech.grid_columnconfigure(1, weight=1)
    framech.grid_columnconfigure(2, weight=1)
    framech.grid_columnconfigure(3, weight=1)

    # Paradigm select frame
    frame3 = Frame(framech)
    frame3.grid(row=1, column=1, padx=4, pady=4)
    frame3.grid_columnconfigure(1, weight=1)

    palbl = Label(frame3, text="Paradigm:")
    palbl.grid(row=0, column=0, padx=6, pady=6, sticky=W)

    modes = [("Continuous", "c"), ("Multi Integration", "m"),]
    i = 1
    for text, mode in modes:
        box = Radiobutton(frame3, text=text, variable=paradigm_mode, value=mode)
        box.grid(row=i, column=0, padx=4, pady=4, sticky=W)
        i = i+1

    # Shutter select frame
    frame4 = Frame(framech)
    frame4.grid(row=1, column=2, padx=4, pady=4)
    frame4.grid_columnconfigure(1, weight=1)

    shlbl = Label(frame4, text="Shutter:")
    shlbl.grid(row=0, column=0, padx=6, pady=6, sticky=W)

    modes = [("Blue", Blue_Shutter), ("Green", Green_Shutter),]
    i = 1
    for text, mode in modes:
        box = Radiobutton(frame4, text=text, variable=shutter_mode, value=mode)
        box.grid(row=i, column=0, padx=4, pady=4, sticky=W)
        i = i+1

    # Plot select frame
    frame7 = Frame(framech)
    frame7.grid(row=1, column=3, padx=4, pady=4)
    frame7.grid_columnconfigure(1, weight=1)

    num_plots = 3
    plots = []
    for i in range(num_plots):
        plots.append(IntVar())
        check = Checkbutton(frame7, text="plot"+str(i+1), variable=plots[i])
        check.grid(row=i, column=0, padx=4, pady=4, sticky=W)

    # Error message
    erlbl = Message(framech, width=360, textvariable=error_msg, font=(None, 11))
    erlbl.grid(row=2, column=1, columnspan=3)

    #Button frames
    frame5 = Frame(root)
    frame5.grid(row=3, column=0)
    frame6 = Frame(root)
    frame6.grid(row=4, column=0)

    but1 = Button(frame5, text="Setup", command=Begin_Test)
    but1.grid(row=1, column=1, padx=10, pady=10)
    but2 = Button(frame5, text="Start", command=lambda: wait_var.set(0), state=DISABLED)
    but2.grid(row=1, column=2, padx=10, pady=10)
    but3 = Button(frame6, text="Re-run", command=lambda: wait_var.set(1), state=DISABLED)
    but3.grid(row=1, column=1, padx=10, pady=10)
    but4 = Button(frame6, text="Change", command=lambda: wait_var.set(2), state=DISABLED)
    but4.grid(row=1, column=2, padx=10, pady=10)
    but5 = Button(frame6, text="Close", command=Close_GUI)
    but5.grid(row=1, column=3, padx=10, pady=10)

root.mainloop()
