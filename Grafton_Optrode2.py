# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:57:46 2016
@author: fred
"""

import h5py
import DAQT7_Obj as DAQ
import SeaBreeze_Obj as SB
import ThorlabsPM100_Objective as P100
import numpy as np
import matplotlib.pyplot as plt

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
import timeit
from timeit import default_timer
from multiprocessing import Process, Value, Array, Pool
from Tkinter import * #Tk, Text, BOTH, W, N, E, S, RAISED, Frame, Message, LEFT, TOP, BOTTOM, DISABLED, NORMAL, PhotoImage, StringVar, Toplevel
#from ttk import Button, Style, Label, Entry, Notebook, Scale
from tkFileDialog import askopenfilename
from PIL import Image, ImageTk

from backports.time_perf_counter import perf_counter
import timeit
from timeit import default_timer

time_start =  time.time()

# ######################### Naming the DAQ ports ##########################
# FIO0 = shutter of the green laser and FIO1 is the shutter of the blue laser
# FIO2 = is the green laser and the FIO3 is the blue laser
Green_Laser = "FIO1"
Green_Shutter = "FIO3"
Blue_Laser = "FIO0"
Blue_Shutter = "FIO2"
Green_Shutter_CloseDelay = 0.03  #Delay in seconds for the shutter to close
Blue_Shutter_CloseDelay = 0.010  #Delay in seconds for the shutter to close

PhotoDiod_Port = "AIN0"
Spectrometer_Trigger_Port = "DAC0"

# ####################### Interrupt like delays (s) ####################### '''
# Usage Ex: Px = Process(target=Timer_Multi_Process, args=(Timer_time,))
# Px.start() and in your code constantly check for "Timer_Is_Done"
# Generally start the process and use a while loop, with the condition statement checking if the timer_is_done value is 1

def Timer_Multi_Process(Time_In_Seconds):
    if Timer_Is_Done.value is 1:
        print 'Error: This timer can be run one at a time. Either the previous timer is still running, or Timer_Is_Done bit is reset from previous timer run'
    wait = Time_In_Seconds
    timer_finish = timeit.default_timer() + wait
    while timeit.default_timer() < timer_finish:
        pass
    Timer_Is_Done.value = 1

def Timer_Multi_Process2(Time_In_Seconds):
    if Timer_Is_Done2.value is 1:
        print 'Error: This timer can be run one at a time. Either the previous timer is still running, or Timer_Is_Done bit is reset from previous timer run'
    wait = Time_In_Seconds
    timer_finish = timeit.default_timer() + wait
    while timeit.default_timer() < timer_finish:
        pass
    Timer_Is_Done2.value = 1

# # A function for initializing the spectrometer (integration time and triggering mode '''
def SB_Init_Process(Spec_handle,Integration_time, Trigger_mode):
    print 'Spectrometer is initialized'
    SB.Init(Spec_handle,Integration_time, Trigger_mode)

# ########## A function for reading the spectrometer intensities ########### '''
def SB_Read_Process(Spec_handle):
    #print 'Spectrumeter is waiting'
    #Correct_dark_counts = True
    #Correct_nonlinearity = True
    Intensities, Time = SB.Read(Spec_handle, True, True)
    #print Intensities
    SB_Current_Record[:] = Intensities[min_wave_index:max_wave_index]
    Spec_Time[Spec_Sampl_Index.value] = Time
    #print(Spec_Sampl_Index.value)
    SB_Is_Done.value = 1
    #print "Intensities are read"
    return

# ######## A function for reading the DAQ analogue inpute on AINX ########
def DAQ_Read():
    results = DAQ.AIN_Read(DAQ_handle, PhotoDiod_Port)
    read_signal[DAC_Sampl_Index.value] = results[0]
    return results[0], time.time()

# ######## A function for reading the Powermeter readings ########
def Power_Read_Process(No_Power_Sample):
    '''
    A function for reading the Power meter -- readPower method returns the Power meter value and the time at which the value was taken.
    Hence, the power_signal Array stores the values and power_time Array stores the times.
    '''
    for i in range(power_index.value, power_index.value + No_Power_Sample):
        power_signal[i], power_time[i] = Power_meter.readPower()
        power_index.value = power_index.value + 1

def Continuous_Paradigm():
    #Spec_Integration_Time = 20000                       # Integration time for free running mode
    P1 = Process(target=SB_Init_Process, args=(Spec_handle, continuous_integration_time*1000, 0))
    P1.start()
    P1.join()
    time.sleep(0.1)

    #Starts reading the powermeter
    if (Power_meter.Error == 0):
        Pros_Power = Process(target=Power_Read_Process, args=(No_Power_Sample, ))
        Pros_Power.start()

    # ## The main loop for recording the spectrometer and the photodiode ##
    #Takes a few photodiode readings before opening the shutter
    P_Timer = Process(target=Timer_Multi_Process, args=(0.1,)) # keep the laser on before opening the shutter
    P_Timer.start()
    while Timer_Is_Done.value == 0:
        DAC_Sampl_Index.value += 1
        read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
    Timer_Is_Done.value = 0
    DAQ.Digital_Ports_Write(DAQ_handle, Shutter_Port, 1)

    # The first integration time starts here, when laser is off
    P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
    P2.start()

    SB_Is_Done.value = 0
    while SB_Is_Done.value == 0:
            DAC_Sampl_Index.value += 1
            read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()

    start = time.time()
    DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)  #Turn the laser on right after first intergration is over, and start recording    

    while Spec_Sampl_Index.value < No_Spec_Sample:
        Start_time = time.time()
        SB_Is_Done.value = 0
        P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
        P2.start()
        while SB_Is_Done.value == 0:
            DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
            read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()

        SB_Full_Records[:,Spec_Sampl_Index.value] = SB_Current_Record[:]
        #print(Spec_Sampl_Index.value)
        #print SB_Full_Records[0,Spec_Sampl_Index]
        #print(Spec_Sampl_Index.value)
        Spec_Sampl_Index.value = Spec_Sampl_Index.value + 1

    DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 0)
    DAQ.Digital_Ports_Write(DAQ_handle, Shutter_Port, 0)
    end = time.time()
    print("Duration: {}".format(end-start))

    P_Timer = Process(target=Timer_Multi_Process, args=(0.1,))
    P_Timer.start()
    while Timer_Is_Done.value == 0:
        DAC_Sampl_Index.value += 1
        read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()

    DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)

    if P2.is_alive():
        P2.terminate()

    if (Power_meter.Error == 0):
        Pros_Power.terminate()

def Multi_Integration_Paradigm():

    if (Power_meter.Error == 0):
        Pros_Power = Process(target=Power_Read_Process, args=(No_Power_Sample, ))
        Pros_Power.start()

    DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)       #Laser is on
    P1 = Process(target=SB_Init_Process, args=(Spec_handle,Integration_base,3))
    P1.start()
    time.sleep(0.1)

    State = 0

    # ## The main loop for recording the spectrometer and the photodiod ##
    while Integration_index.value < len(Integration_list_sec):
        #Start_time = time.time()
        Timer_Is_Done.value = 0
        DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)
        time.sleep(0.02)
        P_Timer = Process(target=Timer_Multi_Process, args=(Integration_marging,)) # keep the laser on before opening the shutter
        P_Timer.start()

        P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
        P2.start()
        time.sleep(0.02)
        DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)       #Laser is on
        DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 5)  # Spec is edge-triggered  and start ~20ms later to acquire
        time.sleep(0.01)
        DAQ.Digital_Ports_Write(DAQ_handle, Shutter_Port, 0)
        time.sleep(0.01)
        Open_delay = time.time()

        while Timer_Is_Done.value == 0:
            DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
            read_signal_ref[DAC_Sampl_Index.value] = State
            read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
        #print 'Elapsed time %f' %(time.time() - Start_time)
        Latch_Laser_Detect = 0

        #if SB_Is_Done.value == 1:
        #    print 'Eroooooooor'
        while SB_Is_Done.value == 1:
            SB_Is_Done.value = 0
            print 'Spectrometer Error, will retrigger the Spectrometer'
            DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)
            time.sleep(0.04)
            P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
            P2.start()
            DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 5)  # Spec is edge-triggered  and start ~20ms later to acquire
            time.sleep(0.04)

        Timer_Is_Done.value = 0
        P_Timer = Process(target=Timer_Multi_Process, args=(Integration_base/float(1000000) - Integration_marging*2 - Integration_list_sec[Integration_index.value],)) # keep the laser on before opening the shutter
        P_Timer.start()
        while Timer_Is_Done.value == 0:
            DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
            read_signal_ref[DAC_Sampl_Index.value] = State
            read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
        DAQ.Digital_Ports_Write(DAQ_handle, Shutter_Port, 1)       #Shutter opens in ~9ms since now
        #time.sleep(0.02)
        CurrentDelay = Integration_list_sec[Integration_index.value] - 0.005
        while SB_Is_Done.value == 0:
            DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
            read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
            read_signal_ref[DAC_Sampl_Index.value] = State
            if  (Latch_Laser_Detect == 0) & (read_signal[DAC_Sampl_Index.value] > 0.3):
                State = 4.5
                print "Integration Time: " + str(CurrentDelay + 0.005)
                Timer_Is_Done.value = 0
                P_Timer = Process(target=Timer_Multi_Process, args=(CurrentDelay,)) # keep the laser on before opening the shutter
                P_Timer.start()
                Latch_Laser_Detect = 1

                while Timer_Is_Done.value == 0:
                    #print 'step 3'
                    DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
                    read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
                    read_signal_ref[DAC_Sampl_Index.value] = State
                State = 0

                DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 0)       #Laser is off

                DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
                read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
                read_signal_ref[DAC_Sampl_Index.value] = State

                DAQ.Digital_Ports_Write(DAQ_handle, Shutter_Port, 0)       #Shutter closes in ~9ms since now

                P_Timer2 = Process(target=Timer_Multi_Process2, args=(Shutter_CloseDelay,)) # keep the laser on before opening the shutter
                P_Timer2.start()
                while Timer_Is_Done2.value == 0:
                    #print 'step 3'
                    DAC_Sampl_Index.value = DAC_Sampl_Index.value + 1
                    read_signal[DAC_Sampl_Index.value], read_time[DAC_Sampl_Index.value] = DAQ_Read()
                    read_signal_ref[DAC_Sampl_Index.value] = State

                Timer_Is_Done2.value = 0
                #time.sleep(0.030)
                DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)       #Laser is on
                #Open_delay[Shutter_Open_Delay_Index] = time.time() - Open_time
                #print ('%f' %time.time(), 'step9'); print (SB_Is_Done.value)

            #read_time_ref[DAC_Sampl_Index] = time.time()
            #read_signal_ref[DAC_Sampl_Index] = Command_State
        SB_Full_Records[:,Spec_Sampl_Index.value] = SB_Current_Record[:]
        Spec_Sampl_Index.value += 1

        Integration_index.value += 1
        print "The whole integration cycle: %f" % (time.time() - Open_delay)
        SB_Is_Done.value = 0
        Timer_Is_Done.value = 0

    DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)
    time.sleep(0.01)
    P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
    P2.start()
    time.sleep(0.01)
    DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 5)  # Spec is edge-triggered  and start ~20ms later to acquire
    time.sleep(0.05)
    while SB_Is_Done.value == 1:
        SB_Is_Done.value = 0
        print 'Spectrometer Error, will retrigger the Spectrometer'
        DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)
        time.sleep(0.01)
        P2 = Process(target=SB_Read_Process, args=(Spec_handle,))
        P2.start()
        DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 5)  # Spec is edge-triggered  and start ~20ms later to acquire
        time.sleep(0.05)
    while SB_Is_Done.value == 0:
        time.sleep(0.1)
        #print 'BackGround'
    SB_Full_Records[:,Spec_Sampl_Index.value] = SB_Current_Record[:]
    DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)

    DAQ.Digital_Ports_Write(DAQ_handle, Laser_Port, 1)       #Laser stays on

    if (Power_meter.Error == 0):
        Pros_Power.terminate()

def Begin_Test():
    '''
    Test if parameters are appropriate, update UI and start the test.
    '''
    MIN_INTEGRATION_time, MAX_INTEGRATION_time = Spec_handle.minimum_integration_time_micros/1000.0, 2000
    MIN_RECORD_time, MAX_RECORD_time = 1, 3600
    MIN_WAVELENGTH, MAX_WAVELENGTH = 100, 10000

    #Device detection errors
    #if Spec_handle.Error == 1:
        #debug("ERROR: Session failed, could not detect spectrometer.")
    #elif DAQ1.Error == 1:
        #debug("ERROR: Session failed, could not detect DAQ, try unplugging and plugging back in.")
    #Potential errors when the user navigates the GUI
    if not is_number(integration_time.get()):
        debug("ERROR: Integration duration is not a number.")
    elif not is_number(record_time.get()):
        debug("ERROR: Recording duration is not a number.")
    elif float(integration_time.get()) < MIN_RECORD_time:
        debug("ERROR: Recording duration is smaller than " + str(MIN_RECORD_time) + ".")
    elif float(integration_time.get()) > MAX_RECORD_time:
        debug("ERROR: Recording duration is greater than " + str(MAX_RECORD_time) + ".")
    elif not is_number(min_length.get()):
        debug("ERROR: Minimum wavelength is not a number.")
    elif float(min_length.get()) < MIN_WAVELENGTH:
        debug("ERROR: Minimum wavelength is smaller than " + str(MIN_WAVELENGTH) + ".")
    elif float(min_length.get()) > MAX_WAVELENGTH:
        debug("ERROR: Minimum wavelength is greater than " + str(MAX_WAVELENGTH) + ".")
    elif not is_number(max_length.get()):
        debug("ERROR: Maximum wavelength is not a number.")
    elif float(max_length.get()) < MIN_WAVELENGTH:
        debug("ERROR: Maximum wavelength is smaller than " + str(MIN_WAVELENGTH) + ".")
    elif float(max_length.get()) > MAX_WAVELENGTH:
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
        root.after(5, Perform_Test)

def Perform_Test():
    # ################# Detecting the spectrometer and the DAQ ###########
    global Laser_Port, Shutter_Port, Shutter_CloseDelay
    if shutter_mode.get() == Green_Shutter:
        Laser_Port = Green_Laser
        Shutter_Port = Green_Shutter
        Shutter_CloseDelay = Green_Shutter_CloseDelay
    else:
        Laser_Port = Blue_Laser
        Shutter_Port = Blue_Shutter
        Shutter_CloseDelay = Blue_Shutter_CloseDelay       

    # ##################### Initializing the variables ###################
    global Wavelengths, min_wave_index, max_wave_index, Integration_base, No_DAC_Sample, No_Spec_Sample, Integration_list_sec, No_Power_Sample
    Wavelengths = Spec_handle.wavelengths()
    min_wave_index = max(bisect.bisect(Wavelengths, float(min_length.get())-1), 0)
    max_wave_index = bisect.bisect(Wavelengths, float(max_length.get()))
    Wavelengths = Wavelengths[min_wave_index:max_wave_index]

    #Integration_list = [8000, 16000, 32000, 64000, 128000, 256000, 512000, 1024000, 2048000]
    #Integration_list_sec = [0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048]
    Integration_list_sec = [] #Integration time for the spectrometer in ms
    for i in range(int(multi_integration_times.get())):
        t = 2**(i+3) #Adds integration times 8, 16, 32, 64, 128... to the list of integration times
        Integration_list_sec.append(t/1000.0)
    #No_Spec_Sample = 500                             #(In seconds) This is the duration before the external edge trigger is given to the spectrometer while the specrumeter started the integration period
    #Integration_base = Integration_list_sec[-1]*1000000 + Integration_marging*2000000    # This is the integration time applied for all the trials
    No_DAC_Sample = 10000 # Number of samples for Photodiod per iteration of the laser exposer. Every sample takes ~0.6 ms.

    if paradigm_mode.get() == "C" or paradigm_mode.get() == "c":
        duration = float(record_time.get())
    else:
        duration = len(Integration_list_sec) * 2.1
    No_Power_Sample = int(duration * 1.2 * 250)
    #print(No_Power_Sample)

    while True:
        global No_Spec_Sample, continuous_integration_time
        durationOfReading = float(record_time.get())*1000
        continuous_integration_time = float(integration_time.get())
        No_Spec_Sample =  int(round(durationOfReading/continuous_integration_time))  

        global SB_Is_Done, SB_Current_Record, Timer_Is_Done, Timer_Is_Done2, SB_Full_Records, Spec_Time
        SB_Is_Done = Value('i', 0)
        SB_Current_Record = Array('d', np.zeros(shape=(len(Wavelengths) ,1), dtype = float ))
        SB_Is_Done.value = 0
        Timer_Is_Done = Value('i', 0)
        Timer_Is_Done.value = 0
        Timer_Is_Done2 = Value('i', 0)
        Timer_Is_Done2.value = 0
        SB_Full_Records = np.zeros(shape=(len(Wavelengths),  No_Spec_Sample), dtype = float )
        Spec_Time = Array("d", np.zeros(shape=(No_Spec_Sample, 1), dtype=float))
     
        global DAC_Sampl_Index, Integration_index, Spec_Sampl_Index, power_index
        DAC_Sampl_Index = Value("i", -1)
        Integration_index = Value("i", 0)
        Spec_Sampl_Index = Value("i", 0)
        power_index = Value("i", 0)

        # Wait for user to press Start
        but2.config(state=NORMAL)
        debug("Ready to start paradigm. Press start to begin.")
        but2.wait_variable(wait_var)
        but2.config(state=DISABLED)

        global read_signal, read_time, Integration_base, Integration_marging, power_signal, power_time
        if paradigm_mode.get() == "C" or paradigm_mode.get() == "c":
            Integration_marging = 0.3          
            Integration_base = 32000
            read_signal = np.zeros(No_DAC_Sample*Integration_base)
            read_time   = np.zeros(No_DAC_Sample*Integration_base)
            power_signal = Array("d", np.zeros(shape=((No_Power_Sample), 1), dtype=float))
            power_time = Array("d", np.zeros(shape=((No_Power_Sample), 1), dtype=float))
            Continuous_Paradigm()

        else: #Mode is M or m
            global read_signal_ref, read_time_ref
            Integration_marging = 0.2        
            Integration_base = 2*1000000
            read_signal_ref = np.zeros(No_DAC_Sample*len(Integration_list_sec))
            read_time_ref   = np.zeros(No_DAC_Sample*len(Integration_list_sec))
            read_signal = np.zeros(No_DAC_Sample*len(Integration_list_sec))
            read_time   = np.zeros(No_DAC_Sample*len(Integration_list_sec))

            power_signal = Array("d", np.zeros(shape=((No_Power_Sample), 1), dtype=float))
            power_time = Array("d", np.zeros(shape=((No_Power_Sample), 1), dtype=float))
            Multi_Integration_Paradigm()

        # ########### Saving the recorded signals in HDF5 format ############
        #Path_to_Fred_Codes = os.path.abspath(os.path.join( os.getcwd(), os.pardir)) + "/Fred"
        #os.chdir(Current_Path)

        Current_Path = os.path.abspath(os.path.join( os.getcwd()))
        Path_to_Records = os.path.abspath(os.path.join( os.getcwd(), os.pardir)) + "/Records"
        os.chdir(Path_to_Records)


        #Filename = "water_4_" + str('%i' %time.time())+ ".hdf5"
        #Filename = "power_check_" + str('%s' %datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H-%M-%S'))+ ".hdf5"
        #Filename = "Opterode_Recording_At" + str('%i' %time.time())+ ".hdf5"

        if filename.get() == "":
            filename_middle = "OptrodeData"
        else:
            filename_middle = filename.get()

        if shutter_mode.get() == Blue_Shutter:
            filename_prefix = "Grafton Blue "
        else:
            filename_prefix = "Grafton Green "
        if paradigm_mode.get() == "c":
            filename_prefix += "Cont "
        else:
            filename_prefix += "Multi "
        if is_suff.get() == 1: #If a suffix is desired, a timestamp is added to the name of the output data
            filename_suffix = str('%s' %datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H-%M-%S'))
            Filename = filename_prefix + "- " + filename_middle + " - " + filename_suffix + ".hdf5"
        else:
            Filename = filename_prefix + "- " + filename_middle + ".hdf5"

        f = h5py.File(Filename, "w")
        Spec_sub1 = f.create_group("Spectrumeter")
        Spec_specification = Spec_sub1.create_dataset("Spectrumeter", (10,), dtype='f')
        Spec_specification.attrs['Serial Number'] = np.string_(Spec_handle.serial_number)
        Spec_specification.attrs['Model'] = np.string_(Spec_handle.model)
        Spec_wavelength = f.create_dataset('Spectrumeter/Wavelengths', data = Wavelengths)

        read_signal2 = np.zeros(DAC_Sampl_Index.value)
        read_time2   = np.zeros(DAC_Sampl_Index.value)
        read_signal2[:] = read_signal[0:DAC_Sampl_Index.value]
        read_time2 = np.asarray(read_time[0:DAC_Sampl_Index.value])
        read_time2 -= read_time2[0]
        Spec_Time2 = np.asarray(Spec_Time[:])
        Spec_Time2 -= Spec_Time2[0]
        Spec_intensities = f.create_dataset('DAQT7/DAC_Readings', data = read_signal2)
        Spec_intensities = f.create_dataset('DAQT7/DAC_Time_Stamps', data = read_time2)
        Spec_intensities = f.create_dataset('Spectrometer/Intensities', data = SB_Full_Records[:, :Spec_Sampl_Index.value])
        Spec_intensities = f.create_dataset('Spectrometer/Spectrometer_Time_Stamps', data = Spec_Time2[:Spec_Sampl_Index.value])
        if paradigm_mode.get() == "m" or paradigm_mode.get() == "M":
            Spec_Time2 = np.asarray(Spec_Time[:len(Integration_list_sec)+1])
            read_signal_ref2 = np.zeros(DAC_Sampl_Index.value)
            read_signal_ref2[:] = read_signal_ref[0:DAC_Sampl_Index.value]
            #Spec_intensities = f.create_dataset('DAQT7/DAC_Command_Signal', data = read_signal_ref2)

        if (Power_meter.Error == 0):
            power_time2 = np.asarray(power_time[:power_index.value])
            power_time2 -= power_time2[0]
            Spec_intensities = f.create_dataset("PM100_PowerMeter/Power_Readings", data=power_signal[:power_index.value])
            Spec_intensities = f.create_dataset("PM100_PowerMeter/Power_Time_Stamps", data=power_time2)
        f.close()
        #print(power_index.value)

        print 'File %s is saved in %s' %(Filename ,  Path_to_Records)

        #Path_to_Fred_Codes = os.path.abspath(os.path.join( os.getcwd(), os.pardir)) + "/Fred"
        os.chdir(Current_Path)

        #SB.Close(Spec_handle)
        #DAQ.Close(DAQ_handle)
        # ######### Plotting the spectrumeter and the photodiod recordings ########
        plt.figure()
        #plt.plot(Spec_handle.wavelengths()[1:1044],SB_Full_Records[1:1044,:])
        for i in range(No_Spec_Sample):
            plt.plot(Wavelengths, SB_Full_Records[:, i])
        plt.title('Spectrometer recordings')
        plt.xlabel('Wavelength (nano meter)')
        plt.ylabel('Intensity')
        plt.show()

        if plots[0].get() == 1:
            plt.figure()
            read_time_index = read_time2 - read_time2[0]
            plt.plot(read_time_index,read_signal2, label = "Photo Diode")
            if paradigm_mode.get() == "c" or paradigm_mode.get() == "C":
                pass
            else:
                read_time_index_ref2 = read_time_ref - read_time_ref[0]
                plt.plot(read_time_index,read_signal_ref2, label = "Command Signal")
            plt.legend()
            plt.title('Photo diode')
            plt.xlabel('Time (s)')
            plt.ylabel('Voltage (v)')
            plt.show()

        if plots[1].get() == 1:
            plt.figure()
            plt.subplot(2, 1, 1)
            read_delay = np.zeros(len(read_time_index))
            for i in range(len(read_time_index)-1):
                read_delay[i] = read_time_index[i+1]-read_time_index[i]
            plt.plot(read_delay)
            plt.title("DAQ Latency")
            plt.xlabel("Iterations index")
            plt.ylabel("Time (s)")

            plt.subplot(2, 1, 2)
            Spec_Delay = np.zeros(len(Spec_Time2))
            for i in range(len(Spec_Time2)-1):
                #Spec_Delay[i] = (SB_Full_Records[0,i+1]-SB_Full_Records[0,i])
                Spec_Delay[i] = Spec_Time2[i+1]-Spec_Time2[i]
            plt.plot(Spec_Delay[:len(Spec_Time2)])
            plt.title('Spectrometer Latency')
            plt.xlabel('Iterations index')
            plt.ylabel('Time (s)')
            plt.tight_layout()
            plt.show()

        if plots[2].get() == 1 and Power_meter.Error == 0:
            plt.figure()
            plt.subplot(2, 1, 1)
            plt.plot(power_time2, power_signal[:power_index.value])  
            plt.title("Powermeter Readings")
            plt.ylabel("Power (W)")
            plt.xlabel("Time (s)")

            plt.subplot(2, 1, 2)
            power_latency = np.zeros(len(power_time2))
            for i in range(len(power_time2)-1):
                power_latency[i] = power_time2[i+1] - power_time2[i]
            plt.plot(power_latency[:len(power_time2)])
            plt.title("Powermeter Latency")
            plt.xlabel("Iterations index")
            plt.ylabel("Time (s)")
            plt.tight_layout()
            plt.show()          
            #pm

        # Wait for user to decide what to do next -- Re-run, change or quit
        but3.config(state=NORMAL)
        but4.config(state=NORMAL)
        debug("Re-run test with same parameters, change parameters, or quit.")
        but3.wait_variable(wait_var)
        but3.config(state=DISABLED)
        but4.config(state=DISABLED)
        # If we clicked the 'change' button, quit loop, otherwise keep going.
        if wait_var.get() == 2:
            but1.config(state=NORMAL)
            Disable_UI(root, False)
            break

def Close_GUI():
    '''
    Closes GUI.
    '''
    time.sleep(0.1)
    DAQ.Close(DAQ_handle)
    SB.Close(Spec_handle)
    root.destroy()
    exit()

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
    Spec_handle = SB.Detect()
    DAQ_handle = DAQ.Init()
    # ############## All the ports are off at the beginning ##############
    DAQ.DAC_Write(DAQ_handle, Spectrometer_Trigger_Port, 0)
    DAQ.Digital_Ports_Write(DAQ_handle, Blue_Laser, 1)       #Laser is on
    DAQ.Digital_Ports_Write(DAQ_handle, Green_Laser, 1)       #Laser is on
    DAQ.Digital_Ports_Write(DAQ_handle, Green_Shutter, 0)       #Shutter is close
    DAQ.Digital_Ports_Write(DAQ_handle, Blue_Shutter, 0)       #Shutter is close
    DAQ.DAC_Write(DAQ_handle, 'DAC1', 0)
    print("")
    Power_meter = P100.DetectPM100D()

    # Creating GUI
    root = Tk()
    root.geometry("380x475+150+150")
    root.grid_columnconfigure(0, weight=1)

    # Setting variables -- These are the default values that can be changed
    # You may want to change these to the default values for your measurements (e.g record_time, shutter_mode, integration_time, wavelengths etc)
    filename = StringVar(value="")      # Filename of output data
    is_suff = IntVar(value=0)           # 0 or 1 depending on whether or not to add a suffix to filename
    integration_time = StringVar(value="20")    # Integration time of spectrometer (ms)
    multi_integration_times = StringVar(value="8")
    record_time = StringVar(value="10")    # Recording duration of spectrometer (s)
    min_length = StringVar(value="500")    # Minimum wavelength to record (nm)
    max_length = StringVar(value="750")    # Maximum wavelength to record (nm)
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

    title = Label(frame1, text="Linear Flow Fluorescence", font=(None, 13))
    title.grid(row=0, column=1, padx=10, pady=10)

    # Parameter frame
    frame2 = Frame(root, relief=RAISED, borderwidth=1)
    frame2.grid(row=1, column=0, padx=8, pady=4, sticky=W+E)

    row_filename, row_suffix, row_recording_duration, row_integration_time, row_multi_times, row_wavelength_range = 1, 2, 3, 4, 5, 6

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

    milbl = Label(frame2, text="Num Multi Integrations", font=(None, 11))
    milbl.grid(row=row_multi_times, column=1, padx=6, pady=6, sticky=E)
    milbl = Entry(frame2, textvariable=multi_integration_times, width=large_entry)
    milbl.grid(row=row_multi_times, column=2, padx=6, pady=6, columnspan=3)

    wrlbl = Label(frame2, text="Wavelength range (nm)", font=(None, 11))
    wrlbl.grid(row=row_wavelength_range, column=1, padx=6, pady=6, sticky=E)
    wrtxt = Entry(frame2, textvariable=min_length, width=small_entry)
    wrtxt.grid(row=row_wavelength_range, column=2, padx=6, pady=6, sticky=W)
    wrlbl = Label(frame2, text="to", font=(None, 11))
    wrlbl.grid(row=row_wavelength_range, column=3, padx=6, pady=6)
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

    #Plot select frame
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
    but5 = Button(frame6, text="Close", command=lambda: Close_GUI())
    but5.grid(row=1, column=3, padx=10, pady=10)

root.mainloop()

