##############################################################
#
#   ECE 588 Final Project PSU Winter 2024
#
#   Python script to automate graphing of the speedup data 
#   from benchmark tests.  Parses speedup data from textfile
#   and graphs against number of processors
#
#
#   Alexander Maso (Main code) 
#   Mohamed (added second graph for time vs number of processors) : 03/08/2024
#
##############################################################

import numpy as np
import matplotlib.pyplot as plt
import sys

def read_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)  # Assuming the first row is header
    return data[:, 0], data[:, 1], data[:, 2]  # Parse number of threads, time, and speedup

def plot_data(threads, speedup, title='Processor Speedup'):
    plt.figure(figsize=(10, 6))
    plt.plot(threads, speedup, marker='o', linestyle='-', color='b')
    plt.title(title)
    plt.xlabel('Number of Processors')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.show()

def plot_time(threads, time, title='Processor Time'):
    plt.figure(figsize=(10, 6))
    plt.plot(threads, time, marker='s', linestyle='-', color='r')
    plt.title(title)
    plt.xlabel('Number of Processors')
    plt.ylabel('Time (Seconds)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_data_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    threads, time, speedup = read_data(file_path)
    plot_data(threads, speedup, title=f'Processor Speedup for {file_path}')
    plot_time(threads, time, title=f'Processor Time for {file_path}')
