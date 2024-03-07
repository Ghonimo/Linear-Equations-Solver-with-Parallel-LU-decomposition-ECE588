##############################################################
#
#   ECE 588 Final Project PSU Winter 2024
#
#   Python script to automate graphing of the speedup data 
#   from benchmark tests.  Parses speedup data from textfile
#   and graphs against number of processors
#
#
#   Alexander Maso 
#
##############################################################


import numpy as np
import matplotlib.pyplot as plt
import sys

def read_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)  # Assuming the first row is header
    return data[:, 0], data[:, 2]  # Parse number of threads and speedup

def plot_data(threads, speedup, title='Processor Speedup'):
    plt.figure(figsize=(10, 6))
    plt.plot(threads, speedup, marker='o', linestyle='-', color='b')
    plt.title(title)
    plt.xlabel('Number of Processors')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_data_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    threads, speedup = read_data(file_path)
    plot_data(threads, speedup, title=f'Processor Speedup for {file_path}')
