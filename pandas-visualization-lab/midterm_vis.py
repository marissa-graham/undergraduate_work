#!//anaconda/bin/python
# hist.py
"""Midterm
Shane McQuarrie
Math 510
28 October 2015

Instructions:
    Run this file from the terminal to see all answers at once, or import the
    file in IPython and follow the directions to see one answer at a time.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def visualize():
    """View the data provided in 'midterm_data.csv'."""

    # Extract the Data
    data = pd.read_csv("midterm_data.csv")
    data.hist(grid=True)
    plt.suptitle("Summary of the Midterm data")
    plt.show()

    data.plot(kind='hist')
    plt.show()

if __name__ == '__main__':
    visualize()