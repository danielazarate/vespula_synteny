#!/bin/python

# A script for calculating slope across variable windows.
# 
# functools.cmp_to_key : transform an old-style comparison function to a key function.
# key function : a callable that returns a value used for sorting or ordering.
# itertools.groupby() : Makes an iterator that returns consecutive keys and groups from the iterable, 
# you can use a key function with this that computes a key value for each element. 

import pandas as pd
import numpy as np
import argparse 
from functools import cmp_to_key
from itertools import groupby


# Set up a helpful argparser ecosystem:
parser = argparse.ArgumentParser(description='Calculates slope across windows in a genetic map.')
parser.add_argument('--input', help='input file')
parser.add_argument('--output', help='filename for the output file')
args = parser.parse_args()

# Function breakpoint: 
# Input parameters: a and b. 
# Returns boolean of True or False depending on threshold comparison.

def breakpoint(a, b, threshold=20):
    a, b = a[1], b[1]
    return (b - a) > threshold

# Function window:
# Determines window based on total distance covered by window.
# Dependent on a threshold comparison.

def window(a, b, threshold=10):
    a, b = a[1], b[1]
    return b > (a + threshold)

# Function calculate_slope:
# Input: a list
# Performs an if-else loop across windows: If there is just one marker per window \
# then None is returned. Else, set up the X and Y variables for the window and calculate the slope.

def calculate_slope(window):
    group = list(window)
    if len(group) == 1:
        return [None]
    else:
        (Y1, X1), *_, (Y2, X2) = group
        return [(Y2 - Y1) / (X2 - X1)] * len(group)

# Function single_slope:
# Same as calculate_slope except it only prints out one value per window.

def single_slope(window):
    A
    group = list(window)
    if len(group) == 1:
        return [None]
    else:
        (Y1, X1), *_, (Y2, X2) = group
        return [(Y2 - Y1) / (X2 - X1)]

# Function calculate_midpoint:
# Calculates midpoint of position values for a window.

def calculate_midpoint(window):
    group = list(window)
    if len(group) == 1:
        return [None]
    else:
        (Y1, X1), *_, (Y2, X2) = group
        return [(X1 + X2)/2]

# Function pickfirst:
# Picks the first element from a tuple. 

def pickfirst(tup):
    return tup[1]

# Function calculate_mean:
# Calculates mean of position values for a window.

def calculate_mean(window):
    group = list(window)
    if len(group) == 1:
        return [None]
    else:
        first = list(map(pickfirst, group))
        return [np.mean(first)]


# Function windowSlope:
# Calculates slope per pairwise marker comparison in a stepwise fashion

def slidingStep(window):
    longlist = list(window)
    slope_list = []
    if len(longlist) == 1:
        return [None]
    else:
        for i in range(len(longlist)):
            if i + 2  > len(longlist):
                break
            else:
                Y1 = longlist[i][0]
                #print(Y1)
                Y2 = longlist[i+1][0]
                #print(Y2)
                X1 = longlist[i][1]
                #print(X1)
                X2 = longlist[i+1][1]
                #print(X2)
                try: 
                    slope = ((Y2 - Y1) / (X2 - X1))
                except: 
                    slope = 0 
                slope_list.append(slope)
                #print(slope_list)
    return [np.mean(slope_list)]

## MAIN ##

df = pd.read_csv(args.input, delimiter=r"\s+")


# Use enumerate and groupby nested structure to slice windows based on key.
# Use a nested for loop to run through all windows. 
# Use enumerate to get a counter and the value from the iterable at the same time!
# Groupby splits all records from dataset into different categories or groups. 
# Use zip to create tuples of information from the centiMorgan and Position columns. 
# For the same window, use the function window to cut into windows. 
# Create a list "StepList" to place all values. 

stepList = [
   (i, v)
        for i, (dummy, g) in enumerate(
        groupby(zip(df["centiMorgan"], df["Position"]), key=cmp_to_key(window))
        )
        for v in slidingStep(g)
]

# Convert stepList into pandas dataframe: 
stepDF = pd.DataFrame(stepList)

# Add a header to the columns of the dataframe:
stepDF.columns = ['index', 'slope']

# Print out some messages for the user:
print("\n---\n")
print("Printing the calculated average slope per window in a dataframe...")
print(stepDF)

# Calculate midpoint values from the same window:

midList = [
   (i, v)
        for i, (dummy, g) in enumerate(
        groupby(zip(df["centiMorgan"], df["Position"]), key=cmp_to_key(window))
        )
        for v in calculate_midpoint(g)
]
midDF = pd.DataFrame(midList)
midDF.columns = ['index', 'midpoint']
#print("\n---\n")
print("Calculating raw midpoints calculated from windows of position values...") 
#print(midDF)

# Calculate mean values from the same window:

meanList = [
   (i, v)
        for i, (dummy, g) in enumerate(
        groupby(zip(df["centiMorgan"], df["Position"]), key=cmp_to_key(window))
        )
        for v in calculate_mean(g)
]
meanDF = pd.DataFrame(meanList)
meanDF.columns = ['index', 'mean']
#print("\n---\n")
print("Calculating raw means calculated from windows of position values...")
#print(meanDF)


# Combine the slopes, mean and midpoint values into one text file and write to file:
print("\n---\n")
print("Printing slope, mean, and midpoint values calculated from windows of position and centiMorgan values...")
combinedDF = pd.concat([stepDF["slope"], meanDF["mean"], midDF["midpoint"]], axis = 1)
print(combinedDF)
print("\n---\n")
print("The combined mean, midpoint, and slopes file was written to tripleStats.txt.")
combinedDF.to_csv(args.output + "tripleStats.txt", sep='\t', header=True, index=False, na_rep="NaN")
print("\n---\n")
print("Congratulations! You've successfully created all files needed to plot recombination rate across scaffolds! With great power, comes great responsibility.")
