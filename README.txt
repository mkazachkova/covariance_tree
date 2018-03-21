Authors: Mariya Kazachkova and Michael Mudgett
Course: CIS I
Assignment: Programming Assignment 4
Date: 11/30/2017

This folder contains our code and data for Programming Assignment 4.
Code:
pa4.py - This is the main script used to run the code.
Covariance_Tree.py - Implementation of covariance tree
basic_transforms_piv_calibration.py - This file holds the functions used in our script.
tests.py - Unit tests
compare_output.py - Script that checks user input with provided input

********Script
To run the script, use the command line arguments:
python pa4.py <file prefix> > <output filename>

Where <file prefix> is the start of any of the file names and <output filename> is any made up name for an output file. The script normally prints the data to the terminal but here we can send it into an output file. For example, to run on the debug-a dataset, enter:
python pa4.py PA234_Student_Data/PA4-debug-A > outfile.txt
Note: this assumes the data files are stored in a folder named PA234_Student_Data 

********Unit Tests
To run the unit tests script, simply enter:
python tests.py
And verify that no error messages are produced on the terminal.

********Results Verificiation
To check our output against the provided output.txt files, run the compare_output.py script as follows:
python compare_output.py <provided output file> <user created output file> 
For example:
python compare_output.py PA234_Student_Data/PA4-E-Debug-Output.txt Output/PA4-E-Debug-Result.txt

This will print out the error (difference) between the user output file and the given output file. The error printed is the averaged error between our d_k values and provided d_k values, our c_k values and provided c_k values, and our magnitude of difference and provided magnitude of difference.
