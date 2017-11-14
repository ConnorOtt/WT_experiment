% Script for testing WT_experiment
clear; close all; clc;

dataSet = load('data/ExperimentalLab2_Section3_Group06_Even_Long.csv');
m = WT_experiment(dataSet, 500);
b = WT_experiment(dataSet, 500);
fake = WT_experiment;
c = meanWT(fake, {m, b});
