% Script for testing WT_experiment
clear; close all; clc;

dataSet = load('data/ExperimentalLab2_Section3_Group06_Even_Long.csv');
m = WT_experiment(dataSet, 500);
m.fileName = 'ExperimentalLab2_Section3_Group06_Even_Long.csv';

m.area = 0.0144;
m.chord = 0.06;
m.momDist = 0.07;
plot(m, 'AoA', 'momCoefAtx', 'rx', 'linewidth', 1)


