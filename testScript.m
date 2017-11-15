% Script for testing WT_experiment
clear; close all; clc;

dataSet = load('data/ExperimentalLab2_Section3_Group06_Even_Long.csv');
m = WT_experiment(dataSet, 500);
b = testParse(m, 'V_inf', [15, 25]);

b{2}.isFinite = 0; % infinite
b{2}.chord = 0.1;    % chord

plot(b{2}, 'AoA', 'liftCoef', 'rx')