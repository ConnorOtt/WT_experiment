% Script for testing WT_experiment
clear; close all; clc;

dataSet = load('data/ExperimentalLab2_Section3_Group06_Even_Long.csv');
m = WT_experiment(dataSet, 500);
b = WT_experiment(dataSet, 500);

figure
plot(m, 'AoA', 'drag', 'rx')
title('Angle of Attack vs. Drag')
xlabel('Angle of Attack, [$^{\circ}$]')
ylabel('Drag [N]')