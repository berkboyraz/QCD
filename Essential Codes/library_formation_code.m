clear all
close all
clc

T= 300; % [kelvin]
applied_voltage=0; % [V]

[psic, Ec, z] = SchrodingerPoisson1D_CB_Kane_Main(T, applied_voltage);


