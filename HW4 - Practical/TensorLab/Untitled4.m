clear all; clc; close all;
% Generate pseudorandom factor matrices U0 and their associated full tensor T.
size_tens = [7 8 9]; R = 4;
U = cpd_rnd(size_tens,R);
T = cpdgen(U);
% Compute the CPD of the full tensor T.
Uhat = cpd(T,R);

