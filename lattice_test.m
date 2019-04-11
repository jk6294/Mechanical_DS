%% Prepare Space
clear; clc;


%% Create Lattice
A = zeros(1000);
nInd = 1:10:1000;
A(nInd,nInd) = 1;

% Compatibility



%% Fourier Transform
AF = fft2(A);