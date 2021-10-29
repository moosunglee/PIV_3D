% Run digital image correlation, also called particle image velocimetry.
%
% This file calls the DIC analysis and runs it for the images in the
% current working directory. Data is saved in mat format in the current
% working directory. Main DIC analysis files were written by members of
% Jeff Fredberg's group at Harvard School of Public Health:
%   Xavier Trepat 03/2008
%   Dhananjay Tambe 03/2009
%
% Written by Jacob Notbohm, Harvard School of Public Health, 2013
% Updated by Moosung Lee, KAIST, 2021

clear;
close all;
clc;
cd0 = matlab.desktop.editor.getActiveFilename;
dash = cd0(strfind(cd0,'MAIN')-1);
cd0 = cd0(1:strfind(cd0,'MAIN')-2);
cd(cd0)
addpath(genpath(cd0));

% Initialize arameters
    % Subset size (pix) 
    w0 = 32;
    % Subset spacing (pix) 
    d0 = 6;
    params.blocksizes = [w0 w0 8]; % [iblocksize jblocksize,kblocksize]
    params.overlap = 1-d0/w0;
    params.padding = 16;
    params.N = 2;
    params.Ni = 4;
    params.use_GPU = true;
    TFM_params=TFM_ANALYSIS.get_default_parameters(params);
    TFM=TFM_ANALYSIS(TFM_params);
    
%% Run
load('imgfiles.mat')
Quiver = TFM.PIV_3D(imref,imcur);
%%
figure
subplot(131), quiver(Quiver.X, Quiver.Y, Quiver.U, Quiver.V), axis image,axis off, title('2D Vektor map')
subplot(132), quiver3(Quiver.X, Quiver.Y, Quiver.Z, Quiver.U, Quiver.V,Quiver.W), axis image,axis off, title('3D Vektor map')
subplot(133), imagesc(Quiver.pkh,[0.9 1]), title('My method'), axis image, axis off, colorbar,title('Correlation')
set(gcf,'color','w')


