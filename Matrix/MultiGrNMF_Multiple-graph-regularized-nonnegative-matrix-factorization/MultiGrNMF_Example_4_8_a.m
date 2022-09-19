clc
close all
clear all

load Input4_8
Output = MultiGrNMF(Input, Options);

figure;
imagesc(Output.W);hold on;
title('W');hold on;
