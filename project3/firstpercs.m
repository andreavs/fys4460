
close all;
clear all;

L = 300;

r = rand(L,L);

p = 0.7;
z = r<p; % This generates the binary array
[lw,num] = bwlabel(z,4);

img = label2rgb(lw,'jet','k','shuffle');
image(img);