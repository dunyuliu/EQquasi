function [path,dx,nzz,l] = model_info(mod)
%MODEL_INFO Summary of this function goes here
%   This is a function to load model information that are used by all
%   post-processing scripts. 
%   Created on 06/16/2021. 

if mod == 1
    path = '../bp7-qd-a-10.002/';
elseif mod == 2
    path = '../bp7-qd-a-10.002.053/';
elseif mod == 3
    path = '../bp7-qd-a-10.002.load04/';
end

l = 500/1000*2;
dx = 10/1000;
nzz = fix(l/dx)+1; 

