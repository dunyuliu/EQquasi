function [path,dx,nzz,l] = model_info(mod)
%MODEL_INFO Summary of this function goes here
%   This is a function to load model information that are used by all
%   post-processing scripts. 
%   Created on 06/16/2021. 

if mod == 11
    path = '../res/202106-BP5-productive-run/grace/dliu.8_case2_dx500_xi0.01_load4e-10_newst/';
    l = 60;
    dx = 0.5;
    nzz = fix(l/dx)+1; 
elseif mod == 12
    path = '../res/202106-BP5-productive-run/grace/dliu.9_case2_dx500_xi0.015_load4e-10_newst/';
    l = 60;
    dx = 0.5;
    nzz = fix(l/dx)+1;  
elseif mod == 21
    path = '../res/202106-BP5-productive-run/frontera/dliu.6_case1_dx1000_xi0.015_load4e-10_newst/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1; 
elseif mod == 22
    path = '../res/202106-BP5-productive-run/frontera/dliu.7_case1_dx1000_ksi0.05_load4e-10/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1;   
elseif mod == 31
    path = '../res/202106-BP5-productive-run/frontera/dliu.10_case0_dx2000_xi0.015_load4e-10_newst/';
    l = 60;
    dx = 2;
    nzz = fix(l/dx)+1;    
elseif mod == 2
    path = '../res/202104-preliminary4Junle/case2_dx500/';
    l = 60;
    dx = 0.5;
    nzz = fix(l/dx)+1;
elseif mod == 3
    path = '../res/202106-BP5-productive-run/frontera/case1_dx1000/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1;        
elseif mod == 4
    path = '../res/202106-BP5-productive-run/frontera/case2_dx500_2node/';
    l = 60;
    dx = 0.5;
    nzz = fix(l/dx)+1;
elseif mod == 5
    path = '../res/202106-BP5-productive-run/grace/case1_dx1000_dy500/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1; 
elseif mod == 53
    path = '../res/202106-BP5-productive-run/grace/case2_dx500_ksi0.02_load4e-10/';
    l = 60;
    dx = 0.5;
    nzz = fix(l/dx)+1;         
elseif mod == 6
    path = '../res/202106-BP5-productive-run/grace/case1_dx1000_ymax100_load2.5e-10/';
    dx = 1;
    l = 90;
    nzz = fix(l/dx)+1; 
elseif mod == 7
    path = '../res/202106-BP5-productive-run/grace/case1_dx1000_L1.0/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1; 
elseif mod == 8
    path = '../res/202106-BP5-productive-run/grace/case1_dx1000_rat1.05/';
    l = 60;
    dx = 1;
    nzz = fix(l/dx)+1;     
elseif mod == 9
    path = '../res/202106-BP5-productive-run/grace/case1_dx1000_ymax100_load0/';
    dx = 1;
    l=90;
    nzz = fix(l/dx)+1;     
elseif mod == 10
    path = '../res/202106-BP5-productive-run/grace/case2_dx500/';
    dx = 0.5;
    l=60;
    nzz = fix(l/dx)+1;       
elseif mod == 11
    path = '../res/202106-BP5-productive-run/grace/case2_dx500_ymax100_load0/';
    dx = 0.5;
    l=60;
    nzz = fix(l/dx)+1; 
elseif mod == 12
    path = '../res/202106-BP5-productive-run/frontera/BP4QD_case1_dx1000_ksi0.05_load4e-10/';
    dx = 1;
    l=60;
    nzz = fix(2*l/dx)+1;   
end
end

