%% Enter the file Name as 'optdigit.data' or 'iris.data'

str = input('Enter the file name to be tested : ','s');
testPca(str);


%% Enter the p parameter for the polynomial kernel
p = input('Enter the p value for the polynomial kernel: ');
polyKernel(str,p);


%% Enter the sigma value for the Radial Kernel
sigma = input('Enter the sigma value for the Radial Kernel : ');
radialKernel(str,sigma);