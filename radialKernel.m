%% For the Class of Radial Basis Kernels
%% k(xi,xj) = exp(-(||xi - xj|| ^ 2) / 2 * (sigma ^ 2)); where sigma is a specified value

function [] = radialKernel(inputFile, sigma)
	%% Reading the data file optdigit.data file %%
	if strcmp(inputFile,'optdigit.data')
		X = dlmread('optdigit.data', ' ', 1, 0 );
	elseif strcmp(inputFile,'iris.data')
		X = dlmread('iris.data', ',', 2, 0 );
	end

	%% Finding the dimensions of the matrix X
	[nRows nCols] = size(X);
	X = X(:,1:nCols - 1);

	%% Centralizing the matrix X around its mean
	X = bsxfun(@minus, X, mean(X));
	
	%%Forming the kernel matrix
	Kernel = zeros(nRows , nRows);
	for i = 1 : nRows
		for j = 1 : nRows
			Kernel(i,j) = exp(-(norm(X(i, :) - X(j, :)) / (2 * sigma * sigma)));
		end
	end
	%Kernel
	
	%%Finding the eigenDecomposition of the kernel matrix and sorting in descending order
	[eigVector eigValue] = eig(Kernel);
	[eigValue, order] = sort(diag(eigValue), 'descend');
	eigVector(:,order) = eigVector; 
	
	%----------------Finding the projections---------------------
	projecMatriX = Kernel * eigVector;
	
	h = figure;
	temp1 = zeros(nRows, 2);
	temp1(:,1) = projecMatriX(:,1);
	temp1(:,2) = projecMatriX(:,2);
	
	for i = 1: nRows
		p1(i) = temp1(i,1);
		p2(i) = temp1(i,2);
		
	end
	
	%% Plotting the projections
	plot(p1,p2, 'ro');
	
	if strcmp(inputFile,'optdigit.data')
		imageName = 'radialKernelProjectionOptDigit';
	elseif strcmp(inputFile,'iris.data')
		imageName = 'radialKernelProjectionIris';
	end
	saveas(h, imageName,'png');
	
	%% Plotting the eigen Variations of the Kernel
	[rowKernel colKernel] = size(Kernel);
	xCoordinate2 = 1: 1: colKernel;
	for i = 1: colKernel
		yCoordinate2(i) = eigValue(i);
	end
	 
	h = figure;
	figureEigenVariation2 = plot(xCoordinate2, yCoordinate2);
	
	if strcmp(inputFile,'optdigit.data')
		imageName = 'radialKernelEigenValueOptDigit';
	elseif strcmp(inputFile,'iris.data')
		imageName = 'radialKernelEigenValueIris';
	end
	saveas(h,imageName,'png');
end