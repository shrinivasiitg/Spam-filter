%% This function works on "optdigit.data set" kept in the current matlab directory
%% WARNING: optdigit.data should be kept in the current folder, else the program will crash.

function[] = testPca(inputFile)
	%% Reading the data file optdigit.data file %%
	if strcmp(inputFile,'optdigit.data')
		X = dlmread('optdigit.data', ' ', 1, 0 );
	elseif strcmp(inputFile,'iris.data')
		X = dlmread('iris.data', ',', 2, 0 );
	end

	[nRows, nCols] = size(X);
	X = X(:, 1:nCols - 1);

	%% Centralizing the data around the mean %%
	X = bsxfun(@minus, X, mean(X));

	%% Singular value decomposition 
	[U, S, V] = svd(X);

	%%Principal directions
	dir1 = V(:, 1);
	dir2 = V(:, 2);

	%%First two singular values
	tempFirst2 = S(:,1);
	tempSecond2 = S(:,2);

	principalComponent1 = U * tempFirst2;
	principalComponent2 = U * tempSecond2;
	
	%% Plotting the projected data set
	h= figure;
	plot(principalComponent1,principalComponent2, 'bo');
	saveas(h,'projectedX','png');
	
	%% Plotting the biplots
	%%%%%%%%% plot for beta = 0 %%%%%%%%%%%%%%%%%%%%
	%% if beta = 0 then G = U, H = S * V'
	%take the first two columns of G and first two rows of H
	
	%% Initializing the G matrix with zeros
	 G = zeros(nRows, 2);
	 
	 %Take the first two columns of G
	 G(:,1) = U(:,1);
	 G(:,2) = U(:,2);
	 
	 %Take the first two rows of H
	 H = S * V';
	 colorVector = 'rgbcmykw';
	 index = 1;

	h = figure;
	for i = 1 : nRows
		temp1 = [0, G(i,1)];
		temp2 = [0, G(i,2)];
		
		 if index == 8
			 index = 1;
		 end
		  
		%% Doing the beta plot for G
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end

	for i = 1 : nCols - 1
		temp1 = [0, H(1,i)];
		temp2 = [0, H(2,i)];
		
		 if index == 8
			 index = 1;
		 end
		 
		 %% Doing the beta plot for G
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end
	saveas(h, 'biPlot2Beta_0','png');
	hold off;

	% %%%%%%%%% plot for beta = 0.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%if beta = 0.5 then G = U * (S ^ 0.5), H = (S ^ 0.5) * V'
	%take the first two columns of G and first two rows of H
	 
	 %% Saving the S.^0.5 matrix as a temporary matrix
	 
	 tempS = S.^0.5;
	 G = U * tempS;
	 H = tempS * V';
	 colorVector = 'rgbcmykw';
	 index = 1;
	h = figure; 
	for i = 1 : nRows
		temp1 = [0, G(i,1)];
		temp2 = [0, G(i,2)];
		
		 if index == 8
			 index = 1;
		 end
		 
		%% Doing the beta plot for G 
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end

	for i = 1 : nCols - 1
		temp1 = [0, H(1,i)];
		temp2 = [0, H(2,i)];
		
		 if index == 8
			 index = 1;
		 end
		%% Doing the beta plot for H  
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end

	hold off
	saveas(h, 'biPlot2Beta_Half','png');

	% %%%%%%%%% plot for beta = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% for Beta = 1, G = U * S, H = V'

	G = U * S;
	H = V';
	colorVector = 'rgbcmykw';
	 index = 1;
		h = figure;
	for i = 1 : nCols - 1
		temp1 = [0, G(i,1)];
		temp2 = [0, G(i,2)];
		
		 if index == 8
			 index = 1;
		 end
		%% Doing the beta plot for G  
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end

	for i = 1 : nCols - 1
		temp1 = [0, H(1,i)];
		temp2 = [0, H(2,i)];
		
		 if index == 8
			 index = 1;
		 end
		%% Doing the beta plot for H
		plot(temp1, temp2, colorVector(index)); hold on
		index = index + 1;
	end

	hold off
	saveas(h, 'biPlot2Beta_1','png');
	
	%Reconstruction of the data points
	%Reconstruct using various direction values
	
	xSample2 = 1: 1 : nRows;
	tempS = zeros(nCols - 1, nCols - 1);
	deviationError2 = zeros(1,nRows);
	
	%% Outer loops regenerates the original matrix at an interval of 4
	for i = 1 : nCols - 1
		tempS(:,i) = V(:, i);
			if mod(i,4) == 0 
				%%Reconstructing the data
				ReconstructedX = U * S * tempS';
				%Finding the deviation error for each data point
				deviation2 = (ReconstructedX - X);
				deviationError2 = zeros(1,nRows);
				
				%Find error for each data point
				for j = 1: nRows
					 deviationError2(j) = norm(deviation2(j,:));
				end
				h = figure;
				plot(xSample2, deviationError2);
				fileName = 'deviationErrorX';
				ext = int2str(i);
				imageName = strcat(fileName, ext);
				saveas(h,imageName,'png');
			end
	end
	
	%finding the eigen value variations
	xCoordinate2 = 1: 1: nCols - 1;
	eigValue2 = sum(S.*S,2);
	for i = 1: nCols - 1
		yCoordinate2(i) = eigValue2(i);
	end
	 
	h = figure;
	figureEigenVariation2 = plot(xCoordinate2, yCoordinate2);
	saveas(h,'eigenValueVariation2','png');

end		%End of the function