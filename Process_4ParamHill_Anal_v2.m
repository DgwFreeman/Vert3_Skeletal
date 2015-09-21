%% March 28, 2008, BCWT

%% Process 3-parameter hill data, and return the values, with their
%% variance, standard deviation, and the 
%% cofindence intervals (0.05 and 0.01) for the coeficients.
%% Fits with NLINFIT; analyzes statistics with NLPARCI

%% Returns a matrix (columns are labeled, rows depends upon number of
%% variables)  Here row 1 is maxVal; row2 Slope; row3 MidPoint
%% Cols: 
%% 1) coeff for NLINFIT
%% 2) standard deviation
%% 3) variance
%% 4) confidence interval (CI) based on covariance matrix
%% alpha = 0.05, lower bound
%% 5) CI based on covariance matrix, alpha = 0.05, upper bound
%% 6) CI based on covariance matrix, alpha = 0.01, lower bound
%% 7) CI based on covariance matrix, alpha = 0.01, upper bound
%% 8) coeff for LSQCURVEFIT
%% 9) CI based on jacobian, alpha = 0.05, lower bound
%% 10) CI based on jacobian, alpha = 0.05, upper bound
%% 11) CI based on jacobian, alpha = 0.01, upper bound
%% 12) CI based on jacobian, alpha = 0.01, upper bound


function [header, vals]=Process_4ParamHill_Anal_v2(xdata, ydata)

% Create the 3-parm Hill function that you will be fitting the data to
% disp('For reference: maxval=x(1); slope=x(2); K50=x(3); OFFSET=x(4)')
f=@(x,xdata)( (x(1)./(1+10.^(x(2)*(xdata - x(3))))) + x(4) );

% Fit the data using LSQCURVEFIT
maxGuess=max(ydata);
slopeGuess=1+2.2*rand;
midGuess=5.5+2*rand;
offsetGuess=min(ydata);


% disp('Now analyze the data fits using NLINFIT:')
% In addition, there is a suggestion to examine/use NLINFIT
% from Matlab, when using NLPARCI
% So, we can do all this fitting again with a different function
% In addition this throws back the covariance matrix
[beta,resid,J,COVbeta,MSE] = nlinfit(xdata,ydata,f,[maxGuess slopeGuess midGuess offsetGuess]);

%% MOD for TESTS on April 22, 2009:
% disp('Size Resid')
% size(resid)
% disp('Size Beta')
% size(beta)
% disp('Size J')
% size(J)
% disp('Rank J')
% rank(J)
% From Matlab's website:
% The confidence interval calculation is valid for systems where the
% length of resid exceeds the length of beta and J has full column rank.
% When J is ill-conditioned, confidence intervals may be inaccurate.
%% Example Output:
%% Size Resid =  1   144
%% Size Beta =   1     3 %% Length of resid is > than beta
%% Size J =  144     3
%% Rank J =     3 %% J has full column rank; we should be good.

%disp('Use NLPARCI to find confidence intervals on the parameters (covariance matrix)')
ci_cov_05 = nlparci(beta,resid,'covar',COVbeta);
ci_cov_01 = nlparci(beta,resid,'covar',COVbeta, 'alpha', 0.01);

%% 1) coeff for NLINFIT
%% 2) standard deviation
%% 3) variance
%% 4) confidence interval (CI) based on covariance matrix
%% alpha = 0.05, lower bound
%% 5) CI based on covariance matrix, alpha = 0.05, upper bound
%% 6) CI based on covariance matrix, alpha = 0.01, lower bound
%% 7) CI based on covariance matrix, alpha = 0.01, upper bound

header={'coeff', 'std', 'var', 'CI_05_LB', 'CI_05_UB', 'CI_01_LB', 'CI_01_UB'};

vals=[];
vals(1:4,1)=beta';
vals(1:4,2)=[sqrt(COVbeta(1,1));sqrt(COVbeta(2,2));sqrt(COVbeta(3,3));sqrt(COVbeta(4,4))];
vals(1:4,3)=[COVbeta(1,1);COVbeta(2,2);COVbeta(3,3);COVbeta(4,4)];
vals(1:4,4:5)=ci_cov_05;
vals(1:4,6:7)=ci_cov_01;
