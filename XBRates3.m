function [k01,k10,k12,k21,k20,k02] = XBRates3(xin, thermochem)

  
%   /* get input data matrices */
%   xin=mxGetPr(i1);
%   reach=mxGetPr(i2);
%   kRT=mxGetPr(i3);
%   f1=mxGetPr(i4);
%   eta=mxGetPr(i5);
%   dGtot=mxGetPr(i6);
%   AA=mxGetPr(i7);
%   BB=mxGetPr(i8);
%   CC=mxGetPr(i9);
%   DD=mxGetPr(i10);
%   S1=mxGetPr(i11);
%   S2=mxGetPr(i12);
%   S3=mxGetPr(i13);
  
%% Ugly and I apologize, I'll be working on something better
% there are probably functions to do this automatically but this function
% is a major bottleneck and shaving microseconds matters, at least until I
% get OpenCL stuff up.
f1 = thermochem.f1;
dGtot = thermochem.dGtot;
kRT = thermochem.kRT;
reach = thermochem.reach;
eta = thermochem.eta;
AA = thermochem.AA;
BB = thermochem.BB;
CC = thermochem.CC;
DD = thermochem.DD;
S1 = thermochem.SXB1;
S2 = thermochem.SXB2;
S3 = thermochem.SXB3;

%%  /* preliminary computations */
  G0 = 0.0;
  G1 = f1 * dGtot + kRT * ( (xin - reach) * (xin - reach) );
  G2 = eta * dGtot + kRT * (xin * xin);
  
%   /* setup the output matrix */
%   plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
%   plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
%   plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
%   plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
%   plhs[4]=mxCreateDoubleMatrix(1,1,mxREAL);
%   plhs[5]=mxCreateDoubleMatrix(1,1,mxREAL);
  
% %   /* complete the computations and finish setting up the output */
% %   /* k01 */
%   output = mxGetPr(plhs[0]);
  k01 = (S1 * 2000.0) * sqrt( kRT / (2.0 * pi)) * exp( -1 * kRT * ((xin - reach) * (xin - reach)) / 2) ;
  %/*k01 = k01 * S1[0];*/ 
  %output[0] = k01;
  
  %/* k10 */
  %output = mxGetPr(plhs[1]);
  k10 = k01 * exp(G1 - G0);
  
%   /* k12 */
%   output = mxGetPr(plhs[2]);
  k12 = S2 * ( AA + (100 / sqrt(kRT)) * (1 - tanh(sqrt(kRT) * (xin-reach))) );
%  output[0] = k12;
  
%   /* k21 */
%   output = mxGetPr(plhs[3]);
  k21 = k12 * exp(G2 - G1);
  
%   /* k20 */
%   output = mxGetPr(plhs[4]);
  k20 = S3 * ( (sqrt(kRT) * (sqrt(BB * (xin * xin)) + CC * xin)) + DD );
%  output[0] = k20;
  
%  /* k02 */
%  output = mxGetPr(plhs[5]);
  k02 = k20 * exp((G0 - G2) + dGtot);

% }
% 
% /*****************************************************/
