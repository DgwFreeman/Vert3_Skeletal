function [k01,k10,k12,k21,k20,k02] = XBRates4(xin, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD, S1, S2, S3)


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

xinr = xin - reach;
kRTs = sqrt(kRT);

%  /* preliminary computations */
G0 = 0.0;
G1 = f1 * dGtot + kRT * (xinr * xinr);
G2 = eta * dGtot + kRT * (xin * xin);

G0e = exp(G0);
G1e = exp(G1);
G2e = exp(G2);

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
k01 = (S1 * 2000.0) * kRTs * exp( -1/2 * kRT * (xinr * xinr)) / sqrt(2.0 * pi);
%/*k01 = k01 * S1[0];*/
%output[0] = k01;

%/* k10 */
%output = mxGetPr(plhs[1]);
k10 = k01 * G1e / G0e;

%   /* k12 */
%   output = mxGetPr(plhs[2]);
k12 = S2 * ( AA + (100 / kRTs) * (1 - tanh(kRTs * xinr)) );
%  output[0] = k12;

%   /* k21 */
%   output = mxGetPr(plhs[3]);
k21 = k12 * G2e / G1e;

%   /* k20 */
%   output = mxGetPr(plhs[4]);
k20 = S3 * ((kRTs * (sqrt(BB) * xin + CC * xin)) + DD);
%  output[0] = k20;

%  /* k02 */
%  output = mxGetPr(plhs[5]);
k02 = k20 * exp(dGtot) * G0e / G2e;

% }
%
% /*****************************************************/
