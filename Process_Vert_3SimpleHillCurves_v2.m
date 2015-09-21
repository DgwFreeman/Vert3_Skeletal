%% BCWT June 3, 2009
%% Version Updates:
%% March 28, 2008
%% Examine the Hill fits, and coefficients of the Hill fits.  Refit the
%% data and determine the standard errors.

function []=Process_Vert_3SimpleHillCurves_v2(pCal, OutDirectory, Data, FigureNumber, Set, Tn_Fraction, XB_Fraction)

%% calculate time step
dt=Data{1,2}(2,1);
Raw_N_Bridges=4*3*60; %% (NMYO*N_Thick_Start*XB_Nodes)

%% Gather force-pCa, and plot
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];

%% create column number to look for data--With the Titin added, we will do
%% the SSF (active force) measurement/estimate with the Thin Filament data;
%% Then we will do the Total Force Eastimate with the Thick Filament Data.
iCol=3;%--Thick Force (Total Force)

% Initialize for the Thin Fil force stuff
iCol2=5;%--Thin Fil Force (Active Force)
Anal_Data_Means2=pCa;
Anal_Dat_SD2=pCa;
Anal_Data2=[];

for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
%% Gather Total Force
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)); %% Gather mean for Data
    Anal_Data=[Anal_Data; Data{i,3}(:,iCol)]; %% Gather Data for indexed and fitting
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)); %% Gather STD DEV of means
%% Gather Active Force    
    Anal_Data_Means2(i, 1)=mean(Data{i,3}(:,iCol2)); %% Gather mean for Data
    Anal_Data2=[Anal_Data2; Data{i,3}(:,iCol2)]; %% Gather Data for indexed and fitting
    Anal_Data_SD2(i, 1)=std(Data{i,3}(:,iCol2)); %% Gather STD DEV of means

end

% disp('BEGIN ANALYSIS:: SSForce')
% Do 4 parameter on parameter on Thick Fil
[head, coeff]=Process_4ParamHill_Anal_v2(Anal_pCa, Anal_Data);

% Do three parameter on Thin Fil
[head2, coeff2]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data2);

%% Plot fit
clf(figure(FigureNumber))
subplot(2,3,1)
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko'), hold on
errorbar(pCa, Anal_Data_Means2, Anal_Data_SD2, 'k*')
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1)))))+coeff(4,1) , 'g-')
plot([4:0.05:9], (coeff2(1,1)./(1+10.^(coeff2(2,1)*([4:0.05:9] - coeff2(3,1))))) , 'b-')
legend(['Max=' num2str(coeff(1,1), 5) '; n_H=' num2str(coeff(2,1), 3), '; pCa_{50}=' num2str(coeff(3,1), 3) '; Offset=' num2str(coeff(4,1), 3) ], ['Max=' num2str(coeff2(1,1), 5) '; n_H=' num2str(coeff2(2,1), 3), '; pCa_{50}=' num2str(coeff2(3,1), 3) ], 'location', 'Northwest')
title(['Set ' num2str(Set) '; ' date])
ylabel('Predicted Force (pN)')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSForce_Total.txt', OutDirectory);
OutFile2=sprintf('%sProcessed_HillData_SSForce_Active.txt', OutDirectory);

% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
%% Write Total Force
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file
%% Write Active Force
fid=fopen(OutFile2,'w');	%open outfile--tab delimited text
fprintf(fid, formatHead);	%header for outfile
fprintf(fid, FormatString, coeff2');
fclose(fid);     %close the file


%% Gather Fraction available-pCa, and plot
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];

%% create column number to look for data
iCol=15;
for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)); %% Gather mean for Data
    Anal_Data=[Anal_Data; Data{i,3}(:,iCol)]; %% Gather Dad for indexed
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)); %% Gather STD DEV of means
end

[head, coeff]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data);
%% Plot fit
figure(FigureNumber)
subplot(2,3,2)
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1))))) , 'g-'), hold on
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko')
legend(['Fit: SS_{max}=' num2str(coeff(1,1), 5)], ['pCa_{50}=' num2str(coeff(3,1), 3) '; n_H=' num2str(coeff(2,1), 3) ], 'location', 'Northwest')
title(['Max Tn Fract = ' num2str(Tn_Fraction)])
ylabel('Fraction Available')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')

% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSFAvail.txt', OutDirectory);

% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file


%% Gather ATP per time-step versus pCa, and plot
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];
%% create column number to look for data
iCol=17;
for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)); %% Gather mean for Data
    Anal_Data=[Anal_Data; Data{i,3}(:,iCol)]; %% Gather Dad for indexed
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)); %% Gather STD DEV of means
end
[head, coeff]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data);
%% Plot fit
figure(FigureNumber)
subplot(2,3,3)
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1))))) , 'g-'), hold on
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko')
legend(['Fit: SS_{max}=' num2str(coeff(1,1), 5)], ['pCa_{50}=' num2str(coeff(3,1), 3) '; n_H=' num2str(coeff(2,1), 3) ], 'location', 'Northwest')
title(['dt calculated = ' num2str(dt) ' s'])
ylabel('ATP per time-step')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')
% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSATP_dt.txt', OutDirectory);
% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file

%% Gather Fraction XB Bound versus pCa, and plot
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];
%% create column number to look for data
iCol=19;
for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)); %% Gather mean for Data
    Anal_Data=[Anal_Data; Data{i,3}(:,iCol)]; %% Gather Dad for indexed
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)); %% Gather STD DEV of means
end
[head, coeff]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data);
%% Plot fit
figure(FigureNumber)
subplot(2,3,4)
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1))))) , 'g-'), hold on
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko')
legend(['Fit: SS_{max}=' num2str(coeff(1,1), 5)], ['pCa_{50}=' num2str(coeff(3,1), 3) '; n_H=' num2str(coeff(2,1), 3) ], 'location', 'Northwest')
title(['Raw Calc; N_{Bridges}= ' num2str(Raw_N_Bridges)])
ylabel('Fract XB Bound')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')
% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSFBnd_Raw.txt', OutDirectory);
% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file

%% Gather Fraction XB Bound versus pCa, and plot
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];
%% create column number to look for data
iCol=19;
for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)*(1/XB_Fraction)); %% Gather mean for Data
    Anal_Data=[Anal_Data; (Data{i,3}(:,iCol)*(1/XB_Fraction))]; %% Gather Dad for indexed
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)*(1/XB_Fraction)); %% Gather STD DEV of means
end
[head, coeff]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data);
%% Plot fit
figure(FigureNumber)
subplot(2,3,5)
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1))))) , 'g-'), hold on
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko')
legend(['Fit: SS_{max}=' num2str(coeff(1,1), 5)], ['pCa_{50}=' num2str(coeff(3,1), 3) '; n_H=' num2str(coeff(2,1), 3) ], 'location', 'Northwest')
title(['Scaled Calc; XB Fraction = ' num2str(XB_Fraction)])
ylabel('Scaled XB Bound')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')
% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSFBnd_Scaled.txt', OutDirectory);
% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file



%% Gather ATP per time-step versus pCa, and plot; Scale this to ATPase
%% Scale ATP per Time step -- into ATP*myosin^-1second^-1
pCa=zeros(length(pCal), 1);
Anal_Data_Means=pCa;
Anal_Data_SD=pCa;
Anal_pCa=[];
Anal_Data=[];
%% create column number to look for data
iCol=17;
for i=1:length(pCal)
    pCa(i, 1)=mean(Data{i,3}(:,2)); %% Gather mean pCa
    Anal_pCa=[Anal_pCa; Data{i,3}(:,2)]; %% Gather pCa indexed
    Anal_Data_Means(i, 1)=mean(Data{i,3}(:,iCol)/(dt*Raw_N_Bridges*XB_Fraction)); %% Gather mean for Data
    Anal_Data=[Anal_Data; (Data{i,3}(:,iCol)/(dt*Raw_N_Bridges*XB_Fraction))]; %% Gather Dad for indexed
    Anal_Data_SD(i, 1)=std(Data{i,3}(:,iCol)/(dt*Raw_N_Bridges*XB_Fraction)); %% Gather STD DEV of means
end
[head, coeff]=Process_3ParamHill_Anal_v2(Anal_pCa, Anal_Data);
%% Plot fit
figure(FigureNumber)
subplot(2,3,6)
plot([4:0.05:9], (coeff(1,1)./(1+10.^(coeff(2,1)*([4:0.05:9] - coeff(3,1))))) , 'g-'), hold on
errorbar(pCa, Anal_Data_Means, Anal_Data_SD, 'ko')
legend(['Fit: SS_{max}=' num2str(coeff(1,1), 5)], ['pCa_{50}=' num2str(coeff(3,1), 3) '; n_H=' num2str(coeff(2,1), 3) ], 'location', 'Northwest')
title(['Scaled ATPase'])
ylabel('ATP myosin^{-1} s^{-1}')
xlabel('pCa')
set(gca, 'xdir', 'reverse')
set(gca, 'box', 'off')
% Add path to access the write out the file portion...
OutFile=sprintf('%sProcessed_HillData_SSATPase_Scaled.txt', OutDirectory);
% write data file for parameters from fit
% format the header for writing the output file
formatHead=[];
FormatString=[];
for i=1:length(head)-1
    formatHead=[formatHead, head{1,i}, '\t'];
    FormatString=[FormatString, '%6.4f\t'];
end
formatHead=[formatHead, head{1,end}, '\n'];
FormatString=[FormatString, '%6.4f\n'];
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, formatHead);	%header for outfile
% Write Data File and close file
fprintf(fid, FormatString, coeff');
fclose(fid);     %close the file

%% Write the quick outfile
OutFig=sprintf('%sSSHillFits.fig', OutDirectory);
saveas(figure(FigureNumber), OutFig, 'fig')

%% Quickly plot some time series data
clf(figure(FigureNumber+1))
for i=1:length(pCal)
    %% Plot force versus time at each pCa
    subplot(1,2,1)
    plot(Data{i,2}(:,1), Data{i,2}(:,2), 'k-'), hold on
    %% Plot fraction avail versus time at each pCa
    subplot(1,2,2)
    plot(Data{i,2}(:,1), Data{i,2}(:,10), 'b-'), hold on
end

subplot(1,2,1)
title(['Set ' num2str(Set) '; ' date])
ylabel('Force (pN)')
xlabel('Time (s)')
set(gca, 'box', 'off')

subplot(1,2,2)
title(['Set ' num2str(Set) '; ' date])
ylabel('Fraction Available')
xlabel('Time (s)')
set(gca, 'box', 'off')

%% Write the quick outfile
OutFig=sprintf('%sTimeSeries_Force_FA.fig', OutDirectory);
saveas(figure(FigureNumber+1), OutFig, 'fig')
