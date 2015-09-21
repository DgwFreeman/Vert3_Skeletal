%%          mean thick thin               mean XB strain      ton
%stem3(mean(temp(CompRuns,(9:10)), 2), temp(CompRuns,11), temp(CompRuns,3), 'ko')

Runs=data.data;
%% Now sort out all of our variables, for how many runs to look at for each
%% set.  Here we have our Run Index (col 1) startin over each time
%% we lanuch a new force trace:
%%
NewTraceCol=1; %% column 1
NewTraceIndex=find(Runs(:, NewTraceCol)==1); %% This will give us a handle
%% to sort on each Run as a new trace startes for keeping track of ton for
%% each trace (once it settles out in the end)
ton_raw=zeros(length(NewTraceIndex), 1);
ton_full=zeros(length(NewTraceIndex), 1); %% For Full Cycle
ton_back=zeros(length(NewTraceIndex), 1); %% Did not complete full cycle
N_raw=zeros(length(NewTraceIndex), 1); % Number of Raw
Pct_full=zeros(length(NewTraceIndex), 1); % Percent of Full
Pct_back=zeros(length(NewTraceIndex), 1); %% Percent of backward exit
% and exited out the back door.
XBF_raw=zeros(length(NewTraceIndex), 1); %% XB Force (mean)
XBF_full=zeros(length(NewTraceIndex), 1); %% XB Force Full Cycle
XBF_back=zeros(length(NewTraceIndex), 1); %% XB Force Back
XBNRG_raw=zeros(length(NewTraceIndex), 1); %% XB Energy (mean)
XBNRG_full=zeros(length(NewTraceIndex), 1); %% XB Energy Full Cycle
XBNRG_back=zeros(length(NewTraceIndex), 1); %% XB Energy Back
Total_XBNRG_raw=zeros(length(NewTraceIndex), 1); %% Total XB Energy (sum)
Total_XBNRG_full=zeros(length(NewTraceIndex), 1); %% Total XB Energy Full Cycle
Total_XBNRG_back=zeros(length(NewTraceIndex), 1); %% Total XB Energy Back

%% Some handles
i_ton=3;
iComplete=18;%% The column for completed runs is 18
ifull=8; % If the 8th column is 1, then it went throgh full cycle
%% otherwise == 0 for going out of the back end
iXBstrain=11;
ikXB=16;


%% Loop through each trace for a given Ca2+ level
for i=1:length(NewTraceIndex)
    if i==length(NewTraceIndex)
        RunsThisTrace=Runs([NewTraceIndex(i,1):1:length(Runs)],:); %% Gather All Runs for this trace
        index_CompRuns=find(RunsThisTrace(:,iComplete)==1); %% Get the index for all completed runs this trace
        CompRunsThisTrace=RunsThisTrace(index_CompRuns,:); %% Gather All Completed Runs for this trace
        FullRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==1),:); %% Gather All Completed Runs that finish the cycle
        BackRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==0),:); %% Gather All Completed Runs that exited backward out the cycle
        ton_raw(i,1)=mean(CompRunsThisTrace(:,i_ton)); %% Gather the time on for all completed runs
        ton_full(i,1)=mean(FullRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
        ton_back(i,1)=mean(BackRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
        N_raw(i,1)=length(CompRunsThisTrace); % Number of Raw
        Pct_full(i,1)=length(FullRunsThisTrace)/length(CompRunsThisTrace); %% Percent
        Pct_back(i,1)=length(BackRunsThisTrace)/length(CompRunsThisTrace); %% Percent
        XBF_full(i,1)=mean(FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB));
        XBF_back(i,1)=mean((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)); %% Some error in accounting if not exit out back--not divided sum yet
        XBF_raw(i,1)=mean( [FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB); ...
            (BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)] ); %% Some error in accounting if not exit out back--not divided sum yet
        XBNRG_full(i,1)=mean((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
        XBNRG_back(i,1)=mean(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
        XBNRG_raw(i,1)=mean( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
            ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
        Total_XBNRG_full(i,1)=sum((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
        Total_XBNRG_back(i,1)=sum(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
        Total_XBNRG_raw(i,1)=sum( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
            ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
    else
        RunsThisTrace=Runs([NewTraceIndex(i,1):1:NewTraceIndex(i+1, 1)-1],:); %% Gather All Runs for this trace
        index_CompRuns=find(RunsThisTrace(:,iComplete)==1); %% Get the index for all completed runs this trace
        CompRunsThisTrace=RunsThisTrace(index_CompRuns,:); %% Gather All Completed Runs for this trace
        FullRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==1),:); %% Gather All Completed Runs that finish the cycle
        BackRunsThisTrace=CompRunsThisTrace(find(CompRunsThisTrace(:,ifull)==0),:); %% Gather All Completed Runs that exited backward out the cycle
        ton_raw(i,1)=mean(CompRunsThisTrace(:,i_ton)); %% Gather the time on for all completed runs
        ton_full(i,1)=mean(FullRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
        ton_back(i,1)=mean(BackRunsThisTrace(:,i_ton)); %% Gather the time on for all Full runs
        N_raw(i,1)=length(CompRunsThisTrace); % Number of Raw
        Pct_full(i,1)=length(FullRunsThisTrace)/length(CompRunsThisTrace); %% Percent
        Pct_back(i,1)=length(BackRunsThisTrace)/length(CompRunsThisTrace); %% Percent
        XBF_full(i,1)=mean(FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB));
        XBF_back(i,1)=mean((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)); %% Some error in accounting if not exit out back--not divided sum yet
        XBF_raw(i,1)=mean( [FullRunsThisTrace(:,iXBstrain).*FullRunsThisTrace(:,ikXB); ...
            (BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).*BackRunsThisTrace(:,ikXB)] ); %% Some error in accounting if not exit out back--not divided sum yet
        XBNRG_full(i,1)=mean((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
        XBNRG_back(i,1)=mean(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
        XBNRG_raw(i,1)=mean( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
            ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
        Total_XBNRG_full(i,1)=sum((FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB));
        Total_XBNRG_back(i,1)=sum(((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB));
        Total_XBNRG_raw(i,1)=sum( [(FullRunsThisTrace(:,iXBstrain).^2).*FullRunsThisTrace(:,ikXB); ...
            ((BackRunsThisTrace(:,iXBstrain)./BackRunsThisTrace(:,i_ton)).^2).*BackRunsThisTrace(:,ikXB)] );
    end
end


OutData=[ton_raw, N_raw, ton_full, Pct_full, ton_back, Pct_back, XBF_raw, XBF_full, XBF_back, XBNRG_raw, XBNRG_full, XBNRG_back, Total_XBNRG_raw, Total_XBNRG_full, Total_XBNRG_back];

OutData_Avg=[mean(ton_raw), std(ton_raw),... 
    mean(N_raw), std(N_raw), ...
    mean(ton_full), std(ton_full),...
    mean(Pct_full), std(Pct_full),...
    mean(ton_back), std(ton_back),...
    mean(Pct_back), std(Pct_back), ...
    mean(XBF_raw), std(XBF_raw),...
    mean(XBF_full), std(XBF_full), ...
    mean(XBF_back), std(XBF_back), ...
    mean(XBNRG_raw), std(XBNRG_raw), ...
    mean(XBNRG_full), std(XBNRG_full), ...
    mean(XBNRG_back), mean(XBNRG_back), ...
    mean(Total_XBNRG_raw), std(Total_XBNRG_raw), ...
    mean(Total_XBNRG_full), std(Total_XBNRG_full), ...
    mean(Total_XBNRG_back), std(Total_XBNRG_back)];

clf(figure(10))
subplot(2,2,1)
errorbar(1, mean(ton_raw), std(ton_raw), 'ko', 'markerfacecolor', 'k');, hold on
errorbar(2, mean(ton_full), std(ton_full), 'ko', 'markerfacecolor', 'g');
errorbar(3, mean(ton_back), std(ton_back), 'ko', 'markerfacecolor', 'm');
set(gca, 'box', 'off')
ylabel('t_{on} (ms)')
xlim([0,4])
set(gca, 'xtick', [1,2,3])
set(gca, 'xticklabel', {'All Events', 'Full Cycle', 'Back Detach'})
title(['Set: ' num2str(SetNumber) ', pCa=' pCa_string])

subplot(2,2,2)
%errorbar(1, mean(ton_raw), std(ton_raw), 'ko', 'markerfacecolor', 'k');, hold on
errorbar(2, mean(Pct_full), std(Pct_full), 'ko', 'markerfacecolor', 'g');, hold on
errorbar(3, mean(Pct_back), std(Pct_back), 'ko', 'markerfacecolor', 'm');
set(gca, 'box', 'off')
ylabel('Num. Events (percent)')
xlim([0,4])
set(gca, 'xtick', [1,2,3])
set(gca, 'xticklabel', {'All Events', 'Full Cycle', 'Back Detach'})
title(['N Total/trace ~' num2str(mean(N_raw)) '     ' date])

subplot(2,2,3)
errorbar(1, mean(XBF_raw), std(XBF_raw), 'ko', 'markerfacecolor', 'k');, hold on
errorbar(2, mean(XBF_full), std(XBF_full), 'ko', 'markerfacecolor', 'g');
errorbar(3, mean(XBF_back), std(XBF_back), 'ko', 'markerfacecolor', 'm');
set(gca, 'box', 'off')
ylabel('XB Force (pN)')
xlim([0,4])
set(gca, 'xtick', [1,2,3])
set(gca, 'xticklabel', {'All Events', 'Full Cycle', 'Back Detach'})

subplot(2,2,4)
errorbar(1, mean(XBNRG_raw), std(XBNRG_raw), 'ko', 'markerfacecolor', 'k');, hold on
errorbar(2, mean(XBNRG_full), std(XBNRG_full), 'ko', 'markerfacecolor', 'g');
errorbar(3, mean(XBNRG_back), std(XBNRG_back), 'ko', 'markerfacecolor', 'm');
set(gca, 'box', 'off')
ylabel('XB Energy (pN nm)')
xlim([0,4])
set(gca, 'xtick', [1,2,3])
set(gca, 'xticklabel', {'All Events', 'Full Cycle', 'Back Detach'})
