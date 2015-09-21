%% BCWT, June 4 2009
%% Try to make a function to calculate time-on for binding on the fly
%% Bingos bangos bongos

%% Inputs:
% Current) thin filament binding state
% Previous) thin filament binding state
% Total) overall set of binding events--start to stop
% XB1) Sub-binding events within an overall total event, only XB state 1
% (lower force bearing)
% XB2) Sub-binding events within an overall total event, only XB state 2
% (higher force bearing)

%% Finalized Total Binding Matrix
% %% Column Headers for Total Binding Event Matrix
%      1 Overall Numerical Index for the run
%      2 Starting simulation time step for this run
%      3 Total Number of time-steps XB bound
%      4 Number of time steps in XB 1-State
%      5 Number of time steps in XB 2-State
%      6 N steps XB1-->XB2
%      7 N steps XB2-->XB1
%      8 N steps XB2-->XB0
%      9 Mean position of XB Node during run
%     10 Mean position of Thin Fil Node during run
%     11 Mean overall bridge strain during run
%     12 Mean XB1 Bridge strain during run
%     13 Mean XB2 Bridge strain during run
%     14  Thick filament node number for event
%     15  Thin filament node number for the even
%     16  kxscaler
%     17  reach (xb0)
%     18 1 or 0; 1 if completed the run
%     19 Variance for position of XB Node during run
%     20 Variance for position of Thin Fil Node during run
%     21 Variance for overall bridge strain during run
%     22 Variance for XB1 Bridge strain during run
%     23 Variance for XB2 Bridge strain during run

function [temp_Prior, Total] = Time_On_Calculator(X, Current, temp_Prior, Total, SS_cutoff, current_step, xb0, kxscaler)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Process the t-on stuff
%%%%%%%%%%%%% Process the t-on stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Now process the Run Index counts
I = find(Current(:,2)>0); %% Grab the rows
Now_Bound_Length=length(I);
%temp_Now, columns:
%Act Node, XB State, XB, Run index, Flag_Checked
temp_Now=[I, Current(I, 2:3), zeros(Now_Bound_Length, 2)];

%% Get some sizing laid out for our columns:
%% Number of columns for Total
BinderCols=23;
%% Number of column for our subset arrays
%SubXBCols=18;

%% Do we have any bridges?  If so, we need to account for them.
if not(isempty(I))
    % If we never had briges before we need to initialize some
    % variables.
    if (isempty(Total))
        %% Update the big matrix for all the current binding
        for iCount=1:Now_Bound_Length  %% Update our Run Matrix
            %            RunCount=RunCount+1; %% Yeah, we have a bridge--update

            %% Make a new row for the Total:
            Total=[Total; zeros(1, BinderCols)];
            Total(iCount,1)=iCount;  %%Index of Run
            temp_Now(iCount, 4)=Total(iCount, 1); %% Update index of run in Temp_Now
            Total(iCount,2)=current_step;  %% current time-step = start time
            %% If this guy was just bound, and not a new binder after the
            %% SS cutoff time, then we don't want to fill the rest of the
            %% array.  This way, we can sort out whether this was a valid
            %% and complete run from start to finish.  We never got the
            %% start of this guy.
            if current_step>SS_cutoff
                Total(iCount,3)=1;   %%Duration in total
                Total(iCount,6)=0; %% Step XB1-->XB2
                Total(iCount,7)=0; %% Step XB2-->XB1
                Total(iCount,8)=0; %% Step XB2-->XB0
                Total(iCount,9)=X(temp_Now(iCount, 3), 1); %% Position of XB Node
                Total(iCount,10)=X(temp_Now(iCount,1), 1); %% Position of Thin Fil Node
                if temp_Now(iCount, 2)==1
                    Total(iCount,4)=1; %% XB 1-State
                    Total(iCount,11)=(X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0; %% Total Bridge Strain
                    Total(iCount,12)=(X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0; %% XB1 Bridge Strain
                else
                    Total(iCount,5)=1; %% XB 2-State
                    Total(iCount,11)=X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1); %% Total Bridge Strain
                    Total(iCount,13)=X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1); %% XB2 Bridge Strain
                end
                Total(iCount,14)=temp_Now(iCount, 1); %% Bridge Number
                Total(iCount,15)=temp_Now(iCount, 3); %% Node Number
                Total(iCount,16)=kxscaler;
                Total(iCount,17)=xb0;
                %% NOTE: Total column 18 will be our 'complete' toggle
                %% Now keep track of our variances; which means keeping
                %% track of the square of the value for now.
                Total(iCount,19)=X(temp_Now(iCount, 3), 1)^2; %% Position of XB Node
                Total(iCount,20)=X(temp_Now(iCount,1), 1)^2; %% Position of Thin Fil Node
                if temp_Now(iCount, 2)==1
                    Total(iCount,21)=( (X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0 )^2; %% Total Bridge Strain
                    Total(iCount,22)=( (X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0 )^2; %% XB1 Bridge Strain
                else
                    Total(iCount,21)=( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% Total Bridge Strain
                    Total(iCount,23)=( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% XB2 Bridge Strain
                end
            end
        end
        temp_Prior=temp_Now; %% Update the now for our prior

        % Otherwise, we had briges before, and we need to account for them.
        % We need to update our accounting.
    else
        %% We have bound bridges
        %% Loop through the ones that are currently bound
        for iCount=1:Now_Bound_Length  %% Update our Run Matrix
            %% If the Now is consistent with the actin node
            %% found in the Prior, then grab the Big Total Index
            i_In=find(temp_Prior(:,1)==temp_Now(iCount, 1));
            %% If i_In exists in for being currently bound, and was not one
            %% of the original brigdges that started prior to SS_cutoff
            if (not(isempty(i_In))&&(Total(temp_Prior(i_In,4),2)>SS_cutoff))
                iRun=temp_Prior(i_In,4);
                temp_Now(iCount, 4)=Total(iRun, 1); %% Update in Temp_Now
                temp_Prior(i_In, 5)=1; %% Update finding this guy
                %% Now update our Total matrix
                Total(iRun,3)=Total(iRun,3)+1;   %%Duration in total
                Total(iRun,9)=Total(iRun,9)+X(temp_Now(iCount, 3), 1); %% Position of XB Node
                Total(iRun,10)=Total(iRun,10)+X(temp_Now(iCount,1), 1); %% Position of Thin Fil Node
                %% Also update our variances
                Total(iRun,19)=Total(iRun,19)+(X(temp_Now(iCount, 3), 1)^2); %% Position of XB Node
                Total(iRun,20)=Total(iRun,20)+(X(temp_Now(iCount,1), 1)^2); %% Position of Thin Fil Node
                if temp_Now(iCount, 2)==1
                    Total(iRun,4)=Total(iRun,4)+1; %% XB 1-State
                    Total(iRun,12)=Total(iRun,12)+((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0); %% XB1 Bridge Strain
                    Total(iRun,11)=Total(iRun,11)+((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0); %% Total Bridge Strain
                    %% Also update our variances
                    Total(iRun,22)=Total(iRun,22)+(((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0) )^2; %% XB1 Bridge Strain
                    Total(iRun,21)=Total(iRun,21)+(((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0) )^2; %% Total Bridge Strain
                    %% If we were previously XB2 and now XB1; update
                    %% counter XB2-->XB1
                    if temp_Prior(i_In, 2)==2
                        Total(iRun,7)=Total(iRun,7)+1; %% Step XB2-->XB1
                    end
                else
                    Total(iRun,5)=Total(iRun,5)+1; %% XB 2-State
                    Total(iRun,13)=Total(iRun,13)+( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) ); %% XB2 Bridge Strain
                    Total(iRun,11)=Total(iRun,11)+( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) ); %% Total Bridge Strain
                    %% Also update our variances
                    Total(iRun,23)=Total(iRun,23)+( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% XB2 Bridge Strain
                    Total(iRun,21)=Total(iRun,21)+( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% Total Bridge Strain
                    %% If we were previously XB1 and now XB2; update
                    %% counter XB1-->XB2
                    if temp_Prior(i_In, 2)==1
                        Total(iRun,6)=Total(iRun,6)+1; %% Step XB1-->XB2
                    end
                end
            else
                %% Else this guy is going to be a new entry and we
                %% must update the RunCount and grab a new
                [RunCount, N_wide]=size(Total);
                RunCount=RunCount+1; %% Yeah, we have a bridge--update
                %% Make a new row for the Total:
                Total=[Total; zeros(1, BinderCols)];
                iRun=RunCount;
                Total(iRun,1)=iRun;  %%Index of Run
                temp_Now(iCount, 4)=Total(iRun, 1); %% Update in Temp_Now
                Total(iRun,2)=current_step;  %% Starting Time-step of simulation
                Total(iRun,3)=1;   %%Duration in total
                Total(iRun,6)=0; %% Step XB1-->XB2
                Total(iRun,7)=0; %% Step XB2-->XB1
                Total(iRun,8)=0; %% Step XB2-->XB0
                Total(iRun,9)=X(temp_Now(iCount, 3), 1); %% Position of XB Node
                Total(iRun,10)=X(temp_Now(iCount,1), 1); %% Position of Thin Fil Node
                %% Also update our variances
                Total(iRun,19)=(X(temp_Now(iCount, 3), 1)^2); %% Position of XB Node
                Total(iRun,20)=(X(temp_Now(iCount,1), 1)^2); %% Position of Thin Fil Node
                if temp_Now(iCount, 2)==1
                    Total(iRun,4)=1; %% XB 1-State
                    Total(iRun,12)=(X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0; %% XB1 Bridge Strain
                    Total(iRun,11)=(X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0; %% Total Bridge Strain
                    %% Also update our variances
                    Total(iRun,22)=(((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0) )^2; %% XB1 Bridge Strain
                    Total(iRun,21)=(((X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1))-xb0) )^2; %% Total Bridge Strain

                else
                    Total(iRun,5)=1; %% XB 2-State
                    Total(iRun,13)=X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1); %% XB2 Bridge Strain
                    Total(iRun,11)=X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1); %% Total Bridge Strain
                    %% Also update our variances
                    Total(iRun,23)=( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% XB2 Bridge Strain
                    Total(iRun,21)=( X(temp_Now(iCount, 3), 1)-X(temp_Now(iCount, 1), 1) )^2; %% Total Bridge Strain

                end
                Total(iRun,14)=temp_Now(iCount, 1); %% Bridge Number
                Total(iRun,15)=temp_Now(iCount, 3); %% Node Number
                Total(iRun,16)=kxscaler;
                Total(iRun,17)=xb0;
            end
        end

        %% We have looped through our Now, but we need to make sure
        %% that we haven't missed anything for our bridges being off
        i_In=find(temp_Prior(:,5)==0);
        for iOff=1:length(i_In) %% Process our mean values and finialize
            %% Only update for our bridges that were later than the start
            %% time
            if (Total(temp_Prior(i_In,4),2)>SS_cutoff)
                iRun=temp_Prior(i_In(iOff,1),4);
                %% Now update our Total matrix (make the mean values)
                %% Total(iRun,3) is a count of our durations (or N for the
                %% mean)
                Total(iRun,9)=Total(iRun,9)/Total(iRun,3); %% Position of XB Node
                Total(iRun,10)=Total(iRun,10)/Total(iRun,3); %% Position of Thin Fil Node
                Total(iRun,12)=Total(iRun,12)/(Total(iRun,7)+1); %% XB1 Bridge Strain -- divide by the number of XB1 occurances
                Total(iRun,11)=Total(iRun,11)/Total(iRun,3); %% Total Bridge Strain
                Total(iRun,13)=Total(iRun,13)/Total(iRun,6); %% XB2 Bridge Strain  -- divide by the number of XB2 occurances
                %% Process our variances
                %% On first pass this will be <x^2>-<x>^2; p613 in Numerical
                %% Recipies in C.  However, I think this is only for large
                %% sample sizes.  We might best follow this up with a better
                %% algorithm in the future. Our mean was calculated just above.
                Total(iRun,19)=( (Total(iRun,19)^2)/Total(iRun,3) )-(Total(iRun,9)^2); %% Position of XB Node
                Total(iRun,20)=( (Total(iRun,20)^2)/Total(iRun,3) )-(Total(iRun,10)^2); %% Position of Thin Fil Node
                Total(iRun,22)=( (Total(iRun,22)^2)/(Total(iRun,7)+1) )-(Total(iRun,12)^2); %% XB1 Bridge Strain  -- account for XB1 occurances
                Total(iRun,21)=( (Total(iRun,21)^2)/Total(iRun,3) )-(Total(iRun,11)^2); %% Total Bridge Strain
                Total(iRun,23)=( (Total(iRun,23)^2)/Total(iRun,6) )-(Total(iRun,13)^2); %% XB2 Bridge Strain  -- account for XB2 occurances
                %% If we were previously XB2 and now we are off
                %% counter XB2-->XB0
                if temp_Prior(i_In(iOff,1), 2)==2
                    Total(iRun,8)=1; %% Step XB2-->XB0
                end
                %% Also update that we have completed this run
                Total(iRun,18)=1;
            end
        end
        temp_Prior=temp_Now; %% Update the now for our prior

    end

else % There are no bridges this time, but we need to keep track of if we ever had any

    %% If we had briges before, we need to do something with the
    %% details to account for no bridges this time
    if not(isempty(Total))

        %% We have looped through our Now, but we need to make sure
        %% that we haven't missed anything for our bridges being off
        i_In=find(temp_Prior(:,5)==0);
        for iOff=1:length(i_In) %% Process our mean values and finialize
            %% Only update for our bridges that were later than the start
            %% time
            if (Total(temp_Prior(i_In,4),2)>SS_cutoff)
                iRun=temp_Prior(i_In(iOff,1),4);
                %% Now update our Total matrix (make the mean values)
                %% Total(iRun,3) is a count of our durations (or N for the
                %% mean)
                Total(iRun,9)=Total(iRun,9)/Total(iRun,3); %% Position of XB Node
                Total(iRun,10)=Total(iRun,10)/Total(iRun,3); %% Position of Thin Fil Node
                Total(iRun,12)=Total(iRun,12)/Total(iRun,3); %% XB1 Bridge Strain
                Total(iRun,11)=Total(iRun,11)/Total(iRun,3); %% Total Bridge Strain
                Total(iRun,13)=Total(iRun,13)/Total(iRun,3); %% XB2 Bridge Strain
                %% Process our variances
                %% On first pass this will be <x^2>-<x>^2; p613 in Numerical
                %% Recipies in C.  However, I think this is only for large
                %% sample sizes.  We might best follow this up with a better
                %% algorithm in the future. Our mean was calculated just above.
                Total(iRun,19)=( (Total(iRun,19)^2)/Total(iRun,3) )-(Total(iRun,9)^2); %% Position of XB Node
                Total(iRun,20)=( (Total(iRun,20)^2)/Total(iRun,3) )-(Total(iRun,10)^2); %% Position of Thin Fil Node
                Total(iRun,22)=( (Total(iRun,22)^2)/Total(iRun,3) )-(Total(iRun,12)^2); %% XB1 Bridge Strain
                Total(iRun,21)=( (Total(iRun,21)^2)/Total(iRun,3) )-(Total(iRun,11)^2); %% Total Bridge Strain
                Total(iRun,23)=( (Total(iRun,23)^2)/Total(iRun,3) )-(Total(iRun,13)^2); %% XB2 Bridge Strain
                %% If we were previously XB2 and now we are off
                %% counter XB2-->XB0
                if temp_Prior(i_In(iOff,1), 2)==2
                    Total(iRun,8)=1; %% Step XB2-->XB0
                end
                %% Also update that we have completed this run
                Total(iRun,18)=1;
            end
        end
        temp_Prior=temp_Now; %% Update the now for our prior
    end

    %% Otherwise, it is just fine have no bridges now, and
    %% we never had bridges before, so just do nothing
    temp_Prior=temp_Now; %% Update the now for our prior
end

