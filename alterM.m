function M = alterM(M, TFstate, kxscaler)

%%%*****************************************************************

% this block modifies the matrix which encodes the connections
% between myosin and actin nodes along XB. (M)
% For any force bearing XB/ACTIN pair (state 2==strong bound high force, and state 1
% == strongly bound, with low force).
% The matrix will be loaded.S
temp=find(TFstate(:,2)>0);

for i=1:length(temp)
    %         if i==1
    %             disp([temp, TFstate(temp, 2), TFstate(temp, 3), BigBoundState(temp, (iStep)), BigBoundPair(temp, (iStep))])
    %         end
    % NOTE on what we have.
    %temp_site=temp(i);
    %temp_brdg=(TFstate(temp(i),3));
    % Fill the actin diagonal
    M(temp(i),temp(i)) = M(temp(i),temp(i)) + kxscaler;
    % Fill the myosin diagonal
    M(TFstate(temp(i),3),TFstate(temp(i),3)) = M(TFstate(temp(i),3),TFstate(temp(i),3)) + kxscaler;
    % Fill the cross term on the actin row
    M(temp(i),TFstate(temp(i),3)) = M(temp(i),TFstate(temp(i),3)) - kxscaler;
    % Fill the cross term on the myosin row
    M(TFstate(temp(i),3),temp(i)) = M(TFstate(temp(i),3),temp(i)) - kxscaler;
end