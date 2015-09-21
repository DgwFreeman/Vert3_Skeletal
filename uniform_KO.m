function [mat] = uniform_KO(s1,s2,frac)

mat=ones(s1, s2); %% Fill as though no KOs
%[~, jTest]=find(mat==1); %% Gather row and col indicies to
%% filter on the KO
%for iii=1:s2 %% Loop down each row, looking at first half
    %% then the second half to fill with zero if s2<(1-mat)
    %% Now the matrix is zeros, update with the Tn that will actually
    %% exist.
%    if(iii<=(s2/2)) %% First Half
for iii=1:s2/2
        if(2*iii/s2<(1-frac) )
            mat(:, iii)=0;
        end
end
for iii=floor(s2/2)+1:s2
%    if(iii>(s2/2)) %% 2nd Half
        if(2*iii/s2-1 < (1-frac))%((iii-s2/2)/(s2/2))<(1-frac) )
            mat(:, iii)=0;
        end
%    end
end