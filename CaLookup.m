%% BCWT; August 18, 2010
%% Gather our Thick facing Thin information.  Final outputs that we care
%% about will be still labeled:
%% LAUREL;  However, LAUREL now has size (NumBridges, N_Thick_Start), or
%% 244 X 4 for this new fly arrangement.
%% HARDY; Still has the target sites for the bridges, but HARDY now is a
%% matrix size (ThinRow_Index, NumSites/3) or 24 by 30 for the current fly.
%% Tn;--our Tn matrix didn't change (not just yet anyhow)

%% HARDY will be binding sites or actin node numbers:
%HARDY=zeros(ThinRow_Index, floor(NumSites/3) );
HARDY=[
% For actin node a)
2:3:91;
3:3:91;
4:3:91;
% For actin node b)
93:3:182;
94:3:182;
95:3:182;
% For actin node c)
184:3:273;
185:3:273;
186:3:273;
% For actin node d)
275:3:364;
276:3:364;
 277:3:364;
% For actin node e)
366:3:455;
367:3:455;
368:3:455;
% For actin node f)
457:3:546;
458:3:546;
459:3:546;
% For actin node g)
548:3:637;
549:3:637;
550:3:637;
% For actin node h)
639:3:728;
640:3:728;
641:3:728];

%% HOWERVER--we now need to keep track of all the angular assignments for
%% these rows of HARDY--making up the 8 thin filaments.
%% For our thin filaments (8 of them--we have the alpha, beta, gamma on
%% each).  This will be labled with a to h, with each row being the
%% direction that alpha, beta, gamma points.  This is the direction the
%% binding site will point--so really opposite (180) from the pointing of
%% the co-linear cross-bridge head direction.
a_Points=mod((120+[0;120;240]), 360); %[120; 240; 0];
b_Points=mod((300+[0;120;240]), 360); %[120; 60; 180];
c_Points=mod((120+[0;120;240]), 360);
d_Points=mod((300+[0;120;240]), 360);
e_Points=mod((120+[0;120;240]), 360);
f_Points=mod((300+[0;120;240]), 360);
g_Points=mod((120+[0;120;240]), 360);
h_Points=mod((300+[0;120;240]), 360);

%% Now we can concat our thin filament directions, where the row will
%% represent the targed row for the thick filaments to adress; having
%% sites for our new LAUREL and HARDY.
Thin_Points=[a_Points; b_Points; c_Points; d_Points; ...
    e_Points; f_Points; g_Points; h_Points];

%% Now assign each thick filament options, for co-linear facing thin
%% filament rows--index of the row of HARDY:
A_faces=[1; 5; 9; 10; 14; 18];
B_faces=[21;23;7;12;2;4];
C_faces=[3;11;13;24;20;16];
D_faces=[15;17;19;6;8;22];
%Thick_facing_Thins=[A_faces; B_faces; C_faces; D_faces];

%% Note the thick filament in the end--row 61=0;  This is the M-Line nodes.
Thick_Points=Angle_Crowns*(0:1:(NumBridges-1))'; %% The angle our crown will point--just the Start#1 to begin
Thick_Points=mod([Thick_Points, Thick_Points+Angle_Thick_Start, Thick_Points+(2*Angle_Thick_Start)], 360); %Now the rest
LAUREL=zeros(MYONodes, N_Thick_Start); %Thick_Thin pairs--where i,j represents node#,Start#; content is the target row in HARDY

[iRow, jCol]=size(Thick_Points);
%% Now loop through and get to where the cross-bridge will point, given
%% cutoff angle 20 degrees
Cutoff_Angle=20;
for i_Thick=1:NMYO %% Loop over each thick filament
    switch i_Thick
        case 1
            Thin_Number=A_faces;
            Thin_Angles=Thin_Points(A_faces, 1);
        case 2
            Thin_Number=B_faces;
            Thin_Angles=Thin_Points(B_faces, 1);
        case 3
            Thin_Number=C_faces;
            Thin_Angles=Thin_Points(C_faces, 1);
        case 4
            Thin_Number=D_faces;
            Thin_Angles=Thin_Points(D_faces, 1);
    end
    
    for i=1:iRow-1 %% Each thick filament will have same points, but we will have to loop on NMYO
        for j=1:jCol %% As well as looping on each Start# within a filament
            
            %% Compare the thick filament angle to nearest thin filament angle.
            %% sort the value on the index; example absolute value = 10 degrees
            %% at the third index.
            [val,index]=min(abs(Thin_Angles-Thick_Points(i,j)));
            if val <= Cutoff_Angle
                %LAUREL(i, j)=Thin_Number(index,1); %% Test for first thick
                LAUREL(((i_Thick-1)*NumBridges)+i, j)=Thin_Number(index,1); %% Something here to loop over the NumBridges index and multiple filaments
            end
        end
    end
end

%% The Tn is built same as the old HARDY/LAUREL COMBO

% Build the equivalent of HARDY for the Troponin Nodes.
% This becomes a matrix with NACT rows and NTn col.  Where
% each half row (1:NTn/2) and NTn/2:NTn represent the halves
% of the rotating helix of actin.  This is currently built such
% that the Handel Nodes (1, 92, ...) are kept out.  The first half
% ranges 2:4:91, then 3:4:91 to fill (as would be for the first line
% and the first Actin filament).  What we are assuming here is that
% each Tn will control the actin it corresponds to (being a super
% actin) and two actins down the way.  This is then staggered with Tn
% being sprikled in alternating pairs.  Such as Tn is at index 2, and 3, 
% for helix 1 and 2 respectively, controlling nodes 4 and 5, 
% respectively.  This will then be the span of the Tm.  From
% one Tn node of the helix until the next.  The way this works
% out the Tm then spans about 37 nm.
% This then repeats to follow at Tn being 6 and 7 ...
% and does so starting at the handel for each actin filament.

% There is a little bit of wierd at the free end, because the two Tn
% at the free end only control themselves (90 and 91) because there
% are none further down the way.  So this is a good start.

NTn=46; % Number of Troponin per filament
for i=1:NACT

    Tn(i, 1:46)=[(((i-1)*NumSites)+2):4:(i*NumSites), (((i-1)*NumSites)+3):4:(i*NumSites)];

end
% So easy as that we have a matrix of where the Tn are sprinkled on the 
% rows of actin filaments, both devided into the helix they will cover.


clear a_Points b_Points c_Points d_Points e_Points f_Points g_Points h_Points ...
    Thin_Points Thick_Points iRow jCol

% leaving LAUREL, HARDY, MIRROR_ON, and set_of_if_nodes

% A_faces
% B_faces
% C_faces
% D_faces
% 
% 
% 
% Thick_Points
% Thick_Points
% Thick_Thin=zeros(244,4); %% Now LAUREL
