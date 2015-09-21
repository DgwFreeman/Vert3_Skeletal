function LAUREL = MakeLAUREL(Angle_Crowns, NumBridges, Angle_Thick_Start, MYONodes, N_Thick_Start, NMYO)

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