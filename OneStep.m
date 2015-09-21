%% AMW
%% A single simulation step

%% Note for AMW: Matlab copies arguments, so sideffecting on copies is ok but probably inefficient

%% Now that dependencies are mostly taken care of it's time to think about optimization.
% occupied changes size every iteration (occupied(end+1) adds a point).
% LAUREL is totally constant - can we kill that branch?
% To be specific, I think LAUREL has zero elements only on certain rows,
% rows that are a multiple of NumBridges (eg LAUREL(61,:) == [0,0,0] rn)
% Reasonably sure it has NMYO such rows.
% XBKO is also semi-constant, only different between runs. Obviously it's
% not quite as patterned as LAUREL either - though with uniform KO...
% Definitely can simplify the inner loop here! Right after the isempty(find
% it computes transition probabilities, but branches on the start state to
% do. That's probably unnecessary. Actually I'll change a bit rn.

function [CountOff, TFstate, ATPcycle] = OneStep(dt, HARDY, LAUREL, X, TFstate, N_Thick_Start, MYONodes, ACTNodes, XBKO, thermochem)

ATPcycle=0;

CountOff=0; % Reverse transitions this step

%%%*****************************************************************
%%% Now that we are RU state, check about the XB Binding

% occupied=0; %this will set occupied to 0, and a 1 by 1 matrix.
occupied = zeros(ACTNodes);

for iXB = 1:MYONodes
    for iStart=1:N_Thick_Start
        % loop over every myosin node--looking across the N_Thick_Start
        % Because LAUREL is no-longer a vector, we will have to shuffel
        % this to work with an Array--which will work out just fine now
        if LAUREL(iXB, iStart) % 0 implies handle node
            % or no XB pointing to a thick filament so skip those.

            % LAUREL(iXB, iStart) is a number between 1-24--pointing to
            % the row of HARDY that could be co-linearly aligned
            % The XB which is is facing has been stored in the
            % TFstate col. 3 as a number iXB+ACTNodes.


            %% Check whether this XB is knocked out or now
            %% If working = 1; if Knocked Out = 0
            if XBKO(iXB, iStart)==1

                % Where the variable 'node' for the rest of this
                % loop refers to the Actin Node for a given iXB
                % might be looking at.

                % Check here to see whether our XB head is indeed bound
                % Pull back the node (thin filament node we are looking for)
                node=find(TFstate((HARDY(LAUREL(iXB,iStart),:)), 3)==(iXB+ACTNodes),1);
                % This will be an index to use in HARDY, returning
                % the act node to look at.
                if isempty(node)
                    % If it has returned an empty set, then we will need to
                    % search over the binding sites and find the nearest
                    % node for it to bind.

                    % So here--we take the difference between the XB position and all the
                    % potential co-linear sites for our thin filament target zones
                    % where dist is the value of the minimum difference
                    % where node is the index used just below--a bit unclean but oh well
                    [~, node] = min(abs(X(iXB+ACTNodes) - X(HARDY(LAUREL(iXB,iStart),:))));
                    % Set the node we are looking as row in TFstate
                    % here node is the thin filament node
                    node=HARDY(LAUREL(iXB,iStart), node);
                    % Set the distance as a vector, rather than abs
                    % value as above.
                    dist=X(iXB+ACTNodes,1) - X(node,1);
                else
                    % it will already be bound somewhere.  Now we
                    % just need to set the node we are looking at
                    % in the row of TFstate
                    node=HARDY(LAUREL(iXB,iStart), node);
                    dist=X(iXB+ACTNodes,1) - X(node,1);
                end


                % Now we are all set with our node we are looking at and we
                % have calulated the distance, then we
                % can proceed with looking at the binding for a given iXB.
                % The XB state it kept in the TFstate col. 2, and the
                % XB it is bound with is kept in the TFstate col. 3.
                % Case 0 --Unbound XB
                % Case 1 --Stongly Bound XB, with low force into
                % the matrix, such that the energy from ATP
                % hydrolysis remains in the XB head.
                % Case 2 --Strongly Bound XB, with high force into
                % the latice, such that the E from ATP hydrolys is
                % transfered from the myosin molecule into the
                % lattice and there is a displacement.
                %
                % The XB dynamics are set up and defined in XBRates3
                % All Cases--Check for the Regulation TFstate col. 1 to
                % be in state 2.

                % If the node has already been attached, then another XB cannot look at
                % it as well.  So, check this by looking at whether occupied is
                % empty?
                % if ( (isempty(find(occupied==node, 1))) ) %& ( (dist > -12)&(dist<12) ) );
                if (occupied(node) == 0)
                    % The node has not been occupied if we get this
                    % far.  However, we should keep track that is
                    % the node gets occupied, we will fill the free
                    % end of the occupied vector using: occupied(end+1,1)=node;
                    
                    % clear p01 p02 p12 p10 p20 p21 %PERFMOD
                    % Compute forward and reverse transition probabilities.
                    if (TFstate(node,1) == 2) % Node is available.
                        [k01,k10,k12,k21,k20,k02] = XBRates3(dist, thermochem);
                        %  [k01,k10,k12,k21,k20,k02] = mexDooItRates(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD);
                        gx = rand;
                        switch TFstate(node, 2)
                            case 0 %%%%%%%% UNBOUND
                                p01 = dt*k01;
                                p02 = dt*k02;
                                if (gx <= p01) % It goes forward to weak state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    TFstate(node, 2)=1;
                                    TFstate(node, 3)=(iXB+ACTNodes); %% Load thick filament node bound at target zone
                                    occupied(node) = 1;  % The node has not been occupied
                                elseif (gx >= 1.0 - p02) % Reverse rarely happens.
                                    TFstate(node, 2)=2;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    CountOff=CountOff+1;
                                    occupied(node) = 1;  % The node has not been occupied
                                    
                                end
                            case 1  %%%%%%%% WEAK BOUND
                                p12 = dt*k12;
                                p10 = dt*k10;
                                if (gx <= p12) % It goes forward to strong state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    TFstate(node, 2)=2;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    occupied(node) = 1;  % The node has not been occupied
                                elseif (gx >= 1.0 - p10) % Reverse rarely happens.
                                    TFstate(node, 2)=0;
                                    TFstate(node, 3)=0;
                                end
                            case 2  %%%%%%%% STRONG BOUND
                                p20 = dt*k20;
                                p21 = dt*k21;
                                if (gx <= p20) % It goes forward to unbound state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    %CountCycle=CountCycle+1;
                                    ATPcycle=ATPcycle+1;
                                    TFstate(node, 2)=0;
                                    TFstate(node, 3)=0;
                                elseif (gx >= 1.0 - p21) % Reverse rarely happens.
                                    TFstate(node, 2)=1;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    occupied(node) = 1;  % The node has not been occupied
                                end
                        end
                    end
                end
            end % end loop on whether XB is knocked out
            
        end % if over the 'whether' handel? (If over Hardy(iXB==0) to over myosin handles.
    end    % end loop over iStart
end    % end loop over Myo. Heads