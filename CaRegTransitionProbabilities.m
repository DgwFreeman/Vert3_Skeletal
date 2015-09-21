function [transitions] = CaRegTransitionProbabilities(BaseTransitions, CoopPass, neighbors, TFstate)

%% Returns transition probabilities for calcium regulation.
% nodes have three states and this is where they transition between them.
% these probabilities are organized in a 3x2 matrix. conceptually,
% [p00 p01 p02; p10 p11 p12; p20 p21 p23] where pxy is the probability
% of transitioning from x to y. but we cut out p00 p11 p22 to get 3x2.

%% unpack parameters
TF3TF12=CoopPass(1, 1);
XB2TF12=CoopPass(2, 1);
XB3TF12=CoopPass(3, 1);
XB2TF23=CoopPass(4, 1);
XB3TF23=CoopPass(5, 1);
TF3TF23=CoopPass(6, 1);
XB2TF21=CoopPass(7, 1);
XB3TF21=CoopPass(8, 1);
TF3TF21=CoopPass(9, 1);
XB2TF32=CoopPass(10, 1);
XB3TF32=CoopPass(11, 1);
TF3TF32=CoopPass(12, 1);

%% before cooperativity
transitions = BaseTransitions;

%% actual cooperativity changes
%% if there's a neighbor node in state TF3
if max(TFstate(neighbors, 1)==2)>0
    transitions(1,2) = transitions(1,2)*TF3TF12;
    transitions(2,3) = transitions(2,3)*TF3TF23;
    transitions(2,1) = transitions(2,1)*TF3TF21;
    transitions(3,2) = transitions(3,2)*TF3TF32;
end
%% if there's a neighbor node with XB3
if max(TFstate(neighbors, 2)==2)>0
    transitions(1,2) = transitions(1,2)*XB3TF12;
    transitions(2,3) = transitions(2,3)*XB3TF23;
    transitions(2,1) = transitions(2,1)*XB3TF21;
    transitions(3,2) = transitions(3,2)*XB3TF32;
end
%% if there's a neighbor node with XB2
if max(TFstate(neighbors, 2)==1)>0
    transitions(1,2) = transitions(1,2)*XB2TF12;
    transitions(2,3) = transitions(2,3)*XB2TF23;
    transitions(2,1) = transitions(2,1)*XB2TF21;
    transitions(3,2) = transitions(3,2)*XB2TF32;
end

transitions(1,1) = 1 - transitions(1,2) - transitions(1,3);
transitions(2,2) = 1 - transitions(2,1) - transitions(2,3);
transitions(3,3) = 1 - transitions(3,1) - transitions(3,2);