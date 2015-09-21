%% BCWT, May 29, 2009
%% Write a few output files to record the code that will run.

function [OutDir] = Initialize_v1(CurDir, Set)

%% OutDirectory for Data
%% GNU Octave has a bug - mkdir doesn't recurse, unlike matlab's

DataDir = sprintf('%s%sDataFiles', CurDir, filesep);
if(~isdir(DataDir))
    mkdir(DataDir);
end

OutDir = sprintf('%s%sSet%d%s', DataDir, filesep, Set, filesep);
if(~isdir(OutDir))
    mkdir(OutDir);
end

%% OutDirectory for Code
CodeDir = sprintf('%s%sCodeBak', DataDir, filesep);
if(~isdir(CodeDir))
    mkdir(CodeDir);
end
  
%% Write all the m files from the Current Directory directory over
%% to the Code Directory for the simulation
  
copyfile([CurDir filesep '*.m'], CodeDir, 'f');
