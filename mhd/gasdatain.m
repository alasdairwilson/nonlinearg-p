function [ X ] = gasdatain( file,num,offset )
%Load in files from the lax wendroff simulation [X] = gasdatain( file,num,offset)
%  Output is an array containing the output from the simulation, input is
%  the file from which you wish to load data (eg. agro1.out), num the
%  number of timesteps you wish to load eg for a 1000 sec program with a
%  writeout every 40 timesteps num would be 25 for the complete run. Offset
%  is in order to fix inconsistancies in the program, it is the number of
%  lines in the matrix+number of lines between timesteps in the datafile.
%  For most files this should be 203 but try 205/204 if returns error.
clear X;
row = 0;
X =  zeros (num,199,199);
for i=1:num
    row = 2+offset*(i-1);

    range = [row,1,row+198,199];

    X(i,:,:) = dlmread(file,'',range);
    
end
return


