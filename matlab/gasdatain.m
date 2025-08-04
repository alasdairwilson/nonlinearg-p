function [ X ] = gasdatain( file,num,rows,columns )
%Load in files from the lax wendroff simulation [X] = gasdatain( file,num,offset)
%  Output is an array containing the output from the simulation, input is
%  the file from which you wish to load data (eg. agro1.out), num the
%  number of timesteps you wish to load eg for a 1000 sec program with a
%  writeout every 40 timesteps num would be 25 for the complete run. O
row = 0;
X =  zeros (num,rows,columns);
for i=1:num
    start = 1+(i-1)*(rows+3);
    range = [start,0,start+rows-1,columns-1];
    
    X(i,:,:) = dlmread(file,',',range);
    
end
return
