function [ A ] = makemovie( data, numframes, clims)
%function [ A ] = makemovie( data, numframes, clims)makes a matlabl movie file given an array of data.
%   The input arguments are the array of data to be movified, the number of
%   frames in the file and the colour limits. The output is a single data
%   object containing the movie. To save the movie as an avi use movie2avi.


fig1 = figure(1);

winsize = get(fig1,'Position');
winsize(1:2) = [0,0];

A=moviein(numframes,fig1,winsize);
for i=1:numframes
   image(:,:) = (data(i,:,:));
   imagesc(image,clims)
   axis xy
   colorbar
  
   A(:,i)=getframe(fig1,winsize);
 end

end

