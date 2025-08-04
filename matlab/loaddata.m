function [ datain ] = loaddata(numframes, dim, dataname)


parfor i = 1:size(dataname,1);
datain(i,:,:,:) = gasdatain(dataname(i,:),numframes,dim(1),dim(2));
end

end

