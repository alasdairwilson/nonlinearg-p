clear density dataname
dataname(1,:) = 'pldeb.out';
dataname(2,:) = 'plbyb.out';
dataname(2,:) = 'plbxb.out';
% dataname(3,:) = 'agvx1.out';
% dataname(4,:) = 'agvz1.out';
numframes=25;
parfor i = 1:size(dataname,1);
density(i,:,:,:) = gasdatain(dataname(i,:),numframes,203);
end

