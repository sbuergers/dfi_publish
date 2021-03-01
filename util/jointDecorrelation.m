function [ dataJD,data_3D ] = jointDecorrelation(data,ERP)
% Joint decorrelation algoritms
%
% (c) Giulio Degano
%
% De Cheveigné, Alain & Parra, Lucas C., "Joint decorrelation, a versatile tool for multichannel data analysis", Neuroimage, 2014

num_time_points=length(data.time{1,1});
data_concat=zeros(num_time_points*length(data.trial),length(data.label));


disp('Concatenating data...')
counter=1;
for i=1:length(data.trial)
    data_concat(counter:counter+num_time_points-1,:)=data.trial{i}';
    counter=counter+num_time_points;
end
disp('done!')

 

%% JD

 

% Covariance C0
C0=(data_concat'*data_concat);
[P,d] = eig(C0);
D = diag(1./sqrt(diag(d)+1e-18));
W = (P*D)*P';
Z=data_concat*W;



Z_filt=ERP.avg';

 

% Covariance C1
C1=(Z_filt'*Z_filt);
[Q,~] = eig(C1);


data_clean=data_concat*W*Q;

 

%% Rebuilding data struct


disp('Rebuilding data...')
counter=1;
dataJD=data;
dataJD.trial={};
data_3D=zeros(num_time_points,length(data.label),length(data.trial));
for i=1:length(data.trial)
    dataJD.trial{i}=data_clean(counter:counter+num_time_points-1,:)';
    data_3D(:,:,i)=data_clean(counter:counter+num_time_points-1,:);
    counter=counter+num_time_points;
end
disp('done!')

 

end


% // eof