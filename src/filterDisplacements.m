function [ data ] = filterDisplacements(U,V,options)
% filter displacement data to reproduce strain filter smoothing of
% commercial DIC packages

% created by Tijmen Vermeij & Johan Hoefnagels

%%% create filter

if ~isfield(options,'filter')
    options.filter = true;
end
if ~isfield(options,'cutofffraction')
    options.cutofffraction = 1;
end
if ~isfield(options,'filt_std')
    options.filt_std = 3;
    warning('no filter standard deviation given, using 3 pixels');
end

% define filter
filt_std = options.filt_std;
EffectiveWindowSize = 2*sqrt(-2*log(1-options.cutofffraction)*filt_std^2);
if options.cutofffraction < 1
    window = ceil(EffectiveWindowSize)+1-rem(ceil(EffectiveWindowSize),2); 
    imageFilter=fspecial('gaussian',window ,filt_std);
    imageFilter(imageFilter<(1-options.cutofffraction)*max(imageFilter(:)))=0;
else
%         window = 15*ceil(EffectiveWindowSize)+1-rem(ceil(EffectiveWindowSize),2); 
    window = ceil(15*filt_std) + 1 - rem(ceil(15*filt_std),2);
    imageFilter=fspecial('gaussian',window ,filt_std);
end


%%% perform filtering
% find nans
nanU = isnan(U);
nanV = isnan(V);

% blur displacement fields (use nanconv script which can handle nan values well)
Uf = nanconv(U,imageFilter, 'nonanout','edge');
Vf = nanconv(V,imageFilter, 'nonanout','edge');

% replace NaNs
Uf(nanU) = NaN;
Vf(nanV) = NaN;

% store filter
data.imageFilter = imageFilter;
data.U = Uf;
data.V = Vf;



end

