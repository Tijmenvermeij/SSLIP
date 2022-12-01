function crs = coarsegrainDisp(f,g,x,y,coarsegrainstep)
% this function coarse grains the data in f and g
% originally used to coarsegrain images which were meant for DIC

% code written by Jan Neggers

% Ncoarse = options.Ncoarse;
[n m]   = size(f);

% crs graining factor
crs.factor = coarsegrainstep - 1;
% superpixelsize
crs.sps = 2^(crs.factor);


% Crop to integer superpixel size
% =================================

nmod = mod(n,crs.sps);
In = floor(1+nmod/2):(n-(nmod-floor(nmod/2)));

mmod = mod(m,crs.sps);
Im = floor(1+mmod/2):(m-(mmod-floor(mmod/2)));

% crop
f = f(In,Im);
g = g(In,Im);
y = y(In);
x = x(Im);

% new size of the image
[n m] = size(f);

% number of superpixels
if coarsegrainstep ~= 1
    crs.Npx = [n m]./crs.sps ;
else
    crs.Npx = [n m] ;
end



% ============================
% Coarse Graining
% ============================

if coarsegrainstep ~= 1
    
        % left coarse graining permutation matrix
        % (unity matrix with each row duplicated superpixel size times)
        Pn = repmat(eye(crs.Npx(1)),crs.sps,1);
        Pn = reshape(Pn,crs.Npx(1),n)';

        % right coarse graining permutation matrix
        Pm = repmat(eye(crs.Npx(2)),crs.sps,1);
        Pm = reshape(Pm,crs.Npx(2),m)';
        
%         % crs grain the mask
%         crs.m = ( Pn' * mask * Pm ) ./ crs.sps^2;

        % crs grain x and y
        crs.x = ( x * Pm ) ./ crs.sps;
        crs.y = ( y * Pn ) ./ crs.sps;
        
    if ~any(isnan(f(:)) | isnan(g(:)))
        % Coarse grain the image (fast if there are no NaN's)
        
        % crs grain f and g
        crs.f = ( Pn' * f * Pm ) ./ crs.sps^2;
        crs.g = ( Pn' * g * Pm ) ./ crs.sps^2;
    else
        % Coarse grain (slower but can handle NaN's)

        % initialize the image matrices
        crs.f = zeros(crs.Npx(1),crs.Npx(2));
        crs.g = zeros(crs.Npx(1),crs.Npx(2));

        % initialize matrices which count NaN's
        fcnt  = zeros(crs.Npx(1),crs.Npx(2));
        gcnt  = zeros(crs.Npx(1),crs.Npx(2));

        % forloop over all subpixels in a superpixel
        for in = 1:crs.sps
            for im = 1:crs.sps

                % get the same subpixel in each superpixel
                In = in:crs.sps:n;
                Im = im:crs.sps:m;
                ftmp = f(In,Im);
                gtmp = g(In,Im);

                % add that pixel to the previous pixel (only if it is not NaN)
                crs.f(~isnan(ftmp)) = crs.f(~isnan(ftmp)) + ftmp(~isnan(ftmp));
                crs.g(~isnan(gtmp)) = crs.g(~isnan(gtmp)) + gtmp(~isnan(gtmp));

                % count the number of nan's encountered per superpixel
                fcnt  = fcnt + isnan(f(In,Im));
                gcnt  = gcnt + isnan(g(In,Im));
            end
        end

        % devide by the number of pixels (excluding NaN's)
        crs.f = crs.f ./ (crs.sps^2 - fcnt);
        crs.g = crs.g ./ (crs.sps^2 - gcnt);
        
        % superpixels which consited of all NaN pixels should remain NaN
        crs.f( fcnt == crs.sps^2 ) = NaN;
        crs.g( gcnt == crs.sps^2 ) = NaN;
        
        % correction for divide by zero
        crs.f(isinf(crs.f)) = NaN;
        crs.g(isinf(crs.g)) = NaN;
    end

else % highest detail
    crs.f = f;
    crs.g = g;
    crs.x = x;
    crs.y = y;
%     crs.m = mask ;
end



% get the size of a superpixel (in mu)
pixelsize(1) = mean(diff(x));
pixelsize(2) = mean(diff(y));
crs.pixelsize = pixelsize*crs.sps;



% the number of pixels in the ROI
crs.Npx = size(crs.f);

% create position fields
[crs.X crs.Y] = meshgrid(crs.x,crs.y);

% if options.debug
%     assignin('base','crs',crs)
% end

% =======================================================================
end
