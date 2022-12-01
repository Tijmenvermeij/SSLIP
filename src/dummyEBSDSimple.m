function ebsd = dummyEBSDSimple(ori,X,Y)
    % creates a "Dummy" EBSD variable from an orientation and some X and Y
    % data
    
    
    prop.x = X(:);
    prop.y = Y(:);
    
    CS = ori.CS;
    oris = repmat(ori,size(prop.x));
    
    phase = ones(size(prop.x));
    ebsd = EBSD(oris,phase,CS,prop);
    
    % define unitcell, only for square grid
    spacing = mean(diff(X(1,:)));
    ebsd.unitCell = [spacing/2 -spacing/2
        spacing/2 spacing/2
        -spacing/2 spacing/2
        -spacing/2 -spacing/2];
    
    ebsd = ebsd.gridify;
    ebsd.scanUnit = 'um';

end