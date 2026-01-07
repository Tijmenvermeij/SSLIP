function ebsd = dummyEBSDSimple(ori,X,Y)
    % creates a "Dummy" EBSD variable from an orientation and some X and Y
    % data
    
    
    prop.x = X(:);
    prop.y = Y(:);
    
    pos = vector3d(prop.x,prop.y,zeros(size(prop.x)));

    CS = ori.CS;
    oris = repmat(ori,size(prop.x));
    
    phase = ones(size(prop.x));
    % ebsd = EBSD(oris,phase,CS,prop);
    ebsd = EBSD(pos,oris,phase,CS, struct());
    
    % gridify EBSD data
    ebsd = ebsd.gridify;
    ebsd.scanUnit = 'um';

end