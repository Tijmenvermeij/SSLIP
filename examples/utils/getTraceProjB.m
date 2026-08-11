function [traces_d, b_front, b_front_norm,traces] = getTraceProjB(sSLocal,nFrontplane)
% calculate traces, traces with direction pointing down (positive y),
% projected burgers vector, and normalized projected burgers vector

% based on local SlipSystem and frontplane normal

% Authors: Tijmen Vermeij 10-04-2020
% Eindhoven University of Technology, Hoefnagels Group
% Based on work of Tim Ramirez MSc thesis

% ensure b and n are normalized
sSLocal.b = sSLocal.b.normalize;
sSLocal.n = sSLocal.n.normalize;

% calc traces on top surface
traces = trace(sSLocal,nFrontplane);

% calc traces as directional vectors, pointing down in positive y direction
traces_d = traces;
traces_d.antipodal = 0;
traces_d = traces_d .* sign(traces_d.y);

% calc burgers vector projected in the plane, length is changed from 1 to account for projection 
b_front = cross(nFrontplane,cross(sSLocal.b,nFrontplane));
b_front = dot(b_front,sSLocal.b).*b_front;
b_front_norm = b_front.normalize;


end

