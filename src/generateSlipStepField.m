function [U] = generateSlipStepField(ebsd,sSLocal,traceLocation,Magnitude)
%generateSlipStepField: Generate a slip step displacement field, based on
%one slip system and a magnitude, slip step is in 2 directions, i.e. on
%both sides of the trace
%   Input required ebsd data, on which the step will be applied. The step
%   will be applied on the full input area

% Generate a slip step displacement field, based on
% one slip system and a magnitude
%
% Input
%  ebsd   - @EBSD (on which slip system is applied, full area is used
%  sSLocal - @SlipSystem (local slip system (already rotated to spec), used
%    to compute slip step
%  traceLocation - double [XPos Ypos]: a position which should be on the
%  slip trace, in same unit and configuration as ebsd pos data
%  Magnitude - double: slip magnitude to be applied, in same unit as ebsd
%  pos data
%
% Output
%  U  - @vector3d field with displacements
%
% Authors: Tijmen Vermeij 10-04-2020
% Eindhoven University of Technology, Hoefnagels Group
% Based on work of Tim Ramirez MSc thesis

% check to see if there is only one slipsystem in sSLocal
if numel(sSLocal) > 1
    error('Only one slipSystem allowed in sSLocal in this function')
end

%% Determine slip traces of local slip systems

nFrontplane = vector3d.Z;

% get projected trace and b 
[trace_dSys, b_front, bSys] = getTraceProjB(sSLocal,nFrontplane);


% initialize displacement field
U = vector3d(zeros(size(ebsd.x)),zeros(size(ebsd.x)),zeros(size(ebsd.x)));

% convert position field into vector field
pos0 = vector3d(ebsd.x,ebsd.y,zeros(size(ebsd.x)));
% conv traceposition to vector
traceLocVec = vector3d(traceLocation(1),traceLocation(2),0);
% convert posvec field to have tracepos as origin
pos1 = pos0 - traceLocVec;

% perpendicular vector of trace
PerpTrace_dSys = cross(trace_dSys,nFrontplane);


PerpTrace_dSys = sign(dot(PerpTrace_dSys,xvector)) * PerpTrace_dSys;

% calc distance from trace line
dis = dot(PerpTrace_dSys,pos1);

% convert distance to pixels
disPix = dis ./ ebsd.dx;

% create matrix which will be used for the final slip step field (will go
% from -1 to 1, and will be multiplied by the slip direction...)
slipFraction = zeros(size(disPix));

% if distance from trace is larger then half a pixel size, make it 0.5,
% conversely -0.5
slipFraction(disPix > 0.5) = 0.5;
slipFraction(disPix < -0.5) = -0.5;

slipFraction(disPix < 0.5 & disPix > -0.5) = disPix(disPix < 0.5 & disPix > -0.5);

U = slipFraction * bSys * Magnitude;


