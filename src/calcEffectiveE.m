function [Eeff] = calcEffectiveE(Hxx,Hxy,Hyx,Hyy)
% calcEffectiveE calc effective shear strains from displacement gradient tensor
% components
% Author: Tijmen Vermeij
      Eeff = sqrt(0.5*(Hxx-Hyy).^2 + 0.5*(Hxy+Hyx).^2);
%     
end

