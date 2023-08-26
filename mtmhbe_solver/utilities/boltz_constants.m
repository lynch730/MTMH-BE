function c = boltz_constants
    c.KB = 1.380649e-23;       % Boltzmann constant J/K, m^2*kg*s^-2*K^-1
    c.QE = 1.602176634e-19;    % Electron charge, J/eV, C
    c.ME = 9.1093837015e-31;   % Mass of the electron, kg
    c.AMU = 1.66053906660e-27; % Atomic mass unit, kg/Dalton  
%     c.AVO = 6.02214076e23;     % Avogadros Number, N/mole
%     c.CC = 299792458.0;        % Speed of Light in Vacuum, m/s 
    c.VMM_PER_TD = 1.0e-21;    % Vm^2 per Townsend
    c.GAMMA = 5.930969584768013e5; % sqrt(2.0 .* c.QE ./ c.ME);  % Factor for converting E to V
    c.EPS0 = 8.8541878128e-12; % Vacuum permittivity F/m, s^4 A^2 kg^-1 m^-3 
end