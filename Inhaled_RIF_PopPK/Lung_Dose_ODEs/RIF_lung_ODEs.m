function dC = RIF_lung_ODEs(~, C, kdiss, kr, kF, effRB, effRA, br_frac, CL, fR, phys, pt)
% RIF_LUNG_ODES - Sets up the lung dose ODEs for rifampin.
%
% INPUTS:
% - C (double array): Contains the state vector of current
%   concentrations/amounts of drug in each compartment
% - kdiss (double): Dissolution constant for rifampin in the lungs (1/h)
% - kr (double): Gut reabsorption rate (1/h)
% - kF (double): Gut transit rate (1/h)
% - effRB (double): Efflux ratio for the bronchi
% - effRA (double): Efflux ratio for the alveoli
% - br_frac (double): Fraction of the drug apportioned to the bronchi for
%   each lung dose
% - CL (double): Systemic clearance (L/h)
% - fR (double): Fractional renal clearance
% - phys (struct): Contains all relevant patient physiological volume and 
%   flow parameters
% - pt (struct): Contains all relevant partition coefficients (not randomly
%   sampled)
%
% OUTPUTS:
% dC (double array): Derivative of the concentration state vector

% Indexing for compartments:
% 1: Venous blood
% 2: Arterial blood
% 3: Lung
% 4: Pleura
% 5: Brain
% 6: Adipose
% 7: Heart
% 8: Muscle
% 9: Skin
% 10: Others
% 11: Bone
% 12: Spleen
% 13: Kidney
% 14: Gut
% 15: Liver
% 16: Lymph Node
% 17: Gut Lumen
% 18: Bronchiolar Epithelial Lining Fluid (bELF)
% 19: Alveolar Epithelial Lining Fluid (aELF)
% 20: Drug Absorption/Dissolution

%% Parameters

Q = phys.Q; % Tissue:blood flow rates
L = phys.L; % Tissue:lymph flow rates
V = phys.V; % Volumes of compartments

dC = zeros(20, 1);

%% Compartment  ----------------------------------------  Terms

% Venous Blood                                          
dC(1) = (1/V.V) * (...
        (Q.Br - L.Br) * C(5)/pt.Tissue.Brain + ...      % brain --> v. blood
        (Q.Ad - L.Ad) * C(6)/pt.Tissue.Adipose + ...    % adipose tissue --> v. blood
        (Q.Hr - L.Hr) * C(7)/pt.Tissue.Heart + ...      % heart --> v. blood
        (Q.Mu - L.Mu) * C(8)/pt.Tissue.Muscle + ...     % muscle --> v. blood
        (Q.Sk - L.Sk) * C(9)/pt.Tissue.Skin + ...       % skin --> v. blood
        (Q.Oth - L.Oth) * C(10)/pt.Tissue.Others + ...  % other tissue --> v. blood
        (Q.Bo) * C(11)/pt.Tissue.Bone + ...             % bone --> v. blood
        (Q.Kd - L.Kd) * C(13)/pt.Tissue.Kidney + ...    % kidney --> v. blood
        (Q.Li - L.Li) * C(15)/pt.Tissue.Liver + ...     % liver --> v. blood
        (L.LN * C(16)/pt.LN) - ...                      % lymph --> v. blood
        (Q.total * C(1)));                              % v. blood --> lungs

% Arterial Blood
dC(2) = (1/V.A) * (...
        (Q.total - L.Lu) * C(3)/pt.Lu - ... % lung --> a. blood
        (Q.total - L.Lu) * C(2));           % a. blood --> all other tissue

% Lung
dC(3) = (1/V.Lu) * (...
        Q.bELF * effRB * C(18) + ...        % bELF --> lung
        Q.aELF * effRA * C(19) - ...        % aELF --> lung
        Q.bELF * C(3) - ...                 % lung --> bELF
        Q.aELF * C(3) + ...                 % lung --> aELF
        Q.total * C(1) - ...                % v. blood --> lung
        (Q.total - L.Lu) * C(3)/pt.Lu - ... % lung --> a. blood
        (L.Lu - Q.Pl) * C(3)/pt.Lu - ...    % lung --> lymph
        Q.Pl * C(3)/pt.Lu);                 % lung --> pleura

% Pleura
dC(4) = (1/V.Pl) * (...
        Q.Pl * C(3)/pt.Lu - ... % lung --> pleura
        Q.Pl * C(4));           % pleura --> lung

% Tissues with lymph flow (Brain, Adipose, Heart, Muscle, Skin, Others)
% Brain
dC(5) = (1/V.Brain) * (...
        Q.Br * C(2) - ...                               % a. blood --> brain
        (Q.Br - L.Br) * C(5)/pt.Tissue.Brain - ...      % brain --> v. blood
        L.Br * C(5)/pt.Tissue.Brain);                   % brain --> lymph

% Adipose tissue
dC(6) = (1/V.Adipose) * (...
        Q.Ad * C(2) - ...                               % a. blood --> adipose
        (Q.Ad - L.Ad) * C(6)/pt.Tissue.Adipose - ...    % adipose --> v. blood
        L.Ad * C(6)/pt.Tissue.Adipose);                 % adipose --> lymph

% Heart
dC(7) = (1/V.Heart) * (...
        Q.Hr * C(2) - ...                               % a. blood --> heart
        (Q.Hr - L.Hr) * C(7)/pt.Tissue.Heart - ...      % heart --> v. blood
        L.Hr * C(7)/pt.Tissue.Heart);                   % heart --> lymph

% Muscle
dC(8) = (1/V.Muscle) * (...
        Q.Mu * C(2) - ...                               % a. blood --> muscle
        (Q.Mu - L.Mu) * C(8)/pt.Tissue.Muscle - ...     % muscle --> v. blood
        L.Mu * C(8)/pt.Tissue.Muscle);                  % muscle --> lymph

% Skin
dC(9) = (1/V.Skin) * (...
        Q.Sk * C(2) - ...                           % a. blood --> skin
        (Q.Sk - L.Sk) * C(9)/pt.Tissue.Skin - ...   % skin --> v. blood
        L.Sk * C(9)/pt.Tissue.Skin);                % skin --> lymph

% Other tissue
dC(10) = (1/V.Others) * (...
        Q.Oth * C(2) - ...                              % a. blood --> others
        (Q.Oth - L.Oth) * C(10)/pt.Tissue.Others - ...  % others --> v. blood
        L.Oth * C(10)/pt.Tissue.Others);                % others --> lymph

% Tissues with no lymph flow (Bone and Spleen)
% Bone 
dC(11) = (1/V.Bone) * (...
        Q.Bo * C(2) - ...               % a. blood --> bone
        Q.Bo * C(11)/pt.Tissue.Bone);   % bone --> v. blood

% Spleen
dC(12) = (1/V.Spleen) * (...
        Q.Sp * C(2) - ...                   % a. blood --> spleen
        Q.Sp * C(12)/pt.Tissue.Spleen);     % spleen --> liver

% Kidney (renal clearance)
dC(13) = (1/V.Kidney) * (...
        Q.Kd * C(2) - ...                               % a. blood --> kidney
        (Q.Kd - L.Kd) * C(13)/pt.Tissue.Kidney - ...    % kidney --> v. blood
        L.Kd *  C(13)/pt.Tissue.Kidney - ...            % kidney --> lymph
        (fR * CL) * C(2));                              % kidney clearance (urine)

% Gut (absorption + reabsorption)
dC(14) = (1/V.Gut) * (...
        Q.Gu * C(2) - ...                           % a. blood --> gut
        (Q.Gu - L.Gu) * C(14)/pt.Tissue.Gut - ...   % gut --> v. blood
        L.Gu *  C(14)/pt.Tissue.Gut + ...           % gut --> lymph
        kr * C(17));                                % gut lumen --> gut

% Liver (hepatic clearance)
dC(15) = (1/V.Liver) * (...
        Q.LA * C(2) + ...                               % a. blood --> liver
        Q.Sp * C(12)/pt.Tissue.Spleen + ...             % spleen --> liver
        (Q.Gu - L.Gu) * C(14)/pt.Tissue.Gut - ...       % gut --> liver
        (Q.Li - L.Li) * C(15)/pt.Tissue.Liver - ...     % liver --> v. blood
        L.Li * C(15)/pt.Tissue.Liver - ...              % liver --> lymph
        (1 - fR) * CL * (Q.LA * C(2) + ...              % liver --> gut lumen
        Q.Sp * C(12)/pt.Tissue.Spleen + ...             % liver --> gut lumen
        (Q.Gu - L.Gu) * C(14)/pt.Tissue.Gut) / Q.Li);   % liver --> gut lumen

% Lymph Node
dC(16) = (1/V.LN) * (...
        L.Br * C(5)/pt.Tissue.Brain + ...       % brain --> lymph
        L.Ad * C(6)/pt.Tissue.Adipose + ...     % adipose --> lymph
        L.Hr * C(7)/pt.Tissue.Heart + ...       % heart --> lymph
        L.Mu * C(8)/pt.Tissue.Muscle + ...      % muscle --> lymph
        L.Sk * C(9)/pt.Tissue.Skin + ...        % skin --> lymph
        L.Oth * C(10)/pt.Tissue.Others + ...    % others --> lymph
        L.Kd * C(13)/pt.Tissue.Kidney + ...     % kidney --> lymph
        L.Gu * C(14)/pt.Tissue.Gut + ...        % gut --> lymph
        L.Li * C(15)/pt.Tissue.Liver + ...      % liver --> lymph
        (L.Lu - Q.Pl) * C(3)/pt.Lu + ...        % lungs --> lymph
        Q.Pl * C(4) - ...                       % pleura --> lymph
        L.LN * C(16)/pt.LN);                    % lymph --> v. blood

% Gut Lumen
dC(17) = (1 - fR) * CL * (...
        (Q.LA * C(2) + ...                                  % a. blood --> gut lumen
        Q.Sp * C(12)/pt.Tissue.Spleen + ...                 % spleen --> gut lumen
        (Q.Gu - L.Gu) * C(14)/pt.Tissue.Gut) / Q.Li) - ...  % gut lumen --> liver
        kr * C(17) - ...                                    % gut lumen --> gut (reabsorption)
        kF * C(17);                                         % gut lumen clearance (feces)

% Bronchi Epithelial Lining Fluid (bELF)
dC(18) = kdiss * C(20) * br_frac + ...  % dissolution --> bELF
        Q.bELF * C(3) - ...             % lung --> bELF
        Q.bELF * effRB * C(18);         % bELF --> lung

% Alveolar Epithelial Lining Fluid (aELF)
dC(19) = kdiss * C(20) * (1 - br_frac) + ...    % dissolution --> aELF
        Q.aELF * C(3) - ...                     % lung --> ELF
        Q.aELF * effRA * C(19);                 % aELF --> lung

% Drug Dissolution
dC(20) = -kdiss * C(20);   % dissolution --> ELF

end