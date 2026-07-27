function root = projectRoot()
%PROJECTROOT  Absolute path of the project root.
%
%   root = projectRoot()
%
% This file lives in <root>/src, so the project root is the parent of the
% folder that contains it. Used by main_risk_scenarios to locate the data/,
% cache/ and figures/ folders regardless of the current working directory.
%
% If the code is kept in a flat layout (this file sitting directly in the
% project root, with no src/ level), the parent would be wrong, so we only
% step up when a sibling data/ or cache/ folder is actually there.

    here = fileparts(mfilename('fullpath'));
    up   = fileparts(here);
    if ~isempty(up) && (exist(fullfile(up,'data'),'dir') || ...
                        exist(fullfile(up,'cache'),'dir') || ...
                        exist(fullfile(up,'src'),'dir'))
        root = up;                 % <root>/src/projectRoot.m  -> <root>
    else
        root = here;               % flat layout: code sits in the root
    end
end
