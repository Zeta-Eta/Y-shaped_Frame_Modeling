close all;
clear;
clc;

addpath(genpath('Functions'));

loadPath = 'hole&plane_JSON_data4test\';

triTopProjFigON = 0; % Triangulation of Top Faces' Orthographic Projection
triBotProjFigON = 0; % Triangulation of Bottom Faces' Orthographic Projection
oriPatchFigON = 1; % Original Patch Plotting (coplanar edges were deleted)
triPatchFigON = 1; % Triangulated Patch Plotting

saveON = 0;
rotyON = 1;
gmON = 1;
saveFileName = 'Y-shaped_Frame_model4test_V1.stl';

lenUnit = 'mm';

%% Size Information Setting
% unit: mm | °(degree) | radian

HHnum = 3; % number of head holders
Hnum = 4; % number of holes in each head holder

% Head Holders
sizeInf.HH.len = 15;
sizeInf.HH.wid = 12;
sizeInf.HH.lenH = 11;
sizeInf.HH.widH = 8;

% Holes
sizeInf.H.r = [2.5, 2.5]; % radii of holes
sizeInf.H.n = 2^5; % edge number of each hole
sizeInf.H.locXYofTopH = [0, 12.5, 0; 0, -7.5, 0]'; % Locations(only need to set X & Y) of Top Holes

% Y-shaped Frame
sizeInf.Y.angle = [16; 0; 8]; % angles between head holders and slant legs
sizeInf.Y.theta = deg2rad(sizeInf.Y.angle); % degree to radian
sizeInf.Y.len = [2; 2; 2]; % additional length
sizeInf.Y.wid = [20; 22; 22]; % width of slant legs
sizeInf.Y.thk = 5 + zeros(HHnum, 1); % thickness
sizeInf.Y.hig = 60 + zeros(HHnum, 1); % height from XY-plane to the top face

sizeInf.Y = struct2table(sizeInf.Y);

%% Data Loading & Preprocessing

% Coordinate Transformation Matrix
T = [
    1,  0,  0;
    0,  0,  1;
    0, -1,  0];

rotY2Z = 0; % degree || CCW: from Y to Z
R = rotx(rotY2Z);
% R = [ 1,  0,  0;
%       0,  cos(deg2rad(rotYZ)), -sin(deg2rad(rotYZ));
%       0,  sin(deg2rad(rotYZ)),  cos(deg2rad(rotYZ))];
T = R*T;

files.hole = dir([loadPath, '*-*.json']);
files.plane = dir([loadPath, 'P*.json']);

initHCP = NaN(3, Hnum, HHnum); % initial Holes' Central Points
for n = 1:length(files.hole)
    tempName = files.hole(n).name;
    tempJSON = jsondecode(fileread([loadPath, tempName]));
    initHCP(:, str2double(tempName(3)), str2double(tempName(1))) = T*tempJSON.markups.controlPoints.position;
end

initPAP = NaN(3, HHnum); % initial Planes' Anchor Points
normVec = NaN(3, HHnum); % Planes' Normal Vectors
for n = 1:length(files.plane)
    tempJSON = jsondecode(fileread([loadPath, files.plane(n).name]));
    initPAP(:, n) = T*tempJSON.markups.center;
    normVec(:, n) = T*tempJSON.markups.normal;
end

% initPAP to PAP = [0, 0, ?] (intercept)
PAP = orthPos(zeros(2, HHnum), initPAP, normVec);

% initHCP to HCP (on the specific plane)
HCP = setHCP(initHCP, PAP, normVec);

%% HCP Fitting
% 'fminsearch': Nelder-Mead Method (derivative-free) 
% * usually has higher precision for high-dimensional optimization problems
% 'fminunc': quasi-Newton Method (derivative-based) 
% * usually has higher speed for high-dimensional optimization problems

fittingMethod = 'fminsearch'; 

RSS = NaN(1, HHnum); % Residual Sum of Squares
fitdHCP = NaN(size(HCP));
baseVec = NaN(3, 2, HHnum);
for n = 1:HHnum
    
    func = @(params) fitHeadHolder(params, HCP(:, :, n), PAP(:, n), normVec(:, n), sizeInf);
    
    initParams = [HCP(1:2, 1, n); ...
        angle(diff(HCP(1, 1:2, n)) + 1i.*diff(HCP(2, 1:2, n)))];
    
    tic;
    if strcmp(fittingMethod, 'fminsearch')
        % Use 'fminsearch'
        options = optimset('fminsearch');
        options.TolFun = 1e-7;
        options.TolX = 1e-7;
        options.MaxFunEvals = 1e5;
        options.Display = 'notify';
        fitdParams = fminsearch(func, initParams, options);
    elseif strcmp(fittingMethod, 'fminunc')
        % Use 'fminunc'
        options = optimoptions('fminunc');
        options.TolFun = 1e-7;
        options.TolX = 1e-7;
        options.MaxFunctionEvaluations = 1e5;
        options.Display = 'notify-detailed';
        fitdParams = fminunc(func, initParams, options);
    end
    toc;
    
    [RSS(n), fitdHCP(:, :, n), baseVec(:, :, n)] = func(fitdParams);
    
end

HHCP = squeeze(mean(fitdHCP, 2)); % Head Holders' Central Points (could be compared to initPAP)

W = [
    -1,  1,  1, -1;
    -1, -1,  1,  1];

HHV = NaN(size(fitdHCP)); % Head Holders' Vertices (top face)
for n = 1:HHnum
    HHV(:, :, n) = HHCP(:, n) + (0.5.*[sizeInf.HH.len, sizeInf.HH.wid].*baseVec(:, :, n))*W; 
end

%% Vertices and Faces Getting
baseP = [mean(HHV(:, 1:2, 1), 2, 'omitnan'), HHV(:, 4, 2), HHV(:, 3, 3)];
cntrP = HHCP;

% Bottom Vertices and Faces
botVnum = 5;
botFace.Bevel = botVnum*((1:HHnum)'-1) + (1:4);
botFace.Level = botVnum*((1:HHnum)'-1) + (4:6);
botFace.Level(end, end) = 1;
botFace.LevelVec = reshape(botFace.Level', 1, []);

botVrtx = NaN(3, botVnum*HHnum);
botVrtxNum = size(botVrtx, 2);
normVec4Bevel = NaN(3, HHnum);
drctVec4Bevel = NaN(3, HHnum);
drctVec4Level = NaN(3, HHnum);
for n = 1:HHnum
    [botVrtx(:, botFace.Bevel(n, :)), normVec4Bevel(:, n), ...
        drctVec4Bevel(:, n), drctVec4Level(:, n)] = ...
        getBevel(baseP(:, n), cntrP(:, n), normVec(:, n), sizeInf.Y(n, :), 0);
end

for n = 1:HHnum
    inPid = mod([n, n+1], HHnum);
    inPid(inPid == 0) = HHnum;
    botVrtx(:, botFace.Level(n, 2)) = getLevel(botVrtx(:, ...
        [botFace.Level(n, 1), botFace.Level(n, end)]), ...
        drctVec4Level(:, inPid));
end

% Top Vertices and Faces
topVnum = botVnum + 4;
topFace.Bevel = topVnum*((1:HHnum)'-1) + (2:7);
topFace.Level = topVnum*((1:HHnum)'-1) + (7:11);
topFace.Level(end, end-1:end) = 1:2;
topFace.LevelVec = reshape(topFace.Level', 1, []);

topVrtx = NaN(3, topVnum*HHnum);
topVrtxNum = size(topVrtx, 2);
topVrtx(:, reshape(topFace.Level(:, 2:4)', 1, [])) = ...
    botVrtx(:, botFace.LevelVec) + sizeInf.Y.thk(1).*[0, 0, 1]';
for n = 1:HHnum
    topVrtx(:, topFace.Bevel(n, 2:end-1)) = ...
        botVrtx(:, botFace.Bevel(n, :)) + sizeInf.Y.thk(1).*normVec4Bevel(:, n);
    topVrtx(:, topFace.Bevel(n, 1)) = ...
        getCorner(topVrtx(:, topFace.Bevel(n, 1) - 1), ...
        [0, 0, 1]', normVec4Bevel(:, n), sizeInf.Y.thk(1));
    topVrtx(:, topFace.Bevel(n, end)) = ...
        getCorner(topVrtx(:, topFace.Bevel(n, end) + 1), ...
        [0, 0, 1]', normVec4Bevel(:, n), sizeInf.Y.thk(1));
end

% Top Vertices and Faces
topVnum = botVnum + 4;
topFace.Bevel = topVnum*((1:HHnum)'-1) + (2:7);
topFace.Level = topVnum*((1:HHnum)'-1) + (7:11);
topFace.Level(end, end-1:end) = 1:2;
topFace.LevelVec = reshape(topFace.Level', 1, []);

topVrtx = NaN(3, topVnum*HHnum);
topVrtxNum = size(topVrtx, 2);
topVrtx(:, reshape(topFace.Level(:, 2:4)', 1, [])) = ...
    botVrtx(:, botFace.LevelVec) + sizeInf.Y.thk(1).*[0, 0, 1]';
for n = 1:HHnum
    topVrtx(:, topFace.Bevel(n, 2:end-1)) = ...
        botVrtx(:, botFace.Bevel(n, :)) + sizeInf.Y.thk(1).*normVec4Bevel(:, n);
    topVrtx(:, topFace.Bevel(n, 1)) = ...
        getCorner(topVrtx(:, topFace.Bevel(n, 1) - 1), ...
        [0, 0, 1]', normVec4Bevel(:, n), sizeInf.Y.thk(1));
    topVrtx(:, topFace.Bevel(n, end)) = ...
        getCorner(topVrtx(:, topFace.Bevel(n, end) + 1), ...
        [0, 0, 1]', normVec4Bevel(:, n), sizeInf.Y.thk(1));
end

% Hole's Vertices and Faces
holeVrtx.botLevel = getHole(sizeInf.H.locXYofTopH, [0, 0, 1]', ...
    [0, 0, sizeInf.Y.hig(1)]', [0, 0, 1]', ...
    [1, 0, 0; 0, -1, 0]', sizeInf.H.r(1), sizeInf.H.n);
holeVrtx.topLevel = getHole(sizeInf.H.locXYofTopH, [0, 0, 1]', ...
    [0, 0, sizeInf.Y.hig(1) + sizeInf.Y.thk(1)]', [0, 0, 1]', ...
    [1, 0, 0; 0, -1, 0]', sizeInf.H.r(1), sizeInf.H.n);

holeVrtx.botBevel = NaN(3, Hnum.*sizeInf.H.n, HHnum);
holeVrtx.topBevel = NaN(3, Hnum.*sizeInf.H.n, HHnum);
for n = 1:HHnum
    holeVrtx.botBevel(:, :, n) = getHole(fitdHCP(:, :, n), normVec(:, n), ...
        botVrtx(:, botFace.Bevel(n, 1)), normVec4Bevel(:, n), ...
        baseVec(:, :, n), sizeInf.H.r(2), sizeInf.H.n);
    holeVrtx.topBevel(:, :, n) = getHole(fitdHCP(:, :, n), normVec(:, n), ...
        topVrtx(:, topFace.Bevel(n, 1)), normVec4Bevel(:, n), ...
        baseVec(:, :, n), sizeInf.H.r(2), sizeInf.H.n);
end

holVrtx = [holeVrtx.botLevel, reshape(holeVrtx.botBevel, 3, []), ...
    holeVrtx.topLevel, reshape(holeVrtx.topBevel, 3, [])];
holVrtxNum = size(holVrtx, 2);
holeNum = holVrtxNum./sizeInf.H.n;

holFace.ID = [zeros(1, size(sizeInf.H.locXYofTopH, 2)), ...
    reshape(repmat(1:HHnum, Hnum, 1), 1, [])]';
% Level = 0 | Bevel_n = n

holFace.all = reshape(1:holVrtxNum, sizeInf.H.n, [])';
holFace.all = [holFace.all, holFace.all(:, 1)];
holFace.bot = holFace.all(1:end/2, :);
holFace.top = holFace.all(end/2 + 1:end, :);

% Combination of Bottom, Top and Holes' Vertices
Vrtx = [topVrtx, botVrtx, holVrtx];

botFace.Bevel = botFace.Bevel + topVrtxNum;
botFace.Level = botFace.Level + topVrtxNum;
botFace.LevelVec = botFace.LevelVec + topVrtxNum;

holFace.all = holFace.all + topVrtxNum + botVrtxNum;
holFace.bot = holFace.bot + topVrtxNum + botVrtxNum;
holFace.top = holFace.top + topVrtxNum + botVrtxNum;

% Side Faces
sidFace.Corner = getSideFace(botFace.Bevel(:, [1, end]), topFace.Bevel(:, [1, end]), 'Corner');
sidFace.Rectangle.Bevel = getSideFace(botFace.Bevel, topFace.Bevel(:, 2:end-1), 'Rectangle');
sidFace.Rectangle.Level = getSideFace(botFace.Level, topFace.Level(:, 2:end-1), 'Rectangle');
sidFace.Rectangle.Hole = getSideFace(flip(holFace.bot, 2), flip(holFace.top, 2), 'Rectangle');
sidFace.all = [sidFace.Corner;
    sidFace.Rectangle.Bevel;
    sidFace.Rectangle.Level;
    sidFace.Rectangle.Hole];

FaceCell = {
    topFace.LevelVec;
    flip(botFace.LevelVec, 2);
    topFace.Bevel;
    flip(botFace.Bevel, 2);
    sidFace.all};

maxFnum = max(cellfun(@(x) size(x, 2), FaceCell));

tempFunc = @(X) [X, NaN(size(X, 1), maxFnum - size(X, 2))];
Face = cell2mat(cellfun(tempFunc, FaceCell, 'UniformOutput', false));

%% Triangulation
TriFace.topLevel = getTriFace4TB(0, 1, Vrtx, topFace.LevelVec, holFace, 1, triTopProjFigON);
TriFace.botLevel = getTriFace4TB(0, 0, Vrtx, botFace.LevelVec, holFace, 1, triBotProjFigON);
TriFace.topBevel = cell(HHnum, 1);
TriFace.botBevel = cell(HHnum, 1);
for n = 1:HHnum
    TriFace.topBevel{n} = getTriFace4TB(n, 1, Vrtx, topFace.Bevel(n, :), holFace, 1, triTopProjFigON);
    TriFace.botBevel{n} = getTriFace4TB(n, 0, Vrtx, botFace.Bevel(n, :), holFace, 1, triBotProjFigON);
end

normVrtx = [
    0,  1,  2,  1;
    0, -2,  0,  2;
    0,  0,  0,  0];
TriFace.sid = getTriFace4S(normVrtx, sidFace.all, 1, 0, 1);

triFace = [
    TriFace.topLevel;
    TriFace.botLevel;
    cell2mat(TriFace.topBevel);
    cell2mat(TriFace.botBevel);
    cell2mat(TriFace.sid)];


%% Original Patch Plotting (coplanar edges were deleted)
if oriPatchFigON
    figure;
    tiledlayout(1, 1, 'Padding', 'tight');
    nexttile;
    
    hold on;
    
    patch('Faces', Face, 'Vertices', Vrtx', 'FaceAlpha', 0.69, ...
        'FaceColor', 0.69*[1 1 1], 'EdgeColor', 0*[1 1 1]);
    view([70, 55]); % [70, 55] | [165, -75]

    for n = 1:HHnum
        patch('XData', HHV(1, :, n), 'YData', HHV(2, :, n), 'ZData', HHV(3, :, n), ...
            'FaceColor', 0.9*[1 1 1], 'EdgeColor', 0*[1 1 1]);
        tempVec = HHCP(:, n) + [-25, 25 + sizeInf.Y.thk(1)].*normVec(:, n);
        plot3(tempVec(1, :), tempVec(2, :), tempVec(3, :), 'r');
    end
    
    hold off;
    
    axis equal;
    xlim([-50 50] + 5*[-1 1]);
    ylim([-50 50] + 5*[-1 1]);
    zlim([0 100] + 5*[-1 1]);
    grid on;

    set(gca, 'XTick', -50:25:50, 'YTick', -50:25:50, 'ZTick', 0:25:100, ...
        'LineWidth', 1, 'FontName', 'Calibri', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('X', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    ylabel('Y', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    zlabel('Z', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    set(gcf, 'Position', [0, 0, 720, 720]);
    
end

%% Triangulated Patch Plotting
if triPatchFigON
    figure;
    tiledlayout(1, 1, 'Padding', 'tight');
    nexttile;
    
    hold on;
    
    patch('Faces', triFace, 'Vertices', Vrtx', 'FaceAlpha', 0.69, ...
        'FaceColor', 0.69*[1 1 1], 'EdgeColor', 0*[1 1 1]);
    view([70, 55]); % [70, 55] | [165, -75]
    
    for n = 1:HHnum
        patch('XData', HHV(1, :, n), 'YData', HHV(2, :, n), 'ZData', HHV(3, :, n), ...
            'FaceColor', 0.9*[1 1 1], 'EdgeColor', 0*[1 1 1]);
        tempVec = HHCP(:, n) + [-25, 25 + sizeInf.Y.thk(1)].*normVec(:, n);
        plot3(tempVec(1, :), tempVec(2, :), tempVec(3, :), 'r');
    end
    
    hold off;
    
    axis equal;
    xlim([-50 50] + 5*[-1 1]);
    ylim([-50 50] + 5*[-1 1]);
    zlim([0 100] + 5*[-1 1]);
    grid on;

    set(gca, 'XTick', -50:25:50, 'YTick', -50:25:50, 'ZTick', 0:25:100, ...
        'LineWidth', 1, 'FontName', 'Calibri', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('X', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    ylabel('Y', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    zlabel('Z', 'FontName', 'Calibri', 'FontSize', 28, 'FontWeight', 'bold', 'FontAngle', 'italic');
    set(gcf, 'Position', [0, 0, 720, 720]);
    
end

%% STL File Saving
if saveON
    outputVrtx = Vrtx;
    if rotyON
        outputVrtx = roty(180)*outputVrtx + [0, 0, sizeInf.Y.hig(1) + sizeInf.Y.thk(1)]';
    end
    stlwrite(triangulation(triFace, outputVrtx'), saveFileName);
end


%% Geometric Measurements
outputPrcsn = '%.4f';
if gmON
    GM = GMof3Dobj(Vrtx, triFace);
    fprintf(['------------------------------\n', ...
        'Area: ', outputPrcsn, ' ', lenUnit, '²\n', ...
        'Volume: ', outputPrcsn, ' ', lenUnit, '³\n', ...
        'Centroid of Area: (', outputPrcsn, ', ', outputPrcsn, ', ', outputPrcsn, ') ', lenUnit, '\n', ...
        'Centroid of Volume: (', outputPrcsn, ', ', outputPrcsn, ', ', outputPrcsn, ') ', lenUnit, '\n', ...
        '------------------------------\n'], GM.A, GM.V, GM.CofA, GM.CofV);
end




