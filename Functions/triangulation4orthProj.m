function triFace = triangulation4orthProj(Pcell, edgeConON, figON)
% triangulation for orthogonal projection
% edgeConON: Edge Constraints/Condition ON

P_3D = cell2mat(Pcell);
P_2D = orthProj(P_3D);

if edgeConON
    edgeConNum = cumsum(cellfun(@(x) size(x, 2), Pcell));
    edgeCon = [1:size(P_2D, 2); 2:size(P_2D, 2)+1]';
    edgeCon(edgeConNum, 2) = [0, edgeConNum(1:end-1)] + 1;
    
    DT = delaunayTriangulation(P_2D', edgeCon);
    triFace = DT.ConnectivityList(isInterior(DT), :);
else
    DT = delaunayTriangulation(P_2D');
    triFace = DT.ConnectivityList;
end

if figON
    figure;
    triplot(triFace, DT.Points(:, 1), DT.Points(: ,2), 'Color', 'k');
    axis image;
    axis off;
    set(gca, 'LineWidth', 1, 'FontName', 'Calibri', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'Position', [0.1, 0.1, 0.8, 0.8]);
    set(gcf, 'Position', [0, 0, 10*diff(xlim), 10*diff(ylim)]);
end


end

