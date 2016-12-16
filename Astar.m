function [path, distance] =  Astar(grid, start, goal, varargin)
%Thetastar calculates a path and its distance between two nodes in a grid graph
%   path = Astar(grid, start, goal), calculates a path between tile 'start'
%   and tile 'goal' in 'grid'. All 8 neighbours of every tile is considered
%   connected. 'grid' must be a real matrix where every element represent the
%   cost of moving there. The heuristic used is the diagonal distance.
%
%   A grid has nodes with indices in column major order.
%   For example, tiles in a 3 x 3 size grid has the following indicies.
%
%   -------------------
%   |  1  |  4  |  7  |
%   -------------------
%   |  2  |  5  |  8  |
%   -------------------
%   |  3  |  6  |  9  |
%   -------------------
%
%   In order to convert cartesian coordinates use MATLAB'S function
%   sub2ind.
%
%   [path, distance] = Astar(grid, start, goal) returns the path between
%   tile 'start' and 'goal' in 'grid' and returns the distance of the path.
%
%   [path, distance] = Astar(grid, start, goal, 'heuristic', VALUE) can
%   be used to choose another heuristic. Implemented heuristics are
%
%   'manhattan'     The manhattan distance is used. (This is default for 4 connected neighbours).
%   'diagonal'      The diagonal distance is used (This is default for 8 connected neighbours).
%   'euclidean'     The euclidean distance is used.
%   'djikstra'      No heuristic is used and the Astar becomes the djikstra algorithm.
%
%   [path, distance] = Astar(grid, start, goal, 'cost', C) can be used
%   to modify the cost of moving one tile. C must be a two-dimensional vector
%   where C(1) is the cost of horizontal and vertical movement and C(2)
%   is the cost of diagonal movement.


[heuristic_function, connectedNghbrs, cost] = parse_inputs(grid, start, goal, varargin{:});

gridSize = size(grid);
parent = zeros(gridSize, 'uint32');
closed = false(gridSize);
open = false(gridSize);

gscore = Inf(gridSize);
fscore = zeros(gridSize);

open(start) = true;
gscore(start) = 0;

while any(any(open))
    s = lowestFscore;
    open(s) = false;
    
    if s == goal
        [path, distance] = reconstruct_path(start, goal);
        return
    end
    
    closed(s) = true;
    nghbrs = getNeighbours(s);
    for u = 1 : length(nghbrs)
        sprim = nghbrs(u);
        if ~closed(sprim)
            updateVertex(s, sprim)
        end
    end
    
end

path = [];
distance = 0;

    function s = lowestFscore
        [~, i] = min(fscore(open));
        lowest = find(open, i);
        s = lowest(end);
    end

    function updateVertex(s, sprim)
        g_old = gscore(sprim);
        computeCost(s, sprim);
        if gscore(sprim) < g_old
            fscore(sprim) = gscore(sprim) + ...
                heuristic_function(sprim, goal, grid, cost);
            open(sprim) = true;
        end
    end

    function computeCost(s, sprim)
        d = gscore(s) + c(s, sprim);
        if d < gscore(sprim)
            parent(sprim) = s;
            gscore(sprim) = d;
        end
    end

    function d = c(s, sprim)
        d = dist(s, sprim) + grid(sprim); 
    end

    function d = dist(s, sprim)
        [i, j] = ind2sub(gridSize, s);
        [k, l] = ind2sub(gridSize, sprim);
        if i == k || j == l
            d = cost(1);
        else
            d = cost(2);
        end
    end

    function [total_path, d] = reconstruct_path(start, goal)
        total_path = goal;
        d = 0;
        current = goal;
        while current ~= start
            d = d + dist(current, parent(current));
            current = parent(current);
            total_path = [current total_path];
        end
    end

    function neighbours = getNeighbours(n)
        
        [i, j] = ind2sub(gridSize, n);
        neighbours = zeros(1, 8);
        index = 1;
        for p = i - 1 : i + 1
            for q = j - 1 : j + 1
                if ~(p == i && q == j) && p > 0 && q > 0 && ...
                        p <= gridSize(1) && q <= gridSize(2) ...
                        && grid(n) ~= Inf
                    neighbours(index) = sub2ind(gridSize, p, q);
                    index = index + 1;
                end
            end
        end
        neighbours(neighbours == 0) = [];
    end


    function [heuristic, connectedNeighbours, cost] = ...
            parse_inputs(grid, start, goal, varargin)
        
        p = inputParser;
        
        defaultCost = [10 14];
        defaultConnectedNeighbours = 8;
        defaultHeuristic = 'diagonal';
        expectedHeuristics = {'manhattan','diagonal','djikstra','euclidean'};
        
        addRequired(p, 'grid', @checkGrid);
        addRequired(p, 'start', @checkNode);
        addRequired(p, 'goal', @checkNode);
        addOptional(p, 'cost', defaultCost, @checkCost);
        addOptional(p, 'connectedNeighbours', defaultConnectedNeighbours, @checkNeighbours);
        addOptional(p, 'heuristic', defaultHeuristic, @(x) any(validatestring(x,expectedHeuristics)));
        
        parse(p, grid, start, goal, varargin{:});
        
        cost = p.Results.cost;
        heuristic = pairHeuristic(p.Results.heuristic);
        connectedNeighbours = p.Results.connectedNeighbours;
        
        function bool = checkCost(cost)
            if ~isreal(cost)
                error('Cost must be real.');
            elseif length(cost) > 2 || isempty(cost)
                error('length(cost) must be 1 or 2.');
            else
                bool = true;
            end
        end
        
        function bool = checkGrid(grid)
            if ~ismatrix(grid)
                error('Grid must be represented as a matrix.');
            elseif sum(size(grid)) <= 2
                error('Grid must be a matrix with at least two elements.');
            elseif ~isreal(grid)
                error('Grid must be a real matrix.') ;
            else
                bool = true;
            end
        end
        
        function bool = checkNode(node)
            if node < 1 || node > numel(grid)
                error('Node index must be a number between 1 and numel(grid)');
            else
                bool = true;
            end
        end
        
        function bool = checkNeighbours(value)
            if value == 4
                defaultHeuristic = 'manhattan';
            elseif value == 8
                defaultHeuristic = 'diagonal';
            else
                error('Connected neighbours must be 4 or 8.');
            end
            bool = true;
        end
        
        
        
    end

    function heuristic_func = pairHeuristic(heuristic)
        if strcmp(heuristic, 'diagonal')
            heuristic_func = @diagonal_heuristic;
        elseif strcmp(heuristic, 'manhattan')
            heuristic_func = @manhattan_heuristic;
        elseif strcmp(heuristic, 'euclidean')
            heuristic_func = @euclidean_heuristic;
        elseif strcmp(heuristic, 'djikstra')
            heuristic_func = @(current, goal, grid, cost) 0;
        else
            error('Invalid heuristic option.')
        end
    end

end


function h = diagonal_heuristic(current, goal, grid, cost)
[i, j] = ind2sub(size(grid), current);
[i_g, j_g] = ind2sub(size(grid), goal);
dx = abs(i - i_g);
dy = abs(j - j_g);
h = cost(1) * (dx + dy) + (cost(2) - 2 * cost(1)) * min(dx, dy);
end

function h = euclidean_heuristic(current, goal, grid, cost)
[i, j] = ind2sub(size(grid), current);
[i_g, j_g] = ind2sub(size(grid), goal);
h = cost(1) * sqrt((i_g - i)^2 + (j_g - j)^2);
end

function h = manhattan_heuristic(current, goal, grid, cost)
[i, j] = ind2sub(size(grid), current);
[i_g, j_g] = ind2sub(size(grid), goal);
dx = abs(i - i_g);
dy = abs(j - j_g);
h = cost(1) * (dx + dy);
end

