
% To show the figure zoom and tool buttons on v2018+
try
    set(groot,'defaultFigureCreateFcn',@(fig, ~)addToolbarExplorationButtons(fig))
catch 
    [];
end