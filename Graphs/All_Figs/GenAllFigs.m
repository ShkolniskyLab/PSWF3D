function [  ] = GenAllFigs(  )

% This function generates all the plots that appear in the paper. 
% Note : The graphs for L=52 are disabled. In order to enable them,
% uncomment the lines 2,3 below.

fprintf('Starting GenGraphs_ErrBound_AllFigs\n');
t=tic;
GenGraphs_ErrBound_AllFigs(  )
t=toc(t);
fprintf('Done GenGraphs_ErrBound_AllFigs (took %f seconds)\n',t);


fprintf('Starting GenGraphs_ErrBound3_AllFigs\n');
t=tic;
GenGraphs_ErrBound3_AllFigs(  )
t=toc(t);
fprintf('Done GenGraphs_ErrBound3_AllFigs (took %f seconds)\n',t);

fprintf('Starting GenGraphs_ErrBound1_AllFigs\n');
t=tic;
GenGraphs_ErrBound1_AllFigs(  )
t=toc(t);
fprintf('Done GenGraphs_ErrBound1_AllFigs (took %f seconds)\n',t);

fprintf('Starting  GenGraphs_ErrBound2_AllFigs\n');
t=tic;
GenGraphs_ErrBound2_AllFigs(  )
t=toc(t);
fprintf('Done GenGraphs_ErrBound2_AllFigs (took %f seconds)\n',t);

fprintf('Starting GenGramian_Graph_AllFigs\n');
t=tic;
GenGramian_Graph_AllFigs( )
t=toc(t);
fprintf('Done GenGramian_Graph_AllFigs (took %f seconds)\n',t);

fprintf('Starting GenXiLPropList_AllFigs\n');
t=tic;
GenXiLPropList_AllFigs(  )
t=toc(t);
fprintf('Done GenXiLPropList_AllFigs (took %f seconds)\n',t);


end

