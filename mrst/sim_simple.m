function x = sim_simple(cartDims, physDims, tf, verb)
   g = computeGeometry(cartGrid(cartDims, physDims));

   rock = struct('perm', repmat(1, [g.cells.num, 1]), ...
                 'poro', repmat(1, [g.cells.num, 1]));

   T = computeTrans(g, rock);

   fluid = initSimpleFluid('n'  , [   1,   1], ...
                           'mu' , [   1,  30], ...
                           'rho', [1000, 800]);

   gravity reset off

   src = addSource([], [1, g.cells.num], [1, -1], 'sat', [ 1, 0 ; 0, 1]);

   state = initState(g, [], 0, repmat([ 0.5, 0.5 ], [g.cells.num, 1]));
   state = incompTPFA(state, g, T, fluid, 'src', src);

   if nargin < 4, verb = false; end
   state = implicitTransport(state, g, tf, rock, fluid, ...
                             'src', src, 'verbose', verb);

   x = struct('g', g, 'rock', rock, 'T', T, 'fluid', fluid, ...
              'src', src, 'state', state);
end
