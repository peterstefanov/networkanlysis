function ij = kroneckerModel(SCALE, edgefactor) 
 %% Set number of vertices.
  N = 2^SCALE;

  %% Set number of edges.
  M = edgefactor * N;

  %% Set initiator probabilities.
  [A, B, C] = deal (0.51, 0.17, 0.17);

  %% Create index arrays.
  ij = ones (2, M);
  %% Loop over each order of bit.
  ab = A + B;
  c_norm = C/(1 - (A + B));
  a_norm = A/(A + B);

  for ib = 1:SCALE,
    %% Compare with probabilities and set bits of indices.
    ii_bit = rand (1, M) > ab;
    jj_bit = rand (1, M) > ( c_norm * ii_bit + a_norm * not (ii_bit) );
    ij = ij + 2^(ib-1) * [ii_bit; jj_bit];
  end

  %% Permute vertex labels
  p = randperm (N);
  ij = p(ij);

  %% Permute the edge list
  p = randperm (M);
  ij = ij(:, p);

  %% Adjust to zero-based labels.
  ij = ij - 1;
  