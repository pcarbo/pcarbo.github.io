%  CRFLOCALIZE  Estimates the labels on segments of an image using the
%               labels of interest regions. 
%
%    CRFLOCALIZE(SEGS,Y,COORD,SCALE) returns a 1 x K matrix where each
%    entry is the probability that a segment k is labeled positive. K is
%    the number of segments in the image. The input arguments to
%    CRFLOCALIZE are:
%
%      SEGS  H x W matrix, where H is the height of the image and W is
%            its width. Each entry (h,w) = k says that the pixel belongs
%            to segment k, where k ranges from 1 to K.
%      Y     1 x N matrix of label probabilities for the interest
%            regions, where N is the number of interest regions in the
%            image. Each entry is the probability that the interest
%            region is positively labeled (i.e. belongs to the object in
%            question). Generally, these labels will be output from the
%            ssmcmc program.
%      COORD 2 x N matrix of (x,y) coordinates for the interest region 
%            centres.
%      SCALE 1 x N matrix of the character scales (in pixels) of the 
%            interest regions. 

function ys = crflocalize (segs, Y, coord, scale)

  % CRF parameters. You may modify these at your leisure.
  adjacentsegthresh      = 0.01;
  adjacentsegandptthresh = 0.1;
  theta                  = 0.1;
  numsamples             = 10e4;

  % Get the number of segments and the number of interest regions.
  nsegs = max(max(segs));
  npts  = length(Y);
  Y     = [Y; 1-Y];

  % Get the height and width of the image.
  [h w] = size(segs);

  % Compute segment contours and overlap.
  % ------------------------------------ 
  % Calculate the size of the contour for each segment. It is simply the set
  % of pixels that are bordering pixels that don't belong to the segment.
  % Also, compute the sparse weighted adjacency matrix between the
  % segments. The weights are calculated according to how much border is
  % shared between the two segments. The variable "overlap" is the number of
  % pixels shared on the border between any two segments.
  fprintf('Computing segment contours and overlap.\n');
  contours = zeros(nsegs,1);
  segarea  = zeros(nsegs,1);
  overlap  = zeros(nsegs,nsegs);
  
  % Repeat for each image segment.
  for i = 1:nsegs
    Si  = (segs == i);
    Sni = (segs ~= i);
	
    % Calculate the total area of the segment.
    segarea(i) = sum(sum(Si));
    
    % This is the set of indices of segments that are adjacent to the
    % current segment.
    [xi xj] = find((Sni & [zeros(1,w); Si(1:h-1,1:w)]) | ...
		   (Sni & [Si(2:h,1:w); zeros(1,w)]) | ...
		   (Sni & [zeros(h,1) Si(1:h,1:w-1)]) | ...
		   (Sni & [Si(1:h,2:w) zeros(h,1)]));
    contours(i) = length(xi);
	
    % Repeat for each bordering pixel.
    for t = 1:contours(i)
      j            = segs(xi(t),xj(t));
      overlap(i,j) = overlap(i,j) + 1;
    end
  end
      
  % Normalize overlap. Repeat for each pair of segments that have
  % non-zero overlap. Notice that the normalized overlap is a value
  % between 0 and 1.
  [is js] = find(overlap);
  for t = 1:length(is);
    i = is(t);
    j = js(t);
    overlap(i,j) = overlap(i,j) / (2*contours(i)) + ...
         	   overlap(i,j) / (2*contours(j)); 
  end

  % Compute adjacencies matrix segments and the local features.
  % ----------------------------------------------------------
  % The weights are calculated according to how much area is shared
  % between the segment and local feature. Repeat for each segment and
  % local feature.
  fprintf('Computing segment-local feature adjacency matrix.\n');
  overlapfs = zeros(nsegs,npts);
  featarea  = zeros(npts,1);

  % Repeat for each local feature.
  for j = 1:npts
    Si = (segs == i);
	
    % Compute the region occupied by the local feature.
    Sj = zeros(h,w);
    Sj = fillcircle(Sj,coord(1,j),coord(2,j),scale(j));
	
    [x y]       = find(Sj);
    featarea(j) = length(x);
    for t = 1:featarea(j)
      i              = segs(x(t),y(t));
      overlapfs(i,j) = overlapfs(i,j) + 1;
    end
  end

  % Normalize overlap. Repeat for each segment and local feature that have
  % non-zero overlap. Notice that the normalized overlap is a value between
  % 0 and 1.
  [is js] = find(overlapfs);
  for t = 1:length(is)
    i = is(t);
    j = js(t);
    overlapfs(i,j) = overlapfs(i,j) / featarea(j); 
  end

  % Construct the MRF.
  % -----------------
  % Create the adjacency graph.
  a = Sparsegraph(overlap > adjacentsegthresh);
    
  % Build the likelihood potential matrix.
  g = cell(nsegs,1);
  for i = 1:nsegs
    js = find(overlapfs(i,:));
    if length(js)
      g{i} = sum(repmat(overlapfs(i,js),2,1) .* Y(:,js),2);
    else
      g{i} = [0.01 1]';
    end
  end

  % Construct the context potentials.
  m = numedges(a);
  f = cell(m,1);
  for u = 1:m
    [i j] = edge(a,u);
    f{u}  = theta + [overlap(i,j) 0;
		     0 overlap(i,j)];
  end
    
  % Create the set of potentials.
  potentials = Tabular(a,f,g);
  
  % Estimate marginals.
  % ------------------
  fprintf('Running blocked Gibbs sampling.\n');
  tic;
  [bel ans]  = bgsf(a,potentials,numsamples,0);
  fprintf('MCMC took %d seconds.\n', round(toc));

  % Return the final result.
  ys = zeros(1,nsegs);
  for j = 1:nsegs
    ys(j) = bel{j}(1);
  end