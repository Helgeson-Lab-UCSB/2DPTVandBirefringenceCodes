%function [xyzs_id]=track(xyzs,maxdisp,inipos,memory,dim,verbose,goodenough)

%The input xyzs needs to be a matrix of positions and time, where time must
%be in the last column and positions in the first columns. the program will
%sort the data according to the first coluomns, where the number of columns
%taken into account is determined by dim. columns in between will be
%ignored.

%also, the matrix needs to be sorted by time.

function [lub] = trackmem(xyzs,maxdisp,dim,goodenough,memory)
tic
%xyzs=single(xyzs);
%whos xyzs
%sum(xyzs)

dd = length(xyzs(1,:));                   %dd is the no. of colomns in xyzs

% INSTEAD OF THIS PART:
% if not keyword_set(dim) then begin
%    dim = 2 < dd
%    message,' Setting dim = '+strcompress(string(dim))+' by default',/inf
% endif
% if not keyword_set(memory) then memory = 0

%dim=2;                  %set the dimention as 2 by default  !!!!! change later
t = (xyzs(:,dd))';      %this is the time as a row vector
%check the input time vector is ok, i.e. sorted and uniform
st_t=circshift(t,[0,1]);
for i = 2:length(t), st(i-1) = t(i)-st_t(i);, end
%st is now a vector of time differences, so all elements should be equal
%and positive, OR ZERO
if sum(st<0)~=0, warning('ERROR - Time vector out of order!'), end
w = find(st>0);  % index to non-zero time steps, specified at beginning of the time jump
z = length(w);  % could be zero --> case if all times are the same
if z==0, warning('ERROR - All positions are at the same time!'),...        
else, if sum(st(w)-st(w(1)))~=0, warning('WARNING - Time vector gapped or not evenly gridded!'), end, end


% IN IDL THIS WAS Z=Z+1L
z = z+1;
% partition the input data by unique times
res = unq(t,[]);
res = res+1;
res = [1,res,length(t)+1];
% get the initial positions
ngood = res(2) - res(1);
% ngood here is number of features, starting from the beginning of the input features file, 
% which have the same time step 
eyes = [1:ngood];  %and ngood has the indeces of those features

%INSTEAD OF THIS PART:
% if keyword_set( inipos ) then begin
%    pos = inipos(0:dim-1,*)
%    istart = 0L
%    n = n_elements(pos(0,*))
% endif else begin
%    pos = xyzs(0:dim-1,eyes)
%    istart = 1L      ;we don't need to track t=0.
%    n = ngood
% endelse

%DID JUST THIS:
pos = xyzs(eyes,1:dim);
istart = 1;      %we don't need to track t=0.
%IN IDL WAS DEFINED AS 1L, long int 1 --> if keyword inipos set, istart set to 0L, 
% go through loop over z even if z=0 (no non-zero time steps were found); otherwise,
% set to 1, goes through loop only if z >= 1.
n = ngood;

%how long are the 'working' copies of the data?
zspan = 50;
if n>200, zspan = 20;, end
if n>500, zspan = 10;, end

resx = zeros(zspan,n) -1;  
bigresx = zeros(z,n) -1;  % the -1 is a marker, will be replaced by row number in 
% original feature file; presumably if particle can no longer be tracked, -1's 
% remain at end. Change to zeros. !!!
mem = zeros( 1,n );
uniqid = [1:n];    
maxid = n;
olist = [1,1];

if goodenough > 0
   dumphash = zeros( 1,n );
   nvalid = ones( 1,n );
end
% we may not need to track the first time step!

% if not keyword_set(inipos) then begin
%    resx(*,0) = eyes
%    if keyword_set(goodenough) then nvalid = nvalid+1
% endif

resx(1,:) = eyes;     %put the first set of feature indeces in the first row of resx
%set up some nice constants
maxdisq = maxdisp^2;
%SKIPPED: verbose = keyword_set( verbose )

%Use fancy code for large n, small d
notnsqrd = (sqrt(n*ngood) > 200) & (dim < 7);

if notnsqrd
    cube = zeros(3^dim,dim);
%   construct the vertices of a 3x3x3... d-dimensional hypercube    
    for d=1:dim
        numb = 0;
        for j=1:3^(d-1):((3^dim));
            cube(j:j+3^(d-1)-1,d) = numb;
            numb = mod(numb+1,3);
        end
    end
%cube    

    %    calculate a blocksize which may be greater than maxdisp, but which
    %    keeps nblocks reasonably small.
    volume = 1;  % volume in dimensional space.  e.g. dim=2, "volume" is the area
    for d=1:dim
        minn = min(xyzs(w,d));
        maxx = max(xyzs(w,d));
        volume = volume*(maxx-minn);
    end
    
    %maxdisp,volume,ngood,dim
    blocksize = max( maxdisp,(volume/(20*ngood))^(1.0/dim));  %Tailor the factor in bottom for the particular system
end

%   Start the main loop over the frames.
for i=istart+1:z %always starts at 2 (while inipos not implemented) 
    %i=i,StartTime(i)=toc
    ispan = mod (i-1,zspan)+1;

    %   Get the new particle positions.
    m = res(i+1) - res(i);  % number of new particles
    eyes = res(i) + [0:m-1];  % points to the lines in the feature file for the new particles

    if m>0  %  ???? why would we need to do this?  Cannot think of an instance
        xyi=xyzs(eyes,1:dim);     %positions of new particles
        found = zeros(1,m);

        %   THE TRIVIAL BOND CODE BEGINS
        if notnsqrd
            %i=i
            %Use the raster metric code to do trivial bonds

            %   construct "s", a one dimensional parameterization of the space
            %   ( which consists of the d-dimensional raster scan of the volume.)
            abi = floor(xyi./blocksize);
            abpos = floor(pos./blocksize);
            %si = zeros(1,m);
            si=ones(1,m);
            spos = zeros(1,n);
            %size(spos);
            dimm = zeros(1,dim);
            coff = 1;

            for j=1:dim
                minn = min([[abi(:,j)];[abpos(:,j)]]);
                %test=([[abi(:,j)];[abpos(:,j)]]);
                maxx = max([[abi(:,j)];[abpos(:,j)]]);
                abi(:,j) = abi(:,j) - minn;
                abpos(:,j) = abpos(:,j) - minn;
                dimm(j) = maxx-minn+1;
                si = si + abi(:,j)'*coff;
                %test=si(1)
                spos = spos + abpos(:,j)'*coff;
                %size(spos)
                coff = coff*dimm(j);
            end
            nblocks = coff;     % the # of blocks in the volume

            %   trim down (intersect) the hypercube if its too big to fit in the
            %   particle volume. (i.e. if dimm(j) lt 3)
            cub = cube;
            deg = find( dimm < 3 );
            if deg
                for j=1:length(deg);, cub = cub( find( cub( :,deg(j) ) < dimm(deg(j))),: );, end         %STUCK HERE, GO BACK
            end

            %   calculate the "s" coordinates of hypercube (with a corner @ the origin)
            scube = zeros(length(cub(:,1)),1);
            coff = 1;
            for j=1:dim
                %scube,test=cub(:,j),coff
                scube = scube + cub(:,j)*coff;
                coff = coff*dimm(j);
            end

            %   shift the hypercube "s" coordinates to be centered around the origin
            coff = 1;
            for j=1:dim
                if dimm(j) > 3, scube = scube - coff;, end
                coff = coff*dimm(j);
            end
            scube = mod(scube + nblocks , nblocks);

            %   get the sorting for the particles by their "s" positions.
            [nrn,isort] = sort(si);

            %  make a hash table which will allow us to know which new particles
            %   are at a given si.
            strt = zeros(1,nblocks) - 1;
            fnsh = zeros(1,nblocks);
            
            for j=1:m
                if strt(si(isort(j))) == -1;
                    strt(si(isort(j))) = j;
                    fnsh(si(isort(j))) = j;
                else
                    fnsh(si(isort(j))) = j;
                end
            end

            %test1=sum(fnsh>0) , test2=sum(fnsh)    
            %   loop over the old particles, and find those new particles in the 'cube'.
            coltot = zeros( 1,m );
            rowtot = zeros( 1,n );
            which1 = zeros( 1,n );
            for j=1:n
                map = -1; 
                s = mod ( (scube + spos(j)) , nblocks) +1;
% if j==1,s,end                  
                w = find( strt(s) ~= -1); ngood=length(w);
%if j<10, ngood, end                
                if ngood ~= 0
                    s = s(w);
         
                    for k=1:ngood, map = [map,isort(strt(s(k)):fnsh(s(k)))];, end
                    map = map(2:length(map));

                    % find those trivial bonds
                    distq = zeros( length(map) ,1 );
                    for d=1:dim
                        %map,d,j
                        %distq,test1=xyi(map,d),test2=pos(j,d)
                        distq = distq + ( xyi(map,d) - pos(j,d) ).^2;
                    end
                    ltmax = distq < maxdisq;
                    rowtot(j) = sum( ltmax );
                                    
                    if rowtot(j) >= 1;
                        w = find( ltmax );
                        coltot( map(w) ) = coltot( map(w) ) +1;
                        which1(j) = map( w(1) );
                    end
                end
            end

            
            ntrk = floor( n - sum(rowtot == 0));
            w = find( rowtot == 1 );, ngood=length(w);
            %test=w(1:10)
            if ngood ~= 0
                %test=sum(coltot)
                ww = find( coltot( which1(w) ) == 1 );, ngood=length(w);
                if ngood ~= 0
                    %ww=ww
                    resx( ispan,w(ww)) = eyes( which1( w(ww) ));
                    found( which1( w(ww)) ) = 1;
%if j==1,found,end                    
                    rowtot( w(ww) ) = 0;
                    coltot( which1(w(ww)) ) = 0;
                end
            end

            labely = find( rowtot >0);, ngood=length(labely);
            if ngood ~= 0
                labelx = find( coltot > 0 );
                nontrivial = 1;
            else
                nontrivial = 0;
            end
            clear abi,clear abpos,clear fnsh, clear rowtot, clear coltot, clear which1, clear isort
        else
            %   or: Use simple N^2 time routine to calculate trivial bonds

            % let's try a nice, loopless way!
            % don't bother tracking perm. lost guys.
            wh = find( pos(:,1) > 0)';, ntrack=length( wh );
            if ntrack == 0, warning('Warning - No valid particles to track!'), end
            xmat=ones(ntrack,m);
            for ng=1:m, xmat(:,ng)=ng;, end
            ymat=ones(m,ntrack);
            for ng=1:ntrack, ymat(:,ng)=ng;, end
            ymat=ymat';

            for d=1:dim
                x = (xyi(:,d))';
                y = (pos(wh,d));
                if d == 1, dq = (x(xmat)-(y(ymat))).^2;, else, dq = dq + (x(xmat)-(y(ymat))).^2;, end
            end
            ltmax = ( dq < maxdisq );

            % figure out which trivial bonds go with which
            rowtot = zeros(1,n);
            %wh=wh
            rowtot(wh) = sum( ltmax , 2)';
            if ntrack > 1, coltot = sum( ltmax, 1 );, else coltot = ltmax;, end        %STRANGE, LTMAX MIGHT BE A VECTOR OR A MATRIX
            which1 = zeros( 1,n );
            for j = 1:ntrack
                mx = max( ltmax(j,:) );   % max is faster than where
                w = find( ltmax(j,:)==mx );
                if length(w)>1, w=w(1);,end
                %ltmax=ltmax
                %TESTltmax=ltmax(j,:)
                %w=w,test=which1(wh(j))
                which1(wh(j)) = w;
            end

            ntrk = floor( n - sum(rowtot == 0) );
            w = find( rowtot == 1 );, ngood=length(w);
            if ngood ~= 0
                ww = find( coltot( which1(w) ) == 1);, ngood=length(ww);
                if ngood ~= 0
                    %eyes=eyes,w=w,ispan=ispan,ww=ww
                    resx( ispan,w(ww) ) = eyes( which1( w(ww) ));
                    found( which1( w(ww)) ) = 1;
                    %test2=rowtot( w(ww) );
                    rowtot( w(ww) ) = 0;
                    coltot( which1(w(ww)) ) = 0;
                end
            end
            labely = find( rowtot > 0 );, ngood=length(labely);
            if ngood ~= 0
                labelx = find( coltot > 0 );
                nontrivial = 1;
            else
                nontrivial = 0;
            end
        clear rowtot, clear coltot , clear which1 
        end

        %   THE TRIVIAL BOND CODE ENDS

        if nontrivial
            
            xdim = length( labelx );
            ydim = length( labely );

            %   make a list of the non-trivial bonds
            bonds = ones(1,2);
            bondlen = 0;
            for j=1:ydim
                distq = zeros( 1,xdim );
                for d=1:dim
                    %d=d,distq=distq,xyi=xyi,pos=pos
                    %test1=xyi(labelx,d),test2=pos(labely(j),d)
                    %test3=xyi(labelx,d) - pos(labely(j),d)
                    %test4=( xyi(labelx,d) - pos(labely(j),d) ).^2
                    distq = distq + ( xyi(labelx,d) - pos(labely(j),d) )'.^2 ;
                end
                w = find( distq < maxdisq );, ngood=length(w);
                bonds = [bonds ; w', zeros(ngood,1)+j];
                bondlen = [ bondlen, distq( w ) ];
            end
            bonds = bonds(2:end,:);
            bondlen = bondlen(2:end);
            numbonds = size( bonds,1 );
            mbonds = bonds;

            if max([xdim,ydim]) < 4
                nclust = 1;
                maxsz = 0;
                mxsz = xdim;
                mysz = ydim;
                bmap = zeros(1,size( bonds,1 )) - 1;
            else
                %THE SUBNETWORK CODE BEGINS

                lista = ones( 1,numbonds );
                listb = ones( 1,numbonds );
                nclust = 0;
                maxsz = 0;
                thru = xdim;
                %thruF=thru

                %bonds1=bonds
                while thru ~= 0
                    %Wstart=w
                    %BONDS=bonds

                    %    the following code extracts connected sub-networks of the non-trivial
                    %    bonds.  NB: lista/b can have redundant entries due to
                    %    multiple-connected subnetworks.
                    %bonds=bonds
                    w = find( bonds(:,2) > 0 );
                    %bonds=bonds,w=w%,w(1)=w(1)%,test=( bonds( w(1),2 ) )
                    lista(1) = ( bonds( w(1),2 ) );
                    %listaA=lista
                    listb(1) = ( bonds( w(1),1 ) );
                    %bondsB=bonds
                    %nclust=nclust
                    bonds(w(1),:) = -(nclust+1);
                    %bondsA=bonds
                    adda  = 2; addb  = 2;
                    donea = 1; doneb = 1;

                    repeat=1;
                    while repeat==1
                        %BONDStart=bonds
                        %donea=donea,adda=adda
                        if donea ~= adda
                            %bonds=bonds
                            w = find( bonds(:,2) == lista(donea));, ngood=length(w);
                            if ngood ~= 0
                                %bonds=bonds,w=w,listb_B=listb
                                listb(addb:addb+ngood-1) = bonds(w,1);
                                %listb_A=listb
                                %w=w
                                bonds(w,:) = -(nclust+1);
                                addb=addb+ngood;
                            end
                            donea=donea+1;
                        end
                        if doneb ~= addb
                            w = find( bonds(:,1) == listb(doneb));, ngood=length(w);
                            if ngood ~= 0
                                lista(adda:adda+ngood-1) = bonds(w,2);
                                %nclust=nclust
                                bonds(w,:) = -(nclust+1);
                                adda=adda+ngood;
                            end
                            doneb=doneb+1;
                        end
                        if ((donea == adda) & (doneb == addb)), repeat=0;, end
                    end
                    % a thing of beauty is a joy forever.

                    %lista=lista,listb=listb,donea=donea,doneb=doneb
                    %testBList=listb( 1:doneb - 1 )
                    %testBSort=sort(listb)%( 1:doneb))
                    %testAList=lista( 1:donea - 1 )
                    %testASort=sort(lista)%( 1:donea))
                    [tempb,idxsortb]=sort(listb( 1:doneb - 1 ));
                    [tempa,idxsorta]=sort(lista( 1:donea - 1));
                    xsz = length( unq( listb( 1:doneb - 1 ), idxsortb ) );
                    ysz = length( unq( lista( 1:donea - 1 ), idxsorta ) );

                    if xsz*ysz > maxsz
                        maxsz = xsz*ysz;
                        mxsz = xsz;
                        mysz = ysz;
                    end

                    %xsz=xsz
                    thru = thru - xsz;
                    nclust = nclust + 1;
                end
                %bonds=bonds
                bmap = ( bonds(:,1) )';
            end

            %   THE SUBNETWORK CODE ENDS

            % % %               SKIPPED THIS PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            % % %             if verbose then begin
            % % %             if nclust GT 1 then message,strcompress(i)+': '+'Permuting'+$
            % % %                strcompress(nclust)+' networks',/inf $
            % % %                else message,strcompress(i)+': '+'Permuting'+ $
            % % %                   strcompress(nclust)+' network',/inf
            % % %             message,'      Max. network:'+strcompress(mxsz)+' x'+$
            % % %                strcompress(mysz),/inf
            % % %             if keyword_set( add ) then message,'      Tracking'+$
            % % %                strcompress(n)+' particles.',/inf
            % % %          endif

            %   THE PERMUTATION CODE BEGINS
            for nc = 0:nclust-1
               %tocp=toc
                %bmap=bmap,nc=nc
                w = find( bmap == -(nc+1));, nbonds=length(w);
                %mbonds=mbonds,w=w
                bonds = mbonds( w,: );
                lensq = bondlen( w );
                              
                %test1=bonds(:,1)
                %test2=sort( bonds(:,1))
                [temp,sortidx]=sort( bonds(:,1));            
                uold = bonds( unq (bonds(:,1)', sortidx') , 1 );
                nold = length( uold );
                %bonds=bonds
                %test=unq(bonds(:,2)',[])
                unew = bonds( unq( bonds(:,2)',[] ) , 2 );
                nnew = length( unew );

                % check that runtime is not excessive
                if nnew > 5
                    rnsteps = 1;
                    for ii = 1:nnew
                        rnsteps = rnsteps * length( find( bonds(:,2) == unew(ii) ) );
                        %I comment the following command Gao
                        if rnsteps > 5e+4               %5e+4          I
                            %whos
                            warning( ' Warning: difficult combinatorics encountered.')
                            warning( ' Program may not finish- Try reducing maxdisp.')
                        end
                        if rnsteps > 2e+6, warning(' Excessive Combinatorics! Try reducing maxdisp.'), end
                    end
                end
                %nnew=nnew
                st = ones( 1,nnew );, fi = ones( 1,nnew );
                h = ones( 1,nbonds );, ok = ones( 1,nold )+1;
                %nnew=nnew,nold=nold
                if nnew-nold > 0
                   nlost = nnew - nold;
                else
                   nlost=0;
                end

                for ii=1:nold, h( find( bonds(:,1) == uold(ii) ) ) = ii; end

                %fi=fi,nbonds=nbonds
                %nnew=nnew
                st(1) = 1;, fi(nnew) = nbonds;
                %nnew=nnew
                if nnew > 1
                    %nnew=nnew
                    sb = ( bonds(:,2) )';
                    sbr = circshift( sb, [0,1]);
                    sbl = circshift( sb, [0,-1]);
                    %st=st,sb=sb,sbr=sbr
                    st(2:end) = find( sb(2:end) ~= sbr(2:end) ) +1;
                    fi(1:nnew-1) = find( sb(1:nbonds-1) ~= sbl(1:nbonds-1) );
                end

                checkflag = 0;
                while checkflag~=2
                    %st=st
                    pt = st - 1;
                    %pt=pt,st=st,fi=fi
                    lost = zeros( 1,nnew );
                    who = 1; losttot = 0;
                    mndisq = nnew*maxdisq;

                    while who~=0
                    %whowhile=who    
                        %lost=lost
                        %,losttot=losttot,nlost=nlost
                        %if i==13,who1=who,end
                        %who=who
                        %pt=pt
                        if pt(who) ~= fi(who)
                            %h=h
                            %test1=pt(who)+1
                            %test2=fi(who)
                            %okB=ok
                            %test=ok(h(pt(who)+1:fi(who)))
                            %ok=ok,h=h,pt=pt,who=who,fi=fi
                            w = find( ok(h(pt(who)+1:fi(who))));, ngood=length(w);
                            if ngood > 0
                                %who=who,test1=pt(who),test2=st(who),test3=fi(who)
                                %w=w
                                if pt(who) ~= st(who)-1;, ok(h(pt(who))) = 1;, end
                                %pt=pt
                                pt(who) = pt(who) + w(1);
                                %who=who,pt=pt,h=h,ok=ok
                                ok(h(pt(who))) = 0;
                                %test=h(pt(who))
                                %okA=ok
                                %if i==13,whoN=who,nnew=nnew,end
                                if who == nnew
                                    %ww=ww,lost=lost                                    % place #1 calc. tot. sqr. disp.
                                    ww = find( lost == 0 );
                                    %pt=pt,ww=ww
                                    %lensq=lensq
                                    %pt=pt,ww=ww
                                    %losttot=losttot
                                    dsq = sum( lensq(pt(ww)) ) + losttot*maxdisq;
                                    
                                    %mndisq,dsq
                                    %if i==13,dsq,mndisq,end
                                    %if (dsq - mndisq) < -1e-10
                                    if dsq < mndisq
                                        minbonds = pt( ww );
                                        mndisq = dsq;
                                    end
                                else
                                    who = who +1;
                                    %if i==13,who2=who,end
                                    %who2=who
                                end
                            else
                                notlost=-lost(who)-1;
                                %losttot,nlost
                                if mod(notlost,2)==1  & (losttot ~= nlost)             %PAY ATTENTION: used to be not lost (who) in idl
                                    lost(who) = 1;
                                    losttot = losttot +1;
                                    if pt(who) ~= st(who)-1, ok(h(pt(who))) = 1;, end
                                    %who,nnew
                                    if who == nnew
                                        % place #2 calc. tot. sqr. disp.
                                        ww = find( lost == 0 );
                                        dsq = sum( lensq(pt(ww)) )+ losttot*maxdisq;
                                        
                                        %mndisq,dsq
                                        if dsq < mndisq
                                            minbonds = pt( ww );
                                            mndisq = dsq;
                                        end
                                    else
                                        who = who +1;
                                        %if i==13,who3=who,end
                                        %who3=who

                                    end

                                    % Fight the 'c' power- long live IDL!
                                else
                                    if pt(who) ~= st(who) -1, ok(h(pt(who))) = 1;, end
                                    pt(who) = st(who) -1;
                                    if lost(who)
                                        lost(who) = 0;
                                        losttot = losttot -1;
                                    end
                                    who = who -1;
                                    %if i==13,who4=who,end
                                    %who4=who
                                end
                            end
                        else
                            %losttot=losttot,nlost=nlost
                            notlost=-lost(who)-1;
                            if mod(notlost,2)==1 & (losttot ~= nlost)                %PAY ATTENTION: used to be not lost (who) in idl
                                lost(who) = 1;
                                losttot = losttot +1;
                                if pt(who) ~= st(who)-1, ok(h(pt(who))) = 1;, end
                                if who == nnew
                                    % place #3 calc. tot. sqr. disp.
                                    ww = find( lost == 0 );
                                    dsq = sum( lensq(pt(ww)) ) + losttot*maxdisq;

                                    if dsq < mndisq
                                        minbonds = pt( ww );
                                        mndisq = dsq;
                                    end
                                else
                                    who = who +1;
                                    %if i==13,who5=who,end
                                    %who5=who
                                end
                            else
                                if pt(who) ~= st(who)-1, ok(h(pt(who))) = 1;, end
                                pt(who) = st(who) -1;
                                if lost(who)
                                    lost(who) = 0;
                                    losttot = losttot -1;
                                end
                                who = who -1;
                                %if i==13,who6=who,end
                                %who6=who
                            end
                        end
                        %whoE=who
                    end

                    checkflag = checkflag +1;
                    if checkflag == 1
                        %   we need to check that our constraint on nlost is not forcing us away from the minimum id's
                        plost = min([ floor( mndisq/maxdisq ) , nnew-1 ]);
                        if plost > nlost+1, nlost = plost;, else checkflag = 2;, end
                    end
                end

                %   update resx using the minimum bond configuration
                %ispan=ispan,labely=labely,bonds=bonds,minbonds=minbonds,eyes=eyes,labelx
                %(bonds(minbonds,1))
                resx(ispan,labely(bonds(minbonds,2))) = eyes(labelx(bonds(minbonds,1)));
                %found=found,labelx=labelx,bonds=bonds,minbonds=minbonds
                found( labelx(bonds(minbonds,1)) ) = 1;
            end
            %   THE PERMUTATION CODE ENDS
            %found=found
        else
            %if verbose, warning([num2str(i),': Only trivial networks']), end
        end

        %     here we want to update our initial position estimates
        %resx=resx%ispan=ispan;
        w = find( resx(ispan,:) > 0);, nww=length(w);
        if nww > 0
            %resx=resx
            pos(w,:) = xyzs(resx(ispan,w),1:dim);
            if goodenough > 0, nvalid(w) = nvalid(w) + 1; end
        else
            warning('Warning, tracking zero particles!')
        end

        %we need to add new guys, as appropriate.
        newguys = find( found == 0 );, nnew=length(newguys);

        if (nnew > 0)
            newarr = zeros( zspan , nnew ) - 1;
            resx = [resx,newarr];
            %resx=resx,
            %ispan=ispan,n=n,eyes=eyes,newguys=newguys

            resx(ispan,n+1:end) = eyes(newguys);
            pos = [[pos];[xyzs(eyes(newguys),1:dim)]];
            mem = [mem,zeros(1,nnew)];
            %maxid=maxid
            uniqid = [uniqid,[1:nnew]+maxid];
            maxid = maxid + nnew;
            if goodenough > 0
               dumphash = [dumphash,zeros(1,nnew)];
               nvalid = [nvalid,ones(1,nnew)];
            end
            n = n + nnew;
        end

    else
        Warning([' Warning- No positions found for t=',num2str(i)])
    end

    %   update the 'memory' array
    w = find( resx(ispan,:) ~= -1 ); nok=length(w);
    if nok ~= 0, mem( w ) = 0;, end         % guys get reset if they're found
    mem = mem + ( resx(ispan,:) == -1 );

    %  if a guy has been lost for more than memory times, mark him as permanently
    %  lost.  For now, set these guys to pos = ( -maxdisp, -maxdisp, ... ),
    %  so we can never track them again. It would be better to make a smaller
    %  pos, but then we'd have to change 'n', which would be gnarly.
    wlost = find( mem == memory+1 );, nlost=length(wlost);
    if nlost > 0
        pos(wlost,:) = - maxdisp;
        % check to see if we should 'dump' newly lost guys
        if goodenough > 0
           wdump = find( nvalid(wlost) < goodenough); ndump=length(wdump);
           if ndump > 0, dumphash( wlost(wdump) ) = 1; end
        end
    end

    %  we need to insert the working copy of resx into the big copy bigresx
    %  do our house keeping every zspan time steps (dumping bad lost guys)
    if (ispan == zspan) | (i == z)
        %  if a permanently lost guy has fewer than goodenough valid positions
        %  then we 'dump' him out of the data structure- this largely alleviates
        %  memory problems associated with the 'add' keyword and 'noise' particles
        %  To improve speed- do it infrequently.
        % in case we've added some we need to pad out bigresx too
        nold = length( bigresx(1,:) );
        n=n;
        nnew = n - nold;
        if nnew > 0
            newarr = zeros( z, nnew ) - 1;
            bigresx = [bigresx,newarr];
        end

        if goodenough > 0
           if (sum(dumphash) > 0)
        %         if verbose then message,'Dumping bad trajectories...',/inf
              wkeep = find( dumphash == 0); nkeep=length(wkeep);
              resx = resx(:,wkeep);
              bigresx = bigresx(:,wkeep);  % this really hurts runtime
              pos = pos(wkeep,:);
              mem = mem(wkeep);
              uniqid = uniqid(wkeep);
              nvalid = nvalid(wkeep);
              n = nkeep;
              dumphash = zeros(nkeep);
           end
        end


        %  if not verbose then message,strcompress(i+1)+' of'$
        %  +strcompress(z)+' done. Tracking'+strcompress(ntrk)$
        %  +' particles,'+strcompress(n)+' tracks total.',/inf
        %resx=resx
        %i,ispan
        bigresx(i-ispan+1:i,:) = resx(1:ispan,:);
        resx = zeros(zspan,n) - 1;

        %  We should pull permanently lost guys, parse them and concat them
        %  onto the 'output list', along with their 'unique id' number to
        %  make scanning the data files a little easier.  Do infrequently.
        wpull = find( pos(:,1) == -maxdisp ); npull=length(wpull);
        if npull > 0
            %   pos( 0, wpull ) = -2*maxdisp
            lillist = [1,1];
            for ipull = 1:npull;
                %wpull=wpull
                wpull2 = find( bigresx(:,wpull(ipull)) ~= -1 );, npull2=length(wpull2);
                %bigresx=bigresx,wpull=wpull,ipull=ipull,wpul2=wpull2
                %test1=bigresx(wpull(ipull),wpull2),test2=zeros(1,npull2-1),test3=uniqid(wpull(ipull)),npull2=npull2
                lillist = [[lillist];[bigresx(wpull2,wpull(ipull)), zeros(npull2,1)+uniqid(wpull(ipull))]];
                lillist=lillist;
            end
            %olist=olist,lillist=lillist
            olist = [[olist];[lillist(2:end,:)]];
        end

        %     now get rid of the guys we don't need anymore....
        %     but watch out for when we have no valid particles to track!
        wkeep = find( pos(:,1) > 0 )';, nkeep=length(wkeep);
        if nkeep == 0, Warning(' We are going to crash now, no particles....'), end
        resx = resx(:,wkeep);
        %bigresx=bigresx
        bigresx = bigresx(:,wkeep);  % this really hurts runtime
        pos = pos(wkeep,:);
        mem = mem(wkeep);
        uniqid = uniqid(wkeep);
        n = nkeep;
        dumphash = zeros(1,nkeep);
        if goodenough > 0, nvalid = nvalid(wkeep); end
    end

end     % the big loop over z time steps....

%   %  make a final scan for short trajectories that weren't lost at the end.
   if goodenough > 0
      nvalid = sum(bigresx > 0);
      wkeep = find( nvalid >= goodenough); nkeep = length(wkeep);
      if nkeep < n
         bigresx = bigresx(:,wkeep);
         n = nkeep;
         uniqid = uniqid(wkeep);
         pos = pos(wkeep,:);
      end
   end

%  make the final scan to 'pull' everybody else into the olist.

wpull = find( pos(:,1) ~= -2*maxdisp )';, npull = length(wpull);
if npull > 0
    lillist = [1,1];
    for ipull = 1:npull
        wpull2 = find(bigresx(:,wpull(ipull)) ~= -1); npull2=length(wpull2);
        %bigresx=bigresx %,uniqid=uniqid
        %bigresx=bigresx,uniqid=uniqid,wpull=wpull,ipull=ipull,wpull2=wpull2
        lillist = [[lillist];[bigresx(wpull2,wpull(ipull)),zeros(npull2,1)+uniqid(wpull(ipull))]];
        lillist=lillist;
        %olist=olist
    end
    olist = [[olist];[lillist(2:end,:)]];
end
olist = olist(2:end,:);

%  free up a little memory for the final step!
bigresx = 0;
resx = 0;

% need to make up a result array!

% if verbose then message,"Preparing result array...",/inf
nolist = length(olist(:,1));
res = zeros(nolist,dd+1);
%olist=olist
%xyzs=xyzs
%dd=dd

for j = 1:dd
    res(:,j) = xyzs(olist(:,1),j);
end
res(:,dd+1) = olist(:,2);

if size(res,1)>0
    lub=luberize(res);
else
    lub=res;
end

% ***************** end of track.pro

%whos


%done=toc,StartTime








