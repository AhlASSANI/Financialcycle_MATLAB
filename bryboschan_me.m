function [nber,peaks,troughs,index_peaks, index_troughs] = bryboschan_me(x,p)
% Input Values:
%   data = data series
%   p = lenghts of periods
%       1 = monthly (default)
%       2 = quarterly

% Calculates the Bry and Boschan (1971) procedure for programmed
% determination of turning points.

% Output Values:
%   index_peaks = index_peaks are marked with the value 1
%   index_troughs = index_troughs are marked with the value 1

% Copyright by Martin P. Everts, University of Bern, 2005
% Modified by Lola Gadea, 2011

if p==not(or(1,2));
    error('Wrong input value for "p".');
end;

p3=3;
p4=4;
p5=5;
p6=6;
p12=12;
p15=15;

if p==2;
    p3=round(p3/12*4);
    p4=round(p4/12*4);
    p5=round(p5/12*4);
    p6=round(p6/12*4);
    p12=round(p12/12*4);
    p15=round(p15/12*4);
end;
Rx=size(x,1);

% Step I: DETERMINATION OF EXTREMES AND SUBSTITUTION OF VALUES.

xsp=spencer_me(x,p15);
xo=bbout_me(x,xsp);


% Step II: DETERMINATION OF CYCLES IN P12-PERIOD MOVING AVERAGE (EXTREMES
% REPLACED).

if p==1;
    xma=maqe_me(xo,p12);
else
    xma=maqe_me(xo,(p12/2));
end;


% Substep II A: IDENTIFICATION OF POINTS HIGHER (OR LOWER) THAN P5
% PERIODS ON EITHER SIDE.

[bcp2a,bct2a]=dates_me(xma,p5);


% Substep II B: ENFORCEMENT OF ALTERNATION OF TURNS BY SELECTING HIGHEST OF
% MULTIPLE PEAKS (OR LOWEST OF MULTIPLE TROUGHS).

[bcp2b,bct2b]=alter_me(bcp2a,bct2a,x,xma);


% Step III: DETERMINATION OF CORRESPONDING TURNS IN SPENCER CURVE (EXTREMES
% REPLACED).

xspo=spencer_me(xo,p15);


% Substep III A: IDENTIFICATION OF HIGHEST (OR LOWEST) VALUE WITHIN P5
% PERIODS OF SELECTED TURN IN P12-PERIODS MOVING AVERAGE.

[bcp3a,bct3a]=refine_me(bcp2b,bct2b,xspo,p5);
[bcp3a,bct3a]=alter_me(bcp3a,bct3a,x,xspo);


% Substep III B: ENFORCEMENT OF MINIMUM CYCLE DURATION OF P15 PERIODS BY
% ELIMINATING LOWER PEAKS AND HIGHER TROUGHS OF SHORTER CYCLES.

[bcp3b,bct3b]=enf_me(bcp3a,bct3a,xspo,p15);


% Step IV: DETERMINATION OF CORRESPONDING TURNS IN SHORT-TERM MOVING
% AVERAGE OF P3 TO P6 PERIODS, DEPENDING ON PCD (PERIODS OF CYCLICAL
% DOMINANCE).

xsp=spencer_me(x,p15);

if p==1;
    pcd=pcd_me(x,xsp,p12*2);
elseif p==2;
    pcd=pcd_me(x,xsp,p12);
end;

if p>1;
    if pcd<=1;
        xpcd=x;
    else
        xpcd=maqe_me(x,2);
    end;
elseif p==1;
    if pcd<=2;
        xpcd=maqo_me(x,3);
    elseif pcd<=4;
        xpcd=maqe_me(x,4);
    elseif pcd<=6;
        xpcd=maqo_me(x,5);
    else
        xpcd=maqe_me(x,6);
    end;
end;


% Substep IV A: IDENTIFICATION OF HIGHEST (OR LOWEST) VALUE WITHIN P5
% PERIODS OF SELECTED TURN IN SPENCER CURVE.

[bcp4a,bct4a]=refine_me(bcp3b,bct3b,xpcd,p5);
[bcp4a,bct4a]=alter_me(bcp4a,bct4a,x,xpcd);


% Step V: DETERMINATION OF TURNING POINTS IN UNSMOOTHED SERIES.

% Substep V A: IDENTIFICATION OF HIGHEST (OR LOWEST) VALUE WITHIN P4
% PERIODS, OR PCD TERM, WHICHEVER IS LARGER, OF SELECTED TURN IN SHORT-TERM
% MOVING AVERAGE.

span=max(pcd,p4);
[bcp5a,bct5a]=refine_me(bcp4a,bct4a,xpcd,span);
[bcp5a,bct5a]=alter_me(bcp5a,bct5a,x,xpcd);


% Substep V B: ELIMINATION OF TURNS WITHIN P6 PERIODS OF BEGINNING AND END
% OF SERIES.

[bcp5b,bct5b]=enfvb_me(bcp5a,bct5a,x,p6);
[bcp5b,bct5b]=alter_me(bcp5b,bct5b,x);


% Substep V C: ELIMINATION OF PEAKS (OR TROUGHS) AT BOTH ENDS OF SERIES
% WITH ARE LOWER (OR HIGHER) THAN VALUES CLOSER TO END.

[bcp5c,bct5c]=enfvc_me(bcp5b,bct5b,x);
[bcp5c,bct5c]=alter_me(bcp5c,bct5c,x);


% Substep V D: ELIMINATION OF CYCLES WHOSE DURATION IS LESS THAN P15
% PERIODS.

[bcp5d,bct5d]=enf_me(bcp5c,bct5c,x,p15);
[bcp5d,bct5d]=alter_me(bcp5d,bct5d,x);


% Substep V E: ELIMINATION OF PHASES WHOSE DURATION IS LESS THAN P5
% PERIODS.

[bcp5e,bct5e]=enfve_me(bcp5d,bct5d,x,p5);
[bcp5e,bct5e]=alter_me(bcp5e,bct5e,x);


% Step VI: STATEMENT OF FINAL TURNING POINTS.

if Rx==size(bcp5e,1) && Rx==size(bct5e,1)
    index_peaks=bcp5e;
    index_troughs=bct5e;
else
    for n=1:j
        index_peaks(n,1)=NaN;
        index_troughs(n,1)=NaN;
    end;
    index_peaks(j+1:Rx,1)=bcp5e;
    index_troughs(j+1:Rx,1)=bct5e;
end;
peaks=find(index_peaks);
troughs=find(index_troughs);

if isempty(peaks)==1 && isempty(troughs)==1
    nber=NaN(Rx,1);
elseif isempty(peaks)==0 && isempty(troughs)==1
    nber=[zeros(peaks,1);ones(Rx-peaks,1)];
elseif isempty(troughs)==0 && isempty(peaks)==1
    nber=[ones(troughs,1);zeros(Rx-troughs,1)];
else    
s=states(peaks,troughs,length(x));
nber=1-(s+1)/2;
end


%==========================================================================
% FUNCTIONS
%==========================================================================

function xsp = spencer_me(x,p)

% Input Values:
%   x = data series
%   p = amount of terms in Spencer forumla

% Filters a data series useing a p-term two sided Spencer filter.  See
% Kendall and Stuart (1966).

% Output Values:
%   xsp = filtered data series

if p==15;
    s=[-3 -6 -5 3 21 46 67 74 67 46 21 3 -5 -6 -3];
elseif p==5;
    s=[-3 12 17 12 -3];
end;

s=s/sum(s);

Rx=size(x,1);

xpad=[x(1,1)*ones((p-1)/2,1);x;x(Rx,1)*ones((p-1)/2,1)];

xsp=zeros(Rx,1);

for i=1:Rx;
    xsp(i,1)=s*xpad(i:i+p-1,1);
end;


%==========================================================================

function xo = bbout_me(x,xsp)

% Input Values:
%   x = data series
%   xsp = Spencer filtered data series

% Determinates extremes and substitutes the values.

% Output Values:
%   xo = data series without outliers

d=x-xsp;
sd=std(d);
xm=mean(d);

ds=(d-xm)./sd;
dsi=abs(ds) >= 3.5;

xo=x;

if sum(dsi) > 0;
    Ro = size(dsi,1);
    for i=1:Ro;
        if dsi(i,1)==1;
            xo(i,1)=xsp(i,1);
        elseif dsi(i,1)>1;
            error('There is an extreme outlier!');
        end;
    end;
    %     disp('Extreme values were determined and substituted.')
end;


%==========================================================================

function xf = maqe_me(x,q)

% Input Values:
%   x = data series
%   q = order of moving-average (only even data)

% Calculates a symmetric centered moving-average of lenghts q.  This
% calculation only works for even values of q.  For odd values please use
% function maqo_me(x,q).

% Output Values:
%   xf = filtered data series

if rem(q,2)==1;
    error('Value q is not even! Change q or use function maqo_me(x,q).');
end;

hq=q/2;

w=ones(q,1);
s=[0;w]+[w;0];

s=s/sum(s);

Rx=size(x,1);

xpad=[x(1,1)*ones(hq,1);x;x(Rx,1)*ones(hq,1)];

xf=zeros(Rx,1);

for i=1:Rx;
    xf(i,1)=ctranspose(s)*xpad(i:i+q,1);
end;


%==========================================================================

function xf = maqo_me(x,q)

% Input Values:
%   x = data series
%   q = order of moving-average (only odd data)

% Calculates a symmetric centered moving-average of lenghts q.  This
% calculation only works for odd values of q.  For even values please use
% function maqe_me(x,q).

% Output Values:
%   xf = filtered data series

if rem(q,2)==0;
    error('Value q is not odd! Change q or use function maqe_me(x,q).');
end;

hq=(q-1)/2;

w=ones(q,1);

w=w/sum(w);

Rx=size(x,1);

xpad=[x(1,1)*ones(hq,1);x;x(Rx,1)*ones(hq,1)];

xf=zeros(Rx,1);

for i=1:Rx;
    xf(i,1)=ctranspose(w)*xpad(i:i+q-1,1);
end;


%==========================================================================

function [bcp,bct] = dates_me(x,q)

% Input Values:
%   x = data series
%   q = length of periods

% Identifies points that are higher (or lower) than q periods on either side.

% Output Values:
%   bcp = business cycle peaks
%   bct = business cycle troughs

Rx=size(x,1);

bcp=zeros(Rx,1);
bct=zeros(Rx,1);

np=0;
nt=0;
st=0;

for i=q+1:Rx-q;
    if x(i,1)==max(x(i-q:i+q,1));
        if x(i,1)==x(i-1,1);
            st=1;
        else
            bcp(i,1)=1;
            np=1;
        end;
    elseif x(i,1)==min(x(i-q:i+q,1));
        if x(i,1)==x(i+1,1);
            st=1;
        else
            bct(i,1)=1;
            nt=1;
        end;
    end;
end;

if st==1;
    disp('DATA SERIES CONTAINS A STEP!!!');
end;
if np==0;
    disp('There are no peaks in this series.');
end;
if nt==0;
    disp('There are no troughs in this series.');
end;


%==========================================================================

function [bcpn,bctn] = alter_me(bcp,bct,x,xmod);

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs
%   xmod = modified data series

% Enforces alternation of turns by selecting highest of multiple peaks (or
% lowest of multiple troughs).

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)

if nargin==3;
    xmod=x;
end;

Rx=size(x,1);

pflag=0;
tflag=0;

bcpn=bcp;
bctn=bct;

for j=1:Rx;
    if and(bcp(j,1),bct(j,1))==1;
        disp('PEAK AND TROUGH AT THE SAME TIME!!!  BOTH VALUES WERE DELETED!!!')
        bcpn(j,1)=0;
        bctn(j,1)=0;
    elseif bcp(j,1)==1;
        if pflag==0;
            pflag=1;
            tflag=0;
            ps=j;
        elseif pflag==1;
            if x(j,1)>x(ps,1);
                bcpn(ps,1)=0;
                ps=j;
            elseif x(j,1)<x(ps,1);
                bcpn(j,1)=0;
            elseif xmod(j,1)>xmod(ps,1);
                bcpn(ps,1)=0;
                ps=j;
            elseif xmod(j,1)<xmod(ps,1);
                bcpn(j,1)=0;
            else;
                bcpn(j,1)=0;
                disp('Two peaks at same hight! Chose first peak!');
            end;
        end;
    elseif bct(j,1)==1;
        if tflag==0;
            tflag=1;
            pflag=0;
            ts=j;
        elseif tflag==1;
            if x(j,1)<x(ts,1);
                bctn(ts,1)=0;
                ts=j;
            elseif x(j,1)>x(ts,1);
                bctn(j,1)=0;
            elseif xmod(j,1)<xmod(ts,1);
                bctn(ts,1)=0;
                ts=j;
            elseif xmod(j,1)>xmod(ts,1);
                bctn(j,1)=0;
            else
                bctn(ts,1)=0;
                ts=j;
                disp('Two troughs at same hight! Chose last trough!');
            end;
        end;
    end;
end;


%==========================================================================

function [bcpn,bctn] = refine_me(bcp,bct,x,q);

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs
%   q = lenght of periods

% Identifies the highest (or lowest) value within q periods of selected
% turn.

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)

pd=find(bcp(:,1));
td=find(bct(:,1));

Rx=size(x,1);

bcpn=zeros(Rx,1);
bctn=zeros(Rx,1);

for i=1:size(pd,1);
    ps=max(1,(pd(i,1)-q));
    pe=min(Rx,(pd(i,1)+q));
    for j=ps:pe;
        if x(j,1)==max(x(ps:pe,1));
            if sum(bcpn(ps:pe,1))==0;
                bcpn(j,1)=1;
            else
                %                 disp('There was a step around the peak.')
            end;
        end;
    end;
    clear('ps','pe');
end;

for m=1:size(td,1);
    ts=max(1,(td(m,1)-q));
    te=min(Rx,(td(m,1)+q));
    for n=ts:te;
        if x(te+ts-n,1)==min(x(ts:te,1));
            if sum(bctn(ts:te,1))==0;
                bctn(te+ts-n,1)=1;
            else
                %                 disp('There was a step around the trough.')
            end;
        end;
    end;
    clear('ts','te');
end;


%==========================================================================

function [bcpn,bctn] = enf_me(bcp,bct,x,q)

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs
%   q = lenghts of period

% Enforces a minimum cycle duration of q periods by eliminating lower
% peaks and higher troughs of shorter cycles.

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)


pd=find(bcp(:,1));
td=find(bct(:,1));

Rx=size(x,1);

bcpn=bcp;
bctn=bct;

for i=1:size(pd,1)-1;
    if (pd(i+1,1)-pd(i,1))<q;
        if x(pd(i+1,1),1)>x(pd(i,1),1);
            bcpn(pd(i,1),1)=0;
        elseif x(pd(i+1,1),1)<x(pd(i,1),1);
            bcpn(pd(i+1,1),1)=0;
        else
            bcpn(pd(i,1),1)=0;
            warning('TWO PEAKS OF SAME HIGHT!');
        end;
    end;
end;

for i=1:size(td,1)-1;
    if (td(i+1,1)-td(i,1))<q;
        if x(td(i+1,1),1)>x(td(i,1),1);
            bctn(td(i+1,1),1)=0;
        elseif x(td(i+1,1),1)<x(td(i,1),1);
            bctn(td(i,1),1)=0;
        else
            bctn(td(i+1,1),1)=0;
            disp('TWO TROUGHS OF SAME HIGHT!');
        end;
    end;
end;


%==========================================================================

function pcd = pcd_me(x,xsp,h)

% Input Values:
%   x = data series
%   xsp = spencer curve of data series
%   h = limit (choose double the size of larges alternative moving average)

% Calculated the periods of cyclical dominance of a series.

% Output Values:
%   pcd = period of cyclical dominance

xo=x-xsp;

pcdv=zeros(h,1);
pcd=0;
k=1;

Rx=size(x,1);

for j=1:h;
    num=sum(abs(xo(1+j:Rx,1)-xo(1:Rx-j,1)));
    den=sum(abs(xsp(1+j:Rx,1)-xsp(1:Rx-j,1)));
    pcdv(j,1)=num/den;
end;

while pcd==0;
    if pcdv(k,1)<1;
        pcd=k;
    end;
    if k==h;
        error('MCD beyond given h-limit. Increase h!');
    end;
    k=k+1;
end;


%==========================================================================

function [bcpn,bctn] = enfvb_me(bcp,bct,x,q)

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs
%   q = lenghts of periods

% Eliminates turns within q periods of beginning and end of series.

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)

pd=find(bcp(:,1));
td=find(bct(:,1));

Rx=size(x,1);

bcpn=bcp;
bctn=bct;
nv=0;

for i=1:size(pd,1);
    if pd(i,1)<=q;
        bcpn(pd(i,1),1)=0;
        nv=1;
    end;
    if pd(i,1)>=Rx-q;
        bcpn(pd(i,1),1)=0;
        nv=1;
    end;
end;

for i=1:size(td,1);
    if td(i,1)<=q;
        bctn(td(i,1),1)=0;
        nv=1;
    end;
    if td(i,1)>=Rx-q;
        bctn(td(i,1),1)=0;
        nv=1;
    end;
end;

% if nv==1;
%     disp('A violation of proximity to endpoints has been found.');
% end;


%==========================================================================

function [bcpn,bctn] = enfvc_me(bcp,bct,x);

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs

% Eliminates peaks (or troughs) at both ends of series which are lower (or
% higher) than values closer to end.

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)

pd=find(bcp(:,1));
td=find(bct(:,1));

bcpn=bcp;
bctn=bct;
nv=0;

Rx=size(x,1);
Rp=size(pd,1);
Rt=size(td,1);

if and(Rp,Rt)>0;
    if pd(1,1)<td(1,1);
        if x(pd(1,1),1)<max(x(1:pd(1,1)-1,1));
            bcpn(pd(1,1),1)=0;
            nv=1;
        end;
    elseif pd(1,1)>td(1,1);
        if x(td(1,1),1)>min(x(1:td(1,1)-1,1));
            bctn(td(1,1),1)=0;
            nv=1;
        end;
    end;
    if pd(Rp,1)>td(Rt,1);
        if x(pd(Rp,1),1)<max(x(pd(Rp,1)+1:Rx,1));
            bcpn(pd(Rp,1),1)=0;
            nv=1;
        end;
    elseif pd(Rp,1)<td(Rt,1);
        if x(td(Rt,1),1)>min(x(td(Rt,1)+1:Rx,1));
            bctn(td(Rt,1),1)=0;
            nv=1;
        end;
    end;
end;

% if nv==1;
%     disp('Eliminated lower (higher) peaks (troughs) close to endpoints.');
% end;


%==========================================================================

function [bcpn,bctn] = enfve_me(bcp,bct,x,q)

% Input Values:
%   bcp = business cycle peaks (only 0 or 1)
%   bct = business cycle troughs (only 0 or 1)
%   x = data series underlying the peaks and troughs
%   q = lenghts of periods

% Enforces a minimum phase duration of q periods by eliminating peaks
% and troughs of shorter cycles.

% Output Values:
%   bcpn = business cycle peaks (only 0 or 1)
%   bctn = business cycle troughs (only 0 or 1)

pd=find(bcp(:,1));
td=find(bct(:,1));

Rx=size(x,1);
Rp=size(pd,1);
Rt=size(td,1);

bcpn=bcp;
bctn=bct;
nv=0;

if abs(Rt-Rp)>1;
    error('The difference between the amount of peaks and troughs is too high!');
end;

for i=1:min(Rp,Rt);
    if abs(pd(i,1)-td(i,1))<q;
        bcpn(pd(i,1),1)=0;
        bctn(td(i,1),1)=0;
        nv=1;
    end;
end;

if and(Rp,Rt)>0;
    if pd(1,1)<td(1,1);
        for i=1:Rp-1;
            if abs(pd(i+1,1)-td(i,1))<q;
                bcpn(pd(i+1,1),1)=0;
                bctn(td(i,1),1)=0;
                nv=1;
            end;
        end;
    elseif pd(1,1)>td(1,1);
        for i=1:Rt-1;
            if abs(pd(i,1)-td(i+1,1))<q;
                bcpn(pd(i,1),1)=0;
                bctn(td(i+1,1))=0;
                nv=1;
            end;
        end;
    end;
end;

% if nv==1;
%     disp('A violation of minimum phase length has been found');
% end;

%==========================================================================
function[state] = states(bcp,bct,nd)
%1=expansion;-1=contraction
state=zeros(nd,1);
if (ismiss(bcp) && ismiss(bct));
   return
elseif ismiss(bcp);
   state(bct(1)+1:nd)=ones(nd-bct(1),1);
   return
elseif ismiss(bct);
    state(bcp(1)+1:nd)=-ones(nd-bcp(1),1);
    return
else
    if bcp(1) < bct(1);  %index_peaks first
        state(1:bcp(1))=ones(bcp(1),1);
    else                %index_troughs first
        state(1:bct(1))=-1*ones(bct(1),1);
    end;
   
    amat=[[bcp,ones(size(bcp,1),1)] ; [bct,-1*ones(size(bct,1),1)]];
    amat=sortrows(amat,1);
    if amat(size(amat,1),1)<nd;
        amat=[amat;[nd,0]];
        end;

    i=1;
   while i<size(amat,1);
      state(amat(i,1)+1:amat(i+1,1))=-1*amat(i,2)*ones(amat(i+1,1)-amat(i,1),1);
      i=i+1;
    end;
    end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function [D]=ismiss(C)
 [a,b]=size(C);
 D=0;
 for i=1:a
     for j=1:b
         if isempty(C(i,j))==1 D=1;
         end
     end
 end
%==========================================================================
