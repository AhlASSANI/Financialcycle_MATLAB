%----------------------------------------------------------------
%% OBJECT: Turning Point Indicator
%----------------------------------------------------------------
% DESCRIPTIVE:  [nbTP,TPIndex,nb_t] = TPI(Dat,Prob,Thresh,Hmax)
%----------------------------------------------------------------
% INPUTS: 
% Dat    : date of recessions (NBER)                      : T x 1 
% Prob   : estimated probabilities of recessions          : T x 1 
% Thresh : if Prob > Thresh, recession signal             : 1 x 1 
% Hmax   : allowed leads of lags                          : 1 x 1 
%----------------------------------------------------------------
% OUTPUTS: 
% nbTP  : nb of turning points                             : 1 x 1 
% TPIndex : ratio nb of estimated TP / nb of observed TP   : 1 x 1 
% nb_t : % of accurate TP detected with a delay of +/-Hmax : 1 x 1 
%----------------------------------------------------------------


function [nbTP,TPIndex,nb_t] = TPI(Dat,Prob,Thresh,Hmax)

T=size(Dat,1);
nbTP = sum((Dat(2:end)-Dat(1:end-1)).^2);
Dat_est = Prob>Thresh;
nbTP_est = sum((Dat_est(2:end)-Dat_est(1:end-1)).^2);
TPIndex = 100*nbTP_est/nbTP;
nb_t=0;
for t=2+Hmax:T-Hmax
   temp = (Dat_est(t-Hmax:t+Hmax)-Dat_est(t-Hmax-1:t+Hmax-1)) ;
   temp = temp*(Dat(t)-Dat(t-1)) ;
   temp = max(temp);
   nb_t = nb_t+temp;
end

nb_t = nb_t/nbTP;