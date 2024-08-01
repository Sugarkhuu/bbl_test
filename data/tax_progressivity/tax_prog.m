%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file constructs a tax progressivity measure for 1946-2017
% This is an extension on the work of:
% Karel Mertens and Jose Montiel-Olea, ``Marginal Tax Rates and Income, 
%                                    New Time Series Evidence''
% March, 2023 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 

% user_path   = "C:\Users\Jrxz12\Dropbox\Master_Thesis\1_Data\tax_progressivity"
% output_path = "C:\Users\Jrxz12\Dropbox\Master_Thesis"
% fullpath    = user_path + "\AMTR_construction\mertens_olea"
% cd("/data/tax_progressivity");

addpath('auxiliary_files');

%% Step 1  Load Input Data
data_DAGI            = xlsread('data/SOI_AGI_Distributions.xlsx','ALLYEARS'); % Statistics of Income, adjusted gross income
data_MR              = xlsread('data/SOI_Marginal_Rates_by_AGI_and_FS.xlsx','ALLYEARS');
TSERIES              = xlsread('data/TIME_SERIES_DATA.xlsx','SERIES');
TSERIES              = TSERIES(1:72,:);
YEARS                = (1946:2017)';
AMTR_old             = xlsread('data/AMTR_1946_2012_Mertens_Olea.xlsx', 2)

%% Step 2 Fit Distributions of Adjusted Gross Income to Available Data from
% the IRS Statistics of Income
FitDistributions  % This looks into data_DAGI
load D

%% Step 3 Construct AGI Floors for Income Percentiles
quan      = [0 0.99 0.95 0.90];
TAXUNITS  = TSERIES(:,2)*1000;
AGIFLOORS = fun_AGIFLOORS(quan,D,YEARS,TAXUNITS);

%% Step 4 Construct MIITRS using method 1
MIITRs_method1

%% Step 5 Construct the MIITR Series for the entire sample. 
load MIITRS_m1;
AMIITR     = [YEARS TSERIES(:,8:11)];                        % Saez (2004)

%% Step 6: Extend Saez (2004) Series By Regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sel   = (YEARS>1986)&(isnan(MIITRS_m1(:,2))==0)&(isnan(AMIITR(:,2))==0);
Y     = AMIITR(sel,2:end);
XX    = MIITRS_m1(:,2:end);
X     = XX(sel,:); 
b     = [ones(sum(sel),1) X]\Y;
Y2hat = [ones(length(XX),1) XX]*b; 
AMIITR((YEARS>1986)&(isnan(MIITRS_m1(:,2))==0)&(isnan(AMIITR(:,2))==1),2:end)= Y2hat((YEARS>1986)&(isnan(MIITRS_m1(:,2))==0)&(isnan(AMIITR(:,2))==1),:);

AMIITR(1:67, 2) = AMTR_old(:, 3)
%% Step 5: Get the AMPTR
data_SSA             = xlsread('data/SSA_rates_and_earnings.xlsx','ALLYEARS');
data_SSA             = data_SSA(1:72,:);
OASDI_rate_EMPLOYEES = data_SSA(:,4)/100;
OASDI_rate_EMPLOYERS = data_SSA(:,5)/100;
HI_rate_EMPL2        = data_SSA(:,6)/100;
OASDI_rate_SELFEMPL  = data_SSA(:,7)/100;
HI_rate_SELFEMPL     = data_SSA(:,8)/100;

OASDIWAGES  = data_SSA(:,9);  % Total Wages subject to OASDI tax of those with wages below the ceiling
HIWAGES     = data_SSA(:,10); % Total Wages subject to Medicare tax of those with wages below the ceiling
OASDISEE    = data_SSA(:,11); % Total Self Employment Earnings subject to OASDI tax of those with earnings below the ceiling
HISEE       = data_SSA(:,12); % Total Self Employment Earnings subject to Medicare tax of those with earnings below the ceiling
TOTINC      = TSERIES(:,3);

% (A) Weighted by MARKET INCOME (SAEZ Income Concept)
AMPTR = OASDIWAGES./TOTINC.*(OASDI_rate_EMPLOYEES+OASDI_rate_EMPLOYERS)./(1+OASDI_rate_EMPLOYERS+0.5*HI_rate_EMPL2)...
       + HIWAGES./TOTINC.* HI_rate_EMPL2./(1+OASDI_rate_EMPLOYERS+0.5*HI_rate_EMPL2)...
       + OASDISEE./TOTINC.*OASDI_rate_SELFEMPL+HISEE./TOTINC.*HI_rate_SELFEMPL;


%% Step 6: Get the AMTR and ATR
ATR   = xlsread('data/ATR.xlsx',1)
ATR   = ATR(34:end-2, 11)

AMTR = AMPTR + AMIITR(:, 2)
%% Step 7: Construct Progressivity Measure
P = (AMTR - ATR) ./ (1 - ATR)

%% Step 8: Save
csvwrite('../data_final/tax_progressivity.csv', [YEARS P])
