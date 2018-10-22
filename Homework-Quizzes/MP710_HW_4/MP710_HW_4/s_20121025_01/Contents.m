% Quantitative MRI Analysis Package
% Version 0.1a (QMAP)   10-Jan-2012
%__________________________________________________________________________
%     __    __  __          _____  
%   / __ \ |  \/  |   /\   |  __ \ 
%  | |  | || \  / |  /  \  | |__) |
%  | |  | || |\/| | / /\ \ |  ___/  Quantitative MRI Analysis Package
%  | |__| || |  | |/ ____ \| |      University of Wisconsin - Madison
%   \___\_\|_|  |_/_/    \_\_|      (C) 2005-2012 | v0.01a | Jan-2011
%__________________________________________________________________________
%
% QMAP is a collection of quantitative MRI analysis software. 
%
% This file contains a function list.
% Type 'doc qmap' for a full description and copyright notice.
%
%   qmap                 - Quantitative MRI Analysis Package
%
% Files
%   qmap                 - Quantitative MRI Analysis Package
%   qmap_B1_eff          - FUNCTION qmap_B1eff = B1_eff(mtflip, tm, pulse_type, fermipulse)
%   qmap_b1_fit_bs       - FUNCTION b1m = qmap_b1_fit_bs(bs_plus, bs_minus)
%   qmap_check_var       - FUNCTION res = check_var(name);
%   qmap_fa_fit_afi      - FUNCTION fam = qmap_fa_fit_afi(afi_0, afi_1, n)
%   qmap_progressbar     - FUNCTION []  = qmap_progressbar(fractiondone, position)
%   qmap_qMT_fit_spgr    - FUNCTION [fv rnrm] = qmap_qMT_fit_spgr(data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, B0, B1, B1_MT, mode, [data_SPGR | R1], alpha_SPGR, TR_SPGR, t_m_SPGR)
%   qmap_savefig         - FUNCTION [] = qmap_savefig(fname, epsflag)
%   qmap_srsos           - FUNCTION [] = qmap_srsos(ims, nd)
%   qmap_superLorentzian - FUNCTION G  = qmap_superLorentzian(dlt, T2, init_flag, table_update, t_m)
%   qmap_t1_fit_ir       - FUNCTION [pd r1 eff res]      = qmap_t1_fit_ir(data, ti, opts)
%   qmap_t1_fit_spgr     - FUNCTION [pd, r1, pd_0, r1_0] = qmap_t1_fit_spgr(data,alpha,tr,fam,Ni, T1_0)
%   qmap_t1_fit_vafi     - FUNCTION [pd r1 fam]          = qmap_t1_fit_vafi(data_spgr, data_afi, alpha_spgr, alpha_afi, tr_spgr, tr1_afi, tr2_afi, fam1)
%   qmap_viewer          - FUNCTION [] = qmap_viewer()
%__________________________________________________________________________
% Copyright (c) 2005-2012, The Regents of the University of Wisconsin
% All rights reserved. Type 'doc qmap' for full license information.

%==========================================================================
% PROGRAMMERS NOTE:
% This <Contents.m> is the contents file for QMAP, used by
% MATLAB's ver to recover version information.
%   Line1: Toolbox Description
%   Line2: Version xx.xx dd-mmm-yyyy
%==========================================================================
