%{
******************************************************************************
Copyright (c) 2022 SoC Design Laboratory, Konkuk University, South Korea
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met: redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer;
redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution;
neither the name of the copyright holders nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Authors: Youngho Seo (younghoseo@konkuk.ac.kr)

Revision History
2022.03.27: Started by Sunwoo Kim
*******************************************************************************
%}

clc; clear all;
clear; close all;
%% parameters
%resource
%dsp_resource = 2240;
dsp_resource = 2880;
mac_units_resource = dsp_resource / 5;
bram_resource = 2352;
board_freqency = 100000000;
%target
target = 1.00; 
step = 0.005;
target_bw = 1.618079652446946; %690T Multi_CLP
%1.546050530665915; %485T Single_CLP 
%target_bw = 1.521962962027214; %485T Multi_CLP
%target_bw = 1.967700675392983; %690T Single_CLP 
%layer parameters
C_C = [3 48 256 192 192];
M_C = [48 128 192 192 128];
E_C = [55 27 13 13 13];
F_C = [55 27 13 13 13];
R_C = [11 5 3 3 3];
S_C = [11 5 3 3 3];
H_C = [227 31 15 15 15];
W_C = [227 31 15 15 15];
U_C = [4 1 1 1 1];

%comm. latencty params.
td_dram_act = 7;
td_dram_pre = 7;
td_dram_rd0 = 4;
td_dram_rd1 = 4;
td_dram_wr0 = 4;
td_dram_wr1 =14;

td_xbar_m2t = 6;
td_drif_t2b = 5;
td_drif_d2m = 3;

nd_dram_ob = 5; % Maximum number of DRAm bursts in open page
DBL= 8;

BL = 16;
NMO= 2;
MOL= 5;

%% opt struct
A = struct('cycles_clp', {}, 'assign_layers', {}, 'tm', {}, 'tc', {}, 'te', {}, 'tf', {}, 'bram_usages', {}, 'peak_bw', {}, 'required_bw', {});

%% assign layers
CNN_Layer = [ 1 2 3 4 5 ];
N_Layer = [ 1 2 4 6 8 ];
num_clp   = 6;
num_layer = 2*length(CNN_Layer);

Nt = length(N_Layer) ^ num_clp;
idc_1 = repmat(N_Layer, 1, Nt / length(N_Layer));
idc_2 = reshape(repmat(N_Layer, length(N_Layer), 1), 1, []);
idc_2 = repmat(idc_2, 1, Nt / (length(N_Layer)^2));
idc_3 = reshape(repmat(N_Layer, length(N_Layer)^2, 1), 1, []);
idc_3 = repmat(idc_3, 1, Nt / (length(N_Layer)^3));
idc_4 = reshape(repmat(N_Layer, length(N_Layer)^3, 1), 1, []);
idc_4 = repmat(idc_4, 1, Nt / (length(N_Layer)^4));
idc_5 = reshape(repmat(N_Layer, length(N_Layer)^4, 1), 1, []);
idc_5 = repmat(idc_5, 1, Nt / (length(N_Layer)^5));
idc_6 = reshape(repmat(N_Layer, length(N_Layer)^5, 1), 1, []);
idc_6 = repmat(idc_6, 1, Nt / (length(N_Layer)^6));

%arrange = [reshape(idc_1, [], 1) reshape(idc_2, [], 1) reshape(idc_3, [], 1) reshape(idc_4, [], 1)];
%arrange = [reshape(idc_1, [], 1) reshape(idc_2, [], 1) reshape(idc_3, [], 1) reshape(idc_4, [], 1) reshape(idc_5, [], 1)];
arrange = [reshape(idc_1, [], 1) reshape(idc_2, [], 1) reshape(idc_3, [], 1) reshape(idc_4, [], 1) reshape(idc_5, [], 1) reshape(idc_6, [], 1)];
%arrange = [reshape(idc_1, [], 1) reshape(idc_2, [], 1) reshape(idc_3, [], 1) reshape(idc_4, [], 1)];
%arrange = [reshape(idc_1, [], 1)];

arrange_layer(:, 1:num_clp) = arrange(:, 1:num_clp);
idx = find(sum(arrange_layer,2) ~= 10);
arrange_layer(idx, :) = [];
arrange_layer = sort(arrange_layer(:, 1:num_clp-1),2);
arrange_layer(:, num_clp) = num_layer - sum(arrange_layer, 2);
arrange_layer = unique(arrange_layer, 'rows');

for idx1 = 1:size(arrange_layer, 1)
    for idx2 = 1:N_Layer(length(N_Layer))
        idc_num_layer(idx1, idx2) = nnz(arrange_layer(idx1, :) == idx2);
    end
end


layer_order = perms(CNN_Layer);

idc0 = 0;
idc1 = 1;
for idx1 = 1:size(arrange_layer, 1)
    for idx2 = 1:length(layer_order)
        for idx3 = 1:num_clp
            for idx4 = 1:N_Layer(length(N_Layer))
                if (arrange_layer(idx1, idx3) > idc0)
                    idc2 = fix((idc1+1)/2);
                    layer_index((idx1-1)*length(layer_order)+idx2, idx3, idx4) = layer_order(idx2, idc2);
                    idc1 = idc1 + 1;
                else
                    layer_index((idx1-1)*length(layer_order)+idx2, idx3, idx4) = 0;
                end
                idc0 = idc0 + 1;
                if (idc1 > num_layer)
                    idc1 = 1;
                end
            end
            idc0 = 0;
        end
    end
end

%%for debug
% for tmp_idc = 1:length(layer_index)
%     layer_index_0_debug(tmp_idc, :) = layer_index(tmp_idc, 1, :);
%     layer_index_1_debug(tmp_idc, :) = layer_index(tmp_idc, 2, :);
%     layer_index_2_debug(tmp_idc, :) = layer_index(tmp_idc, 3, :);
%     layer_index_3_debug(tmp_idc, :) = layer_index(tmp_idc, 4, :);
%     layer_index_4_debug(tmp_idc, :) = layer_index(tmp_idc, 5, :);
%     layer_index_5_debug(tmp_idc, :) = layer_index(tmp_idc, 6, :);
% end
%% assign layer for S-CLP
%layer_index(1, 1, :) = [1 1 2 2 3 3 4 4 5 5];

%% min_possible_cycles
total_num_mac = 0;
for idx1 = 1:num_layer
    idx2 = fix((idx1+1)/2);
    num_mac(idx2) = C_C(idx2) * M_C(idx2) * E_C(idx2) * F_C(idx2) * R_C(idx2) * S_C(idx2);
    total_num_mac = total_num_mac + num_mac(idx2);
    min_possible_cycles = total_num_mac / mac_units_resource;
end

cand_x = 0;
cand_a = 0;
while cand_a == 0
    %% optimize_compute
    cand_x = 0;
    target
    target_performance = min_possible_cycles / target;
    for idx1 = 1:size(layer_index, 1)
        for idx2 = 1:num_clp
            num_mac_clp(idx1, idx2) = 0;
            for idx3 = 1:N_Layer(length(N_Layer))
                if (layer_index(idx1, idx2, idx3) ~= 0)
                    num_mac_clp(idx1, idx2) = num_mac_clp(idx1, idx2) + num_mac(layer_index(idx1, idx2, idx3));
                    m_c_clp_t(idx2, idx3) = C_C(layer_index(idx1, idx2, idx3)) * M_C(layer_index(idx1, idx2, idx3));
                end
            end
            m_c_clp_max(idx2) = max(m_c_clp_t(idx2,:));
            m_c_clp_min(idx2) = min(m_c_clp_t(idx2,:));
            m_c_clp(idx2) = gcd(m_c_clp_max(idx2), m_c_clp_min(idx2));
            %%m_c_clp_debug(idx1,idx2) = max(m_c_clp_t(idx2,:));
        end
        
        for idx2 = 1:num_clp
            mac_units_ratio_clp(idx2) = num_mac_clp(idx1, idx2) / total_num_mac;
            mac_units_ratio_clp(idx2) = ceil(mac_units_ratio_clp(idx2) * mac_units_resource * 0.99);          
            mac_units_ratio_clp_debug(idx1, idx2) = mac_units_ratio_clp(idx2); %% for debug
            
            rem_min = mac_units_ratio_clp(idx2);
            for idx3 = 1:m_c_clp(idx2)
                if (rem(m_c_clp(idx2), idx3) == 0)
                    if (abs(mac_units_ratio_clp(idx2) - idx3) <= rem_min)
                        rem_min = abs(mac_units_ratio_clp(idx2) - idx3);
                        mac_units_clp(idx1, idx2) = idx3;
                    end
                end
            end
        end
        if mac_units_resource - sum(mac_units_clp(idx1,1:num_clp-1)) > mac_units_clp(idx1,num_clp)
            mac_units_clp(idx1, num_clp) = mac_units_resource - sum(mac_units_clp(idx1,1:num_clp-1));
        end
        clear m_c_clp;
    end
    app_cycles_clp = num_mac_clp ./ mac_units_clp;
    sum_mac_units = sum(mac_units_clp,2);
    
    idx_cand_A_mac_units = find(sum_mac_units <= mac_units_resource);
    idx_cand_A_cycles = find(max(app_cycles_clp(:,:), [], 2) <= target_performance);
    idx_cand_A = intersect(idx_cand_A_mac_units, idx_cand_A_cycles);
    for idx1 = 1:length(idx_cand_A)
        for idx2 = 1:num_clp
            idc_clp(idx1, idx2) = 0;
            tmp_idc = 0;
            for idx_tm = 1:mac_units_clp(idx_cand_A(idx1), idx2)
                idx_tc = fix(mac_units_clp(idx_cand_A(idx1), idx2) / idx_tm);
                idx_tc_r = idx_tc;
                cycles_clp = 0;
                for idx3 = 1:N_Layer(length(N_Layer))
                    if (layer_index(idx_cand_A(idx1), idx2, idx3) ~= 0)
                        cycles(idx3) = (ceil(C_C(layer_index(idx_cand_A(idx1), idx2, idx3))/idx_tc)) * ceil((M_C(layer_index(idx_cand_A(idx1), idx2, idx3))/idx_tm)) * E_C(layer_index(idx_cand_A(idx1), idx2, idx3)) * F_C(layer_index(idx_cand_A(idx1), idx2, idx3)) * R_C(layer_index(idx_cand_A(idx1), idx2, idx3)) * S_C(layer_index(idx_cand_A(idx1), idx2, idx3));
                        cycles_clp = cycles_clp + cycles(idx3);
                    end
                end
                
                if(cycles_clp <= target_performance && idx_tm*idx_tc == mac_units_clp(idx_cand_A(idx1), idx2))
                    idc_clp(idx1, idx2) = idc_clp(idx1, idx2) + 1;                    
                    for idx3 = 1:N_Layer(length(N_Layer))
                        if (layer_index(idx_cand_A(idx1), idx2, idx3) ~= 0)
                            cand_X_cycles_layer(idx1, idx2, idc_clp(idx1, idx2),idx3) = cycles(idx3);
                        end
                    end
                    cand_X_cycles_clp(idx1, idx2, idc_clp(idx1, idx2)) = cycles_clp;
                    cand_X_tm(idx1, idx2, idc_clp(idx1, idx2)) = idx_tm;
                    cand_X_tc(idx1, idx2, idc_clp(idx1, idx2)) = idx_tc;
                    cand_X_assign_layer(idx1, idx2, idc_clp(idx1, idx2), :) = layer_index(idx_cand_A(idx1), idx2, :);
                    cand_x(idx2) = 1;
                end
            end
        end
    end
    
    if cand_x ~= 0
        [row, col] = find(~idc_clp);
        row_idx = sort(row);
        idx = 1;
        while(idx < length(row_idx))
            if row_idx(idx+1,:) == row_idx(idx,:)
                row_idx(idx+1,:) = [];
            else
                idx = idx + 1;
            end
        end
        idc_clp(row_idx, :) = [];
        cand_X_cycles_layer(row_idx, :, :, :) = [];
        cand_X_cycles_clp(row_idx, :, :) = [];
        cand_X_tm(row_idx, :, :) = [];
        cand_X_tc(row_idx, :, :) = [];
        cand_X_assign_layer(row_idx,:,:,:) = [];
    end
    
    
    min_comm_cycles = 4020800;
    if sum(cand_x) == num_clp
        %% Optimize memory
        idc = 0;
        for idxx1 = 1:length(idc_clp)
            cand_idx1 = 1:idc_clp(idxx1, 1);
            cand_idx2 = 1:idc_clp(idxx1, 2);
            cand_idx3 = 1:idc_clp(idxx1, 3);
            cand_idx4 = 1:idc_clp(idxx1, 4);
            cand_idx5 = 1:idc_clp(idxx1, 5);
            cand_idx6 = 1:idc_clp(idxx1, 6);
            Nt = idc_clp(idxx1, 1) * idc_clp(idxx1, 2) * idc_clp(idxx1, 3) * idc_clp(idxx1, 4) * idc_clp(idxx1, 5) * idc_clp(idxx1, 6);
            %Nt = idc_clp(idxx1, 1) * idc_clp(idxx1, 2) * idc_clp(idxx1, 3) * idc_clp(idxx1, 4) * idc_clp(idxx1, 5);% * idc_clp(idxx1, 6);
            
            cand_idx1_n = repmat(cand_idx1, 1, Nt/idc_clp(idxx1, 1));
            cand_idx2_n = reshape(repmat(cand_idx2, idc_clp(idxx1, 1), 1), 1, []);
            cand_idx2_n = repmat(cand_idx2_n, 1, Nt/(idc_clp(idxx1, 1)*idc_clp(idxx1, 2)));
            cand_idx3_n = reshape(repmat(cand_idx3, idc_clp(idxx1, 1)*idc_clp(idxx1, 2), 1), 1, []);
            cand_idx3_n = repmat(cand_idx3_n, 1, Nt/(idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)));
            cand_idx4_n = reshape(repmat(cand_idx4, idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3), 1), 1, []);
            cand_idx4_n = repmat(cand_idx4_n, 1, Nt/(idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)*idc_clp(idxx1, 4)));
            cand_idx5_n = reshape(repmat(cand_idx5, idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)*idc_clp(idxx1, 4), 1), 1, []);
            cand_idx5_n = repmat(cand_idx5_n, 1, Nt/(idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)*idc_clp(idxx1, 4)*idc_clp(idxx1, 5)));
            cand_idx6_n = reshape(repmat(cand_idx6, idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)*idc_clp(idxx1, 4)*idc_clp(idxx1, 5), 1), 1, []);
            cand_idx6_n = repmat(cand_idx6_n, 1, Nt/(idc_clp(idxx1, 1)*idc_clp(idxx1, 2)*idc_clp(idxx1, 3)*idc_clp(idxx1, 4)*idc_clp(idxx1, 5)*idc_clp(idxx1, 6)));
                       
            %cand_idx = [cand_idx1_n; cand_idx2_n; cand_idx3_n; cand_idx4_n; cand_idx5_n];
             
            cand_idx = [cand_idx1_n; cand_idx2_n; cand_idx3_n; cand_idx4_n; cand_idx5_n; cand_idx6_n];
             
            for idx3 = Nt/2:Nt
                for idx2 = 1:num_clp
                    cand_X_cycles_layer_tmp(idx3, idx2, :) = cand_X_cycles_layer(idxx1, idx2, cand_idx(idx2,idx3), :);
                    cand_X_cycles_clp_tmp(idx3, idx2) = cand_X_cycles_clp(idxx1, idx2, cand_idx(idx2,idx3));
                    cand_X_tc_tmp(idx3, idx2) = cand_X_tc(idxx1, idx2, cand_idx(idx2,idx3));
                    cand_X_tm_tmp(idx3, idx2) = cand_X_tm(idxx1, idx2, cand_idx(idx2,idx3));
                    cand_X_assign_layer_tmp(idx3, idx2, :) = cand_X_assign_layer(idxx1, idx2, cand_idx(idx2,idx3), :);
                end
            end
%             
%             
%             for idx2 = 1:num_clp
%                 cand_X_cycles_layer_tmp(idxx1, idx2, :) = cand_X_cycles_layer(idxx1, idx2, idc_clp(idxx1, idx2), :);
%                 cand_X_cycles_clp_tmp(idxx1, idx2) = cand_X_cycles_clp(idxx1, idx2, idc_clp(idxx1, idx2));
%                 cand_X_tc_tmp(idxx1, idx2) = cand_X_tc(idxx1, idx2, idc_clp(idxx1, idx2));
%                 cand_X_tm_tmp(idxx1, idx2) = cand_X_tm(idxx1, idx2, idc_clp(idxx1, idx2));
%                 cand_X_assign_layer_tmp(idxx1, idx2, :) = cand_X_assign_layer(idxx1, idx2, idc_clp(idxx1, idx2), :);
%             end
%         end
%         
        clear bram_resource_bw; %for graph
        tmp_flag(1: 90000) = 0; %for graph
        for idx1 = 1:length(cand_X_tc_tmp)%28%59%
            for idx2 = 1:num_clp
                idc_cand(idx2) = 0;
                for idx3 = 1:N_Layer(length(N_Layer))  
                    idc_idx3(idx2, idx3) = 0;  
                    idc_idx3_f(idx2, idx3) = 0;
                    if (cand_X_assign_layer_tmp(idx1, idx2, idx3)) ~= 0
                        if idx3 == 1
                           compare_layer_idx = 0;
                        else
                           compare_layer_idx = cand_X_assign_layer_tmp(idx1, idx2, idx3-1);
                        end
                        if cand_X_assign_layer_tmp(idx1, idx2, idx3) ~= compare_layer_idx
                            if(cand_X_assign_layer_tmp(idx1, idx2, idx3) == 1)
                                idx_te = 14;
                                idc_te = 4;
                            else
                                idx_te = E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                idc_te = 1;
                            end
%                             idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)) = 0; %for debug
                            while(idx_te > 12)
                                if(cand_X_assign_layer_tmp(idx1, idx2, idx3) == 1)
                                    idx_tf = 19;
                                    idc_tf = 3;
                                else
                                    idx_tf = F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                    idc_tf = 1;
                                end
                                while(idx_tf > 12)
                                    %% Bram Usage
                                    b_ia = ( ((idx_te - 1) * U_C(cand_X_assign_layer_tmp(idx1, idx2, idx3))) + R_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) ) * ( ((idx_tf - 1) * U_C(cand_X_assign_layer_tmp(idx1, idx2, idx3))) + S_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) );
                                    b_wt = R_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) * S_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                    b_oa = idx_te * idx_tf;
                                    
                                    if b_ia > 512
                                        num_bram_ia_t = ceil(b_ia / 512) * cand_X_tc_tmp(idx1, idx2) * 2;
                                    else
                                        num_bram_ia_t = ceil(b_ia / 256) * cand_X_tc_tmp(idx1, idx2);
                                    end
                                    
                                    if b_wt < 10
                                        num_bram_wt_t = 0;
                                    elseif b_wt > 512
                                        num_bram_wt_t = ceil(b_wt / 512) * cand_X_tc_tmp(idx1, idx2) * cand_X_tm_tmp(idx1, idx2) * 2;
                                    else
                                        num_bram_wt_t = ceil(b_wt / 256) * cand_X_tc_tmp(idx1, idx2) * cand_X_tm_tmp(idx1, idx2);
                                    end
                                    num_bram_oa_t = ceil(b_oa / 512) * cand_X_tm_tmp(idx1, idx2) * 2;
                                    bram_usage_t = num_bram_ia_t + num_bram_wt_t + num_bram_oa_t;
                                    
                                    %% Required bandwidth
                                    num_pixels_wt = (cand_X_tm_tmp(idx1, idx2) * cand_X_tc_tmp(idx1, idx2) * R_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) * S_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) * (M_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / cand_X_tm_tmp(idx1, idx2)) * (C_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / cand_X_tc_tmp(idx1, idx2)) * (E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / idx_te) * (F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / idx_tf));
                                    num_pixels_oa = (cand_X_tm_tmp(idx1, idx2) * idx_te * idx_tf * (M_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / cand_X_tm_tmp(idx1, idx2)) * (E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / idx_te) * (F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / idx_tf));
                                    num_pixels_ia = (cand_X_tc_tmp(idx1, idx2) * (U_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) * (idx_te - 1) + R_C(cand_X_assign_layer_tmp(idx1, idx2, idx3))) * (U_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) * (idx_tf - 1) + S_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)))) * ((M_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / (cand_X_tm_tmp(idx1, idx2))) * (C_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / (cand_X_tc_tmp(idx1, idx2))) * (E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / idx_te) * (F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3)) / (idx_tf) ));
                                    num_pixel_t = num_pixels_ia + num_pixels_wt + num_pixels_oa;
                        
                                    required_bw_t =  (num_pixel_t / cand_X_cycles_layer_tmp(idx1, idx2, idx3)) * 32 * board_freqency / 8000000000;
                   
                                    % for debug                    
%                                     idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)) = idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)) + 1;
%                                     required_bw_debug(idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)), cand_X_assign_layer_tmp(idx1, idx2, idx3)) = required_bw;
%                                     bram_usage_debug(idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)), cand_X_assign_layer_tmp(idx1, idx2, idx3)) = bram_usage_t;
%                                     te_debug(idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)), cand_X_assign_layer_tmp(idx1, idx2, idx3)) = idx_te;
%                                     tf_debug(idc_debug(cand_X_assign_layer_tmp(idx1, idx2, idx3)), cand_X_assign_layer_tmp(idx1, idx2, idx3)) = idx_tf;
                                     
                                    if (bram_usage_t <= bram_resource )%&& required_bw_t < target_bw)
                                        idc_idx3(idx2, idx3) = idc_idx3(idx2, idx3) + 1;
                                        idc_idx3_f(idx2, idx3) = idc_idx3(idx2, idx3);
                                        num_pixel(idx2, idx3, idc_idx3(idx2, idx3)) = num_pixel_t;
                                        num_pixel_ia(idx2, idx3, idc_idx3(idx2, idx3)) = num_pixels_ia;
                                        num_pixel_wt(idx2, idx3, idc_idx3(idx2, idx3)) = num_pixels_wt;
                                        required_bw(idx2, idx3, idc_idx3(idx2, idx3)) = required_bw_t;
                                        num_bram_wt(idx2, idx3, idc_idx3(idx2, idx3)) = num_bram_wt_t;
                                        num_bram_ia(idx2, idx3, idc_idx3(idx2, idx3)) = num_bram_ia_t;
                                        num_bram_oa(idx2, idx3, idc_idx3(idx2, idx3)) = num_bram_oa_t;
                                        cand_A_te(idx2, idx3, idc_idx3(idx2, idx3)) = idx_te;
                                        cand_A_tf(idx2, idx3, idc_idx3(idx2, idx3)) = idx_tf;
                                        
                                        MM(idx2, idx3, idc_idx3(idx2, idx3)) = M_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        CC(idx2, idx3, idc_idx3(idx2, idx3)) = C_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        EE(idx2, idx3, idc_idx3(idx2, idx3)) = E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        FF(idx2, idx3, idc_idx3(idx2, idx3)) = F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        RR(idx2, idx3, idc_idx3(idx2, idx3)) = R_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        SS(idx2, idx3, idc_idx3(idx2, idx3)) = S_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        UU(idx2, idx3, idc_idx3(idx2, idx3)) = U_C(cand_X_assign_layer_tmp(idx1, idx2, idx3));
                                        TMM(idx2, idx3, idc_idx3(idx2, idx3)) = cand_X_tm_tmp(idx1, idx2);
                                        TCC(idx2, idx3, idc_idx3(idx2, idx3)) = cand_X_tc_tmp(idx1, idx2);
                                        TEE(idx2, idx3, idc_idx3(idx2, idx3)) = idx_te;
                                        TFF(idx2, idx3, idc_idx3(idx2, idx3)) = idx_tf;
                                        
                                    end  
                                    idc_tf = idc_tf + 1;
                                    idx_tf = ceil(F_C(cand_X_assign_layer_tmp(idx1, idx2, idx3))/ idc_tf);
                                end
                                idc_te = idc_te + 1;
                                idx_te =ceil(E_C(cand_X_assign_layer_tmp(idx1, idx2, idx3))/ idc_te);
                            end
                        else
                            idc_idx3(idx2, idx3) = idc_idx3(idx2, idx3 - 1);
                            idc_idx3_f(idx2, idx3) = 1;
                            num_pixel(idx2, idx3, :) = num_pixel(idx2, idx3 - 1, :);
                            num_pixel_ia(idx2, idx3,:) = num_pixel_ia(idx2, idx3-1,:);
                            num_pixel_wt(idx2, idx3,:) = num_pixel_wt(idx2, idx3-1,:);
                            required_bw(idx2, idx3, :) = required_bw(idx2, idx3-1, :);
                            num_bram_wt(idx2, idx3, :) = num_bram_wt(idx2, idx3 - 1, :);
                            num_bram_ia(idx2, idx3, :) = num_bram_ia(idx2, idx3 - 1, :);
                            num_bram_oa(idx2, idx3, :) = num_bram_oa(idx2, idx3 - 1, :);
                            cand_A_te(idx2, idx3, :) = cand_A_te(idx2, idx3 - 1, :);
                            cand_A_tf(idx2, idx3, :) = cand_A_tf(idx2, idx3 - 1, :);
                            
                            MM(idx2, idx3, :) = MM(idx2, idx3-1, :);
                            CC(idx2, idx3, :) = CC(idx2, idx3-1, :);
                            EE(idx2, idx3, :) = EE(idx2, idx3-1, :);
                            FF(idx2, idx3, :) = FF(idx2, idx3-1, :);
                            RR(idx2, idx3, :) = RR(idx2, idx3-1, :);
                            SS(idx2, idx3, :) = SS(idx2, idx3-1, :);
                            UU(idx2, idx3, :) = UU(idx2, idx3-1, :);
                            TMM(idx2, idx3, :) = TMM(idx2, idx3-1, :);
                            TCC(idx2, idx3, :) = TCC(idx2, idx3-1, :);
                            TEE(idx2, idx3, :) = TEE(idx2, idx3-1, :);
                            TFF(idx2, idx3, :) = TFF(idx2, idx3-1, :);
                        end
                    end
                end
            end
            
            %% select opt A
            if(size(idc_idx3_f,1) == 1)
                cand_idx_t = idc_idx3_f;
            else
                cand_idx_t = [idc_idx3_f(1,:) idc_idx3_f(2,:) idc_idx3_f(3,:) idc_idx3_f(4,:) idc_idx3_f(5,:) idc_idx3_f(6,:)];
                %cand_idx_t = [idc_idx3_f(1,:) idc_idx3_f(2,:) idc_idx3_f(3,:) idc_idx3_f(4,:) idc_idx3_f(5,:)];
                %cand_idx_t = [idc_idx3_f(1,:)];
                
            end
            cand_idx = nonzeros(cand_idx_t);
            [t,nonzero_idx] = find(cand_idx_t);
            clear t;            
            q = floor((nonzero_idx-1)/N_Layer(length(N_Layer))) + 1;
            r = rem(nonzero_idx-1, N_Layer(length(N_Layer))) + 1;
            clear nonzero_idx;
            
            if nnz(cand_idx) == 10
                cand_idx1 = 1:cand_idx(1);
                cand_idx2 = 1:cand_idx(2);
                cand_idx3 = 1:cand_idx(3);
                cand_idx4 = 1:cand_idx(4);
                cand_idx5 = 1:cand_idx(5);
                cand_idx6 = 1:cand_idx(6);
                cand_idx7 = 1:cand_idx(7);
                cand_idx8 = 1:cand_idx(8);
                cand_idx9 = 1:cand_idx(9);
                cand_idx10 = 1:cand_idx(10);
                Nt = cand_idx(1) * cand_idx(2) * cand_idx(3) * cand_idx(4) * cand_idx(5) * cand_idx(6) * cand_idx(7) * cand_idx(8) * cand_idx(9) * cand_idx(10);
                
                clear cand_idx_n;
                cand_idx1_n = repmat(cand_idx1, 1, Nt/cand_idx(1));
                cand_idx_n(1,:) = cand_idx1_n;
                cand_idx2_n = reshape(repmat(cand_idx2, cand_idx(1), 1), 1, []);
                cand_idx_n(2,:) = repmat(cand_idx2_n, 1, Nt/(cand_idx(1)*cand_idx(2)));
                cand_idx3_n = reshape(repmat(cand_idx3, cand_idx(1)*cand_idx(2), 1), 1, []);
                cand_idx_n(3,:) = repmat(cand_idx3_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)));
                cand_idx4_n = reshape(repmat(cand_idx4, cand_idx(1)*cand_idx(2)*cand_idx(3), 1), 1, []);
                cand_idx_n(4,:) = repmat(cand_idx4_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)));
                cand_idx5_n = reshape(repmat(cand_idx5, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4), 1), 1, []);
                cand_idx_n(5,:) = repmat(cand_idx5_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)));
                cand_idx6_n = reshape(repmat(cand_idx6, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5), 1), 1, []);
                cand_idx_n(6,:) = repmat(cand_idx6_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)));
                cand_idx7_n = reshape(repmat(cand_idx7, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6), 1), 1, []);
                cand_idx_n(7,:) = repmat(cand_idx7_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)));
                cand_idx8_n = reshape(repmat(cand_idx8, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7), 1), 1, []);
                cand_idx_n(8,:) = repmat(cand_idx8_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)*cand_idx(8)));
                cand_idx9_n = reshape(repmat(cand_idx9, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)*cand_idx(8), 1), 1, []);
                cand_idx_n(9,:) = repmat(cand_idx9_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)*cand_idx(8)*cand_idx(9)));
                cand_idx10_n = reshape(repmat(cand_idx10, cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)*cand_idx(8)*cand_idx(9), 1), 1, []);
                cand_idx_n(10,:) = repmat(cand_idx10_n, 1, Nt/(cand_idx(1)*cand_idx(2)*cand_idx(3)*cand_idx(4)*cand_idx(5)*cand_idx(6)*cand_idx(7)*cand_idx(8)*cand_idx(9)*cand_idx(10)));
                %cand_idx_n(10,:) = reshape(repmat(cand_idx10, Nt/cand_idx(10), 1), 1, []);
                
                
                for idx2 = 1:Nt%858%1:Nt
                    num_pixel_r(1) = num_pixel(q(1),r(1),cand_idx_n(1, idx2));
                    num_pixel_ia_r(1) = num_pixel_ia(q(1),r(1),cand_idx_n(1, idx2));
                    num_pixel_wt_r(1) = num_pixel_wt(q(1),r(1),cand_idx_n(1, idx2));
                    required_bw_r(1) = required_bw(q(1),r(1),cand_idx_n(1, idx2));                    
                    num_bram_wt_r(1) = num_bram_wt(q(1),r(1),cand_idx_n(1, idx2));
                    num_bram_ia_r(1) = num_bram_ia(q(1),r(1),cand_idx_n(1, idx2));
                    num_bram_oa_r(1) = num_bram_oa(q(1),r(1),cand_idx_n(1, idx2));
                    te_r(1) = cand_A_te(q(1),r(1),cand_idx_n(1, idx2));
                    tf_r(1) = cand_A_tf(q(1),r(1),cand_idx_n(1, idx2));
                                        
                    MM_r(1) = MM(q(1),r(1),cand_idx_n(1, idx2));
                    CC_r(1) = CC(q(1),r(1),cand_idx_n(1, idx2));
                    EE_r(1) = EE(q(1),r(1),cand_idx_n(1, idx2));
                    FF_r(1) = FF(q(1),r(1),cand_idx_n(1, idx2));
                    RR_r(1) = RR(q(1),r(1),cand_idx_n(1, idx2));
                    SS_r(1) = SS(q(1),r(1),cand_idx_n(1, idx2));
                    UU_r(1) = UU(q(1),r(1),cand_idx_n(1, idx2));
                    TMM_r(1) = TMM(q(1),r(1),cand_idx_n(1, idx2));
                    TCC_r(1) = TCC(q(1),r(1),cand_idx_n(1, idx2));
                    TEE_r(1) = TEE(q(1),r(1),cand_idx_n(1, idx2));
                    TFF_r(1) = TFF(q(1),r(1),cand_idx_n(1, idx2));
                    
                    for idx3 = 2:10
                        if q(idx3) == q(idx3 - 1) && cand_X_assign_layer_tmp(idx1,q(idx3),r(idx3)) == cand_X_assign_layer_tmp(idx1,q(idx3-1),r(idx3-1))
                            num_pixel_r(idx3) = num_pixel(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            num_pixel_ia_r(idx3) = num_pixel_ia(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            num_pixel_wt_r(idx3) = num_pixel_wt(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            required_bw_r(idx3) = required_bw(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            num_bram_wt_r(idx3) = num_bram_wt(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            num_bram_ia_r(idx3) = num_bram_ia(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            num_bram_oa_r(idx3) = num_bram_oa(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            te_r(idx3) = cand_A_te(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            tf_r(idx3) = cand_A_tf(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            
                            MM_r(idx3) = MM(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            CC_r(idx3) = CC(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            EE_r(idx3) = EE(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            FF_r(idx3) = FF(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            RR_r(idx3) = RR(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            SS_r(idx3) = SS(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            UU_r(idx3) = UU(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            TMM_r(idx3) = TMM(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            TCC_r(idx3) = TCC(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            TEE_r(idx3) = TEE(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                            TFF_r(idx3) = TFF(q(idx3-1),r(idx3-1),cand_idx_n(idx3-1, idx2));
                        else
                            num_pixel_r(idx3) = num_pixel(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            num_pixel_ia_r(idx3) = num_pixel_ia(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            num_pixel_wt_r(idx3) = num_pixel_wt(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            required_bw_r(idx3) = required_bw(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            num_bram_wt_r(idx3) = num_bram_wt(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            num_bram_ia_r(idx3) = num_bram_ia(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            num_bram_oa_r(idx3) = num_bram_oa(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            te_r(idx3) = cand_A_te(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            tf_r(idx3) = cand_A_tf(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                                                        
                            MM_r(idx3) = MM(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            CC_r(idx3) = CC(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            EE_r(idx3) = EE(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            FF_r(idx3) = FF(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            RR_r(idx3) = RR(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            SS_r(idx3) = SS(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            UU_r(idx3) = UU(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            TMM_r(idx3) = TMM(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            TCC_r(idx3) = TCC(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            TEE_r(idx3) = TEE(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                            TFF_r(idx3) = TFF(q(idx3),r(idx3),cand_idx_n(idx3, idx2));
                        end
                    end
                    
                    pixel2clp = rot90([q; num_pixel_r]);
                    pixel_ia2clp = [q; num_pixel_ia_r];
                    pixel_wt2clp = [q; num_pixel_wt_r];
                    required_bw2clp = [q; required_bw_r];
                    num_bram_wt2clp = [q; num_bram_wt_r];
                    num_bram_ia2clp = [q; num_bram_ia_r];
                    num_bram_oa2clp = [q; num_bram_oa_r];
                    
                    m2clp = [q; MM_r];
                    c2clp = [q; CC_r];
                    e2clp = [q; EE_r];
                    f2clp = [q; FF_r];
                    r2clp = [q; RR_r];
                    s2clp = [q; SS_r];
                    u2clp = [q; UU_r];
                    tm2clp = [q; TMM_r];
                    tc2clp = [q; TCC_r];                    
                    te2clp = [q; te_r];
                    tf2clp = [q; tf_r];
                    
                    idc_layer = 1;
                    for idx3 = 1:num_clp
                        [pixel_t,pixel_i] = find(num_pixel_r(:,1) == idx3);
                        [ia_t,ia_i] = find(pixel_ia2clp(1,:) == idx3); 
                        [wt_t,wt_i] = find(pixel_wt2clp(1,:) == idx3); 
                        [b_t,b_i] = find(required_bw2clp(1,:) == idx3); 
                        [bw_t,bw_i] = find(num_bram_wt2clp(1,:) == idx3); 
                        [bi_t,bi_i] = find(num_bram_ia2clp(1,:) == idx3); 
                        [bo_t,bo_i] = find(num_bram_oa2clp(1,:) == idx3); 
                        
                        
                        [m_t, m_i] = find(m2clp(1,:) == idx3);
                        [c_t, c_i] = find(c2clp(1,:) == idx3);
                        [e_t, e_i] = find(e2clp(1,:) == idx3);
                        [f_t, f_i] = find(f2clp(1,:) == idx3);
                        [r_t, r_i] = find(r2clp(1,:) == idx3);
                        [s_t, s_i] = find(s2clp(1,:) == idx3);
                        [u_t, u_i] = find(u2clp(1,:) == idx3);
                        [tm_t,tm_i] = find(tm2clp(1,:) == idx3);      
                        [tc_t,tc_i] = find(tc2clp(1,:) == idx3);
                        [te_t,te_i] = find(te2clp(1,:) == idx3);      
                        [tf_t,tf_i] = find(tf2clp(1,:) == idx3);   
                        
                        for idx4 = 1:N_Layer(length(N_Layer))
                            if (cand_X_assign_layer_tmp(idx1, idx3, idx4)) ~= 0
                                required_bandwidth_layer(idx3, idx4) = (num_pixel_r(idc_layer) / cand_X_cycles_layer_tmp(idx1, idx3, idx4)) * 32 * board_freqency / 8000000000;
                                
                                idc_layer = idc_layer + 1;
                            end
                        end
                        cand_A_required_bw_ideal(idx3) = max(required_bandwidth_layer(idx3, :));
                        cand_A_num_pixel_ia(idx3) = sum(pixel_ia2clp(2, ia_i));
                        cand_A_num_pixel_wt(idx3) = sum(pixel_wt2clp(2, wt_i));
                        cand_A_required_bw(idx3) = max(required_bw2clp(2, b_i));
                        cand_A_num_bram_wt(idx3) = max(num_bram_wt2clp(2, bw_i));
                        cand_A_num_bram_ia(idx3) = max(num_bram_ia2clp(2, bi_i));
                        cand_A_num_bram_oa(idx3) = max(num_bram_oa2clp(2, bo_i));
                        te(idx3,1:length(te_i)) = te2clp(2, te_i);
                        tf(idx3,1:length(tf_i)) = tf2clp(2, tf_i);                   
                        
                        
                        M(idx3,1:length(m_i)) = m2clp(2, m_i);    
                        C(idx3,1:length(c_i)) = c2clp(2, c_i);    
                        E(idx3,1:length(e_i)) = e2clp(2, e_i);    
                        F(idx3,1:length(f_i)) = f2clp(2, f_i);    
                        R(idx3,1:length(r_i)) = r2clp(2, r_i);    
                        S(idx3,1:length(s_i)) = s2clp(2, s_i);   
                        U(idx3,1:length(s_i)) = u2clp(2, s_i);      
                        
                        TM(idx3,1:length(tm_i)) = tm2clp(2, tm_i);    
                        TC(idx3,1:length(tc_i)) = tc2clp(2, tc_i);    
                        TE(idx3,1:length(te_i)) = te2clp(2, te_i);    
                        TF(idx3,1:length(tf_i)) = tf2clp(2, tf_i);   
                                               
                        clear piex_t; clear pixel_i; clear te_t; clear te_i; clear tf_t; clear tf_i;
                        clear bw_t; clear bw_i; clear bi_t; clear bi_i; clear bo_t; clear bo_i;
                        clear m_t; clear m_i;
                        clear c_t; clear c_i;
                        clear e_t; clear e_i;
                        clear f_t; clear f_i;
                        clear r_t; clear r_i;
                        clear s_t; clear s_i;
                        clear u_t; clear u_i;
                        clear tm_t; clear tm_i;
                        clear tc_t; clear tc_i;
                    end
                    
                    read_pix_tmp = cand_A_num_pixel_ia + cand_A_num_pixel_wt;
                    [read_pix, clp_idx] = sort(read_pix_tmp); %clear read_pix_tmp;
                    num_sec = num_clp;
                    
                    
                    for sec_idx = 1:num_sec
                        if(sec_idx == 1)
                            comm_pix(sec_idx) = read_pix(sec_idx);     
                        else                           
                            comm_pix(sec_idx) = read_pix(sec_idx) - read_pix(sec_idx - 1);
                        end
                        comm_cycles_sec(sec_idx) = (num_sec - (sec_idx-1)) * comm_pix(sec_idx);                        
                    end                    
                    for sec_idx = 1:num_sec
                       find_idx = find(clp_idx == sec_idx);
                       for ratio_idx = 1:find_idx                           
                           comm_cycles(sec_idx, ratio_idx) = comm_cycles_sec(ratio_idx);
                           comp_ratio(sec_idx, ratio_idx) = comm_pix(ratio_idx) / read_pix_tmp(sec_idx);
                           comp_cycles(sec_idx, ratio_idx) = comp_ratio(sec_idx, ratio_idx) * cand_X_cycles_clp_tmp(idx1, sec_idx);
                       end                    
                    end
                    exe_cycles = max(comm_cycles, comp_ratio);
                    exe_cycles_max = max(sum(exe_cycles,2));
                                       
                    peak_bandwidth = sum(cand_A_required_bw);
                    peak_bandwidth_ideal = sum(cand_A_required_bw_ideal);
                    total_bram_usage = sum(cand_A_num_bram_wt) + sum(cand_A_num_bram_ia) + sum(cand_A_num_bram_oa);
                    
                    %% est new ver
                    numCLP = num_clp;
                    PixPerCycle = 5;
                    for idx3_t = 1:num_clp
                        needTotalPass(idx3_t) = 0;
                        for idx4_t = 1:size(M,2)
                            if (M(idx3_t, idx4_t) ~=0)
                                needTotalPass(idx3_t) = needTotalPass(idx3_t) + ((M(idx3_t, idx4_t)./TM(idx3_t, idx4_t)).*(C(idx3_t, idx4_t)./TC(idx3_t, idx4_t)).*ceil(E(idx3_t, idx4_t)./TE(idx3_t, idx4_t)).*ceil(F(idx3_t, idx4_t)./TF(idx3_t, idx4_t)));
                                needCompCycles(idx3_t) = 0;
                            end
                        end
                    end
                    
                    for idx3_t = 1:num_clp
                        idxPass = 1;
                        for idx4_t = 1:size(M,2)
                            if (M(idx3_t, idx4_t) ~=0)
                                tee = TE(idx3_t, idx4_t); tff = TF(idx3_t, idx4_t); tm = TM(idx3_t, idx4_t); tc = TC(idx3_t, idx4_t);
                                e = E(idx3_t, idx4_t); f = F(idx3_t, idx4_t); m = M(idx3_t, idx4_t); c = C(idx3_t, idx4_t);
                                rr = R(idx3_t, idx4_t); s = S(idx3_t, idx4_t); u = U(idx3_t, idx4_t);
                                
                                for idxTE = 0:tee:(e-1)
                                    if (idxTE+tee >= e)
                                        tee = e - idxTE;
                                    end
                                    for idxTF = 0:tff:(f-1)
                                        if (idxTF+tff >= f)
                                            tff = f - idxTF;
                                        end
                                        for idxTM = 0:tm:(m-1)
                                            for idxTC = 0:tc:(c-1)
                                                needPix_IA(idx3_t, idxPass) = tc * (((tee-1)*u)+rr) * (((tff-1)*u)+s);
                                                needPix_Wt(idx3_t, idxPass) = tm * tc * rr * s;
                                                needCompCycles(idx3_t, idxPass+1) = tee * tff * rr * s;
                                                idxPass = idxPass + 1;
                                            end
                                        end
                                        tff = TF(idx3_t);
                                    end
                                    tee = TE(idx3_t);
                                end
                                needPix_IA(idx3_t, idxPass) = 0;
                                needPix_Wt(idx3_t, idxPass) = 0;
                            end
                        end
                    end
                    
                    remnantIA = [0, 0, 0, 0, 0, 0];
                    remnantWt = [0, 0, 0, 0, 0, 0];
                    numCompCyclesRemnant = [0, 0, 0, 0, 0, 0];
                    currentPass = [0, 0, 0, 0, 0, 0];
                    
                    remnantCompCycles = [0, 0, 0, 0, 0, 0];
                    intervalCycles = 0;
                    totalCycles = 0;
                    
                    tmp_cnt = 0;
                    
                    %% Estimation Algorithm
                    idx_cnt = 1;
                    while(1)
                        %% Calculate Number of Remnant Pixels Per Pass
                        for idxCLP = 1:numCLP
                            if (currentPass(idxCLP) <= needTotalPass(idxCLP))
                                if ((remnantIA(idxCLP) == 0) && (remnantWt(idxCLP) == 0) && (remnantCompCycles(idxCLP) == 0))
                                    currentPass(idxCLP) = currentPass(idxCLP) + 1;
                                    
                                    remnantIA(idxCLP) = needPix_IA(idxCLP, currentPass(idxCLP));
                                    remnantWt(idxCLP) = needPix_Wt(idxCLP, currentPass(idxCLP));
                                    remnantCompCycles(idxCLP) = needCompCycles(idxCLP, currentPass(idxCLP));
                                end
                            else
                                remnantCompCycles(idxCLP) = 0;
                            end
                        end
                        %% Estimate Number of Active DMAs
                        remnantRdPix = nonzeros([remnantIA remnantWt]);
                        numActiveDMA = length(remnantRdPix);
                        
                        %% Estimate Number of Interval Cycles
                        intervalCommCycles = min(remnantRdPix) * numActiveDMA / PixPerCycle;
                        intervalCompCycles = min(nonzeros(remnantCompCycles));
                        intervalCycles = min([intervalCompCycles, intervalCommCycles]);
                        
                        if (size(intervalCommCycles) == 0)
                            intervalPix = 0;
                        else
                            intervalPix = round(intervalCycles / numActiveDMA * PixPerCycle);
                        end
                        
                        totalCycles = totalCycles + intervalCycles;
                        
                        %% Update number of remaining pixels/cycles
                        idx_remIA = find((remnantIA - intervalPix) >= 0);
                        remnantIA(idx_remIA) = remnantIA(idx_remIA) - intervalPix;
                        
                        idx_remWt = find((remnantWt - intervalPix) >= 0);
                        remnantWt(idx_remWt) = remnantWt(idx_remWt) - intervalPix;
                        
                        remnantCompCycles = max((remnantCompCycles - intervalCycles), 0);
                        
                        if all(currentPass == needTotalPass+1)
                            exe_cycles_max = totalCycles;
                            totalCycles = 0;
                            idx_cnt = idx_cnt+1;
                            for idxCLP = 1:numCLP
                                currentPass(idxCLP) = 0;
                            end
                            break;
                        end
                            
                    end
                    
                    %% 
                    tol = eps(5);
%                     if(peak_bandwidth <= 9.0)
%                     graph_idx = fix(peak_bandwidth * 10^4);                  
%                     if(tmp_flag(graph_idx) == 0)
%                         tmp_bandwidth(graph_idx) = 9.0;
%                     end
%                     end
%                     if(3000 >= total_bram_usage && peak_bandwidth <= 9.0)
%                         if(peak_bandwidth <= tmp_bandwidth(graph_idx))
%                             tmp_bandwidth(graph_idx) = peak_bandwidth;
%                             bram_resource_bw(graph_idx, :) = [peak_bandwidth total_bram_usage];
%                             bram_resource_bw_ideal(graph_idx, :) = [peak_bandwidth_ideal total_bram_usage];
%                             
%                             for idx3 = 1:num_clp                          
%                                 tm_bw(graph_idx, idx3) = cand_X_tm_tmp(idx1, idx3);
%                                 tc_bw(graph_idx, idx3) = cand_X_tc_tmp(idx1, idx3);
%                                 te_bw(graph_idx, idx3, :) = te(idx3, :);
%                                 tf_bw(graph_idx, idx3, :) = tf(idx3, :);
%                             end                            
%                             bram_bw(graph_idx, 1) = num_bram_wt2clp(2, 1) + num_bram_oa2clp(2, 1) + num_bram_ia2clp(2, 1);
%                             bram_bw(graph_idx, 2) = num_bram_wt2clp(2, 2) + num_bram_oa2clp(2, 2) + num_bram_ia2clp(2, 2);
%                             bram_bw(graph_idx, 3) = num_bram_wt2clp(2, 3) + num_bram_oa2clp(2, 3) + num_bram_ia2clp(2, 3);
%                             bram_bw(graph_idx, 4) = num_bram_wt2clp(2, 5) + num_bram_oa2clp(2, 5) + num_bram_ia2clp(2, 5);
%                             bram_bw(graph_idx, 5) = num_bram_wt2clp(2, 7) + num_bram_oa2clp(2, 7) + num_bram_ia2clp(2, 7);
%                             bram_bw(graph_idx, 6) = num_bram_wt2clp(2, 9) + num_bram_oa2clp(2, 9) + num_bram_ia2clp(2, 9);
%                             tmp_flag(graph_idx) = 1;
%                         end
%                     end
%                     if total_bram_usage == 1238
%                        xxxx= 1; 
%                     end
                    
%                     min_idx = 1;
%                     if (bram_resource >= total_bram_usage)
%                         if(idc == 0)
%                             min_comm_cycles = exe_cycles_max+1;
%                             min_idx = Nt +1;
%                         end
%                     else
%                         min_idx = 1 + min_idx;
%                     end

                    if (bram_resource >= total_bram_usage && min_comm_cycles >= exe_cycles_max)
                        idc = idc + 1;
                        min_comm_cycles = exe_cycles_max;
                        A(idc).comm_cycle= comm_cycles;
                        A(idc).comp_cycles = comp_cycles;
                        A(idc).exe_time = exe_cycles_max;
                        plot_exe_time(idc) = exe_cycles_max;
                        A(idc).cycles_clp = cand_X_cycles_clp_tmp(idx1, :);
                        A(idc).peak_bw = peak_bandwidth;
                        A(idc).bram_usages = total_bram_usage;
                        A(idc).assign_layers(:, :) = cand_X_assign_layer_tmp(idx1, :, :);
                        A(idc).bram_wt = cand_A_num_bram_wt;
                        A(idc).bram_ia = cand_A_num_bram_ia;
                        A(idc).bram_oa = cand_A_num_bram_oa;
                        %% CLPs
                        for idx3 = 1:num_clp
                            A(idc).tm(idx3) = cand_X_tm_tmp(idx1, idx3);
                            A(idc).tc(idx3) = cand_X_tc_tmp(idx1, idx3);
                            A(idc).te(idx3, :) = te(idx3, :);
                            A(idc).tf(idx3, :) = tf(idx3, :);
                            A(idc).required_bw(idx3) = cand_A_required_bw(idx3);
                        end                        
                        cand_a = 1;
                        
                    end
                    
                    clear comm_pix;
                    clear comm_cycles;
                    clear comp_ratio;
                    clear comp_cycles
                end
           end
            
        end  
        
            clear cand_idx1_n;
            clear cand_idx2_n;
            clear cand_idx3_n;
            clear cand_idx4_n;
            clear cand_idx5_n;
            clear cand_idx6_n;
        end
        
        if cand_a == 0
            target = target - step;
        end
    else
        target = target - step;
    end
end

% 
% for idx1111 = 1:size(bram_resource_bw_ideal, 1)
%         find_tmp = find(sum(bram_resource_bw_ideal, 2) == 0);   
% end                  
% bram_resource_bw(find_tmp,:) = [];
% bram_resource_bw_ideal(find_tmp,:) = [];
% tm_bw(find_tmp,:) = [];
% tc_bw(find_tmp,:) = [];
% te_bw(find_tmp,:,:) = [];
% tf_bw(find_tmp,:,:) = [];
% bram_bw(find_tmp,:) = [];
% 
% 
% idx = 1;
% while(idx < size(bram_resource_bw_ideal, 1))
%     if bram_resource_bw_ideal(idx+1,2) > bram_resource_bw_ideal(idx,2)
%         bram_resource_bw_ideal(idx+1,:) = [];
%         bram_resource_bw(idx+1,:) = [];
%         tm_bw(idx+1,:) = [];
%         tc_bw(idx+1,:) = [];
%         bram_bw(idx+1,:) = [];
%         te_bw(idx+1,:, :) = [];
%         tf_bw(idx+1,:, :) = [];
%     else
%         idx = idx + 1;
%     end
% end
% 
% for idx1 = 1:size(te_bw,1)
%    e_f_te_tf(idx1, 1) = ceil(E_C(1) / te_bw(idx1, 1)) * ceil(F_C(1) / tf_bw(idx1, 1));
%    e_f_te_tf(idx1, 2) = ceil(E_C(1) / te_bw(idx1, 2)) * ceil(F_C(1) / tf_bw(idx1, 2));
%    e_f_te_tf(idx1, 3) = ceil(E_C(2) / te_bw(idx1, 3)) * ceil(F_C(2) / tf_bw(idx1, 3));
%    e_f_te_tf(idx1, 4) = ceil(E_C(3) / te_bw(idx1, 4)) * ceil(F_C(3) / tf_bw(idx1, 4));
%    e_f_te_tf(idx1, 5) = ceil(E_C(5) / te_bw(idx1, 5)) * ceil(F_C(5) / tf_bw(idx1, 5));
%    e_f_te_tf(idx1, 6) = ceil(E_C(4) / te_bw(idx1, 6)) * ceil(F_C(4) / tf_bw(idx1, 6));
% end
% 
% 
% b = load('C:\Users\user\Desktop\bram_bw_485.txt');
% 
% y = bram_resource_bw(:,1);
% x = bram_resource_bw(:,2);
% 
% 
% yy = b(:,1);
% xx = b(:,2);
% 
% title('Tradeoff between BRAM and bandwidth');
% xlabel('BRAMs');
% ylabel('Bandwidth');
% axis([0, 2000, 0, 4]);
% grid on
% hold on
% plot(x,y, 'mo', 'LineWidth', 2);
% plot(xx,yy, 'go', 'LineWidth', 2);
% 
% legend('695T (DSP 2880)', '485T (DSP 2240)');





















