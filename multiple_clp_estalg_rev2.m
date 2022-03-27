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
clc; close all
clear all;

CD_select = 4;

%% Initialization
E_ = [55 0;
    55 0;
    27 0;
    13 0;
    13 0;
    13 0;
    13 0;
    13 0;
    ];
F_ = [55 0;
    55 0;
    27 0;
    13 0;
    13 0;
    13 0;
    13 0;
    13 0;
    ];
R_ = [11 0;
     11 0;
     5 0;
     3 0;
     3 0;
     3 0;
     3 0;
     3 0;
     ];
S_ = [11 0;
     11 0;
     5 0;
     3 0;
     3 0;
     3 0;
     3 0;
     3 0;
     ];
U_ = [4 0;
     4 0;
     1 0;
     1 0;
     1 0;
     1 0;
     1 0;
     1 0;
     ];
M_ = [48 0;
     48 0;
     256 0;
     192 0;
     192 0;
     192 0;
     192 0;
     256 0;
     ];
C_ = [3 0;   %1a
     3 0;   %1b
     48 0;  %2
     256 0; %3a
     256 0; %3b
     192 0; %4a
     192 0; %4b
     192 0; %5
     ];
 
if CD_select == 1
 
%% CD 1
% %{
numCLP = 5; Nlayer =2;

map_clp = [ 1,0;
            2,0;
            3,0;
            4,5;
            6,7;
            ];
        
TM = [48 48;
      48 48;
      128 128;
      64 64;
      96 96;
      ];
TC = [1 1;
      1 1;
      2 2;
      2 2;
      1 1;
      ];
TE = [14 14;
      14 14;
      27 13;
      13 13;
      13 13;
      ];
TF = [19 19;
      14 14;
      27 13;
      13 13;
      13 13;
      ];
%}


elseif CD_select == 2
        
%% CD 2

% %{
numCLP = 6; Nlayer =2;

map_clp = [ 1,0;
            2,0;
            3,0;
            4,0;
            5,0;
            6,7;
            ];
        
TM = [48 48;
      48 48;
      128 128;
      64 64;
      64 64;
      96 96;
      ];
TC = [1 1;
      1 1;
      1 1;
      1 1;
      1 1;
      1 1;
      ];
TE = [14 14;
      14 14;
      27 13;
      13 13;
      13 13;
      13 13;
      ];
TF = [19 19;
      14 14;
      27 13;
      13 13;
      13 13;
      13 13;
      ];
%}

elseif CD_select == 3

%% CD 3

% %{
numCLP = 5; Nlayer =2;

map_clp = [ 1,0;
            2,0;
            3,0;
            4,5;
            6,7;
            ];
        
TM = [48 48;
      48 48;
      128 128;
      64 64;
      48 48;
      ];
TC = [1 1;
      1 1;
      2 2;
      2 2;
      2 2;
      ];
TE = [14 14;
      14 14;
      27 13;
      13 13;
      13 13;
      ];
TF = [19 19;
      14 14;
      27 13;
      13 13;
      13 13;
      ];
  %}
  

else
%% CD 4
% %{
numCLP = 6; Nlayer =2;

map_clp = [ 1,0;
            2,0;
            3,0;
            4,5;
            6,0;
            7,0;
            ];
        
TM = [48 48;
      48 48;
      128 128;
      64 64;
      48 48;
      48 48;
      ];
TC = [1 1;
      1 1;
      2 2;
      2 2;
      1 1;
      1 1;
      ];
TE = [14 14;
      14 14;
      27 13;
      13 13;
      13 13;
      13 13;
      ];
TF = [19 19;
      14 14;
      27 13;
      13 13;
      13 13;
      13 13;
      ];
%}
  
end
  
    
E=zeros(numCLP,Nlayer);
F=zeros(numCLP,Nlayer);
M=zeros(numCLP,Nlayer);
C=ones(numCLP,Nlayer);
R=zeros(numCLP,Nlayer);
S=zeros(numCLP,Nlayer);
U=zeros(numCLP,Nlayer);
for i_clp =1:numCLP
    for i_l = 1:Nlayer
        if map_clp(i_clp,i_l) ~= 0
            E(i_clp,i_l) = E_(map_clp(i_clp,i_l), 1);
            F(i_clp,i_l) = F_(map_clp(i_clp,i_l), 1);
            M(i_clp,i_l) = M_(map_clp(i_clp,i_l), 1);
            C(i_clp,i_l) = C_(map_clp(i_clp,i_l), 1);
            R(i_clp,i_l) = R_(map_clp(i_clp,i_l), 1);
            S(i_clp,i_l) = S_(map_clp(i_clp,i_l), 1);
            U(i_clp,i_l) = U_(map_clp(i_clp,i_l), 1);
        end
    end
end
TH = (((TE-1).*U)+R);
TW = (((TF-1).*U)+S);




% numCLP = 6;
PixPerCycle = 2.5;

for idxCLP = 1:numCLP
    needTotalPass(idxCLP) = 0;
    for idxLayer = 1:2
        if (E(idxCLP,idxLayer) ~= 0)
            needTotalPass(idxCLP) = needTotalPass(idxCLP) + (ceil(M(idxCLP,idxLayer)/TM(idxCLP,idxLayer))*ceil(C(idxCLP,idxLayer)/TC(idxCLP,idxLayer))*ceil(E(idxCLP,idxLayer)/TE(idxCLP,idxLayer))*ceil(F(idxCLP,idxLayer)/TF(idxCLP,idxLayer)));
        end
    end
end

%needCompCycles = [0];


%% est new ver
PixPerCycle = 5;

for idx3_t = 1:numCLP
    needTotalPass(idx3_t) = 0;
    for idx4_t = 1:size(M,2)
        if (M(idx3_t, idx4_t) ~=0)
            needTotalPass(idx3_t) = needTotalPass(idx3_t) + ((M(idx3_t, idx4_t)./TM(idx3_t, idx4_t)).*(C(idx3_t, idx4_t)./TC(idx3_t, idx4_t)).*ceil(E(idx3_t, idx4_t)./TE(idx3_t, idx4_t)).*ceil(F(idx3_t, idx4_t)./TF(idx3_t, idx4_t)));
            needCompCycles(idx3_t) = 0;
        end
    end
end


    for idx3_t = 1:numCLP
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

    remnantIA = zeros(1,numCLP);
    remnantWt = zeros(1,numCLP);
    numCompCyclesRemnant = zeros(1,numCLP);
    currentPass = zeros(1,numCLP);

    remnantCompCycles = zeros(1,numCLP);
    intervalCycles = 0;
    totalCycles = 0;

    tmp_cnt = 0;

    %% Estimation Algorithm
    idx_cnt = 1;


    
est_res = zeros(1,length(100:20:400));

PixPC_array = 1.0 : 0.2 : 4.0;

for iip = 1:length(PixPC_array)
    
    PixPerCycle = PixPC_array(iip);
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
    est_res(iip) = exe_cycles_max;
    
%     for idxCLP = 1:numCLP
%         idxPass = 1;
%         for idxLayer = 1:2
%             comp_te = 0;
%             comp_tf = 0;
%             comp_r = 0;
%             comp_s = 0;
%             if (E(idxCLP,idxLayer) ~= 0)
%                 te = TE(idxCLP, idxLayer); tf = TF(idxCLP, idxLayer); tm = TM(idxCLP, idxLayer); tc = TC(idxCLP, idxLayer);
%                 e = E(idxCLP, idxLayer); f = F(idxCLP, idxLayer); m = M(idxCLP, idxLayer); c = C(idxCLP, idxLayer);
%                 r = R(idxCLP, idxLayer); s = S(idxCLP, idxLayer); u = U(idxCLP, idxLayer);
% 
%                 for idxTE = 0:te:(e-1)
%                     if (idxTE+te >= e)
%                         te = e - idxTE;
%                     end
%                     for idxTF = 0:tf:(f-1)
%                         if (idxTF+tf >= f)
%                             tf = f - idxTF;
%                         end
%                         for idxTM = 0:tm:(m-1)
%                             if (idxTM+tm >= m)
%                                 tm = m - idxTM;
%                             end
%                             for idxTC = 0:tc:(c-1)
%                                 if (idxTC+tc >= c)
%                                     tc = c - idxTC;
%                                 end
%                                 needPix_IA(idxCLP, idxPass) = tc * (((te-1)*u)+r) * (((tf-1)*u)+s);
%                                 needPix_Wt(idxCLP, idxPass) = tm * tc * r * s;
%                                 needCompCycles(idxCLP, idxPass) = comp_te * comp_tf * comp_r * comp_s;
% 
%                                 comp_te = te;
%                                 comp_tf = tf;
%                                 comp_r = r;
%                                 comp_s = s;
% 
%                                 idxPass = idxPass + 1;
%                             end
%                         end
%                         tf = TF(idxCLP);
%                     end
%                     te = TE(idxCLP);
%                 end
%                 needPix_IA(idxCLP, idxPass) = 0;
%                 needPix_Wt(idxCLP, idxPass) = 0;
%             end        
%             needCompCycles(idxCLP, idxPass) = comp_te * comp_tf * comp_r * comp_s;
%             idxPass = idxPass + 1;
%         end
%     end

%     remnantIA = zeros(1,numCLP);
%     remnantWt = zeros(1,numCLP);
%     numCompCyclesRemnant = zeros(1,numCLP);
%     currentPass = zeros(1,numCLP);
% 
%     remnantCompCycles = zeros(1,numCLP);
%     intervalCycles = 0;
%     totalCycles = 0;
% 
%     tmp_cnt = 0;
end

if CD_select==1    
    dur_est = est_res;
else
    load dur_est.mat
    dur_est = [dur_est; est_res];
end
save dur_est.mat dur_est
    

%{
%% Estimation Algorithm
idx_cnt = 1;
for PixPerCycle = 2.5:0.5:2.5
    while(1)
        
        if ((remnantIA(4) == 0) || (remnantWt(4) == 0) )
               aa =1; 
        end
        %% Calculate Number of Remnant Pixels Per Pass
        for idxCLP = 1:numCLP
            if (currentPass(idxCLP) <= needTotalPass(idxCLP))
                if ((remnantIA(idxCLP) == 0) && (remnantWt(idxCLP) == 0) && (remnantCompCycles(idxCLP) == 0))
                    currentPass(idxCLP) = currentPass(idxCLP) + 1;
                                       
                    remnantIA(idxCLP) = needPix_IA(idxCLP, currentPass(idxCLP));
                    remnantWt(idxCLP) = needPix_Wt(idxCLP, currentPass(idxCLP));
                    remnantCompCycles(idxCLP) = needCompCycles(idxCLP, currentPass(idxCLP));
                    
                    
                    currentPassTime(idxCLP, currentPass(idxCLP)) = 0;
                    
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
        
        for idx_clp =  1:numCLP
            if(currentPass(idx_clp) == needTotalPass(idx_clp)+1 && tmp_flag(idx_clp) == 0)
                exeCLPCycles(idx_clp) = totalCycles + remnantCompCycles(idx_clp);
                tmp_flag(idx_clp) = 1;
            end
            
%             if (currentPass(idx_clp) == needTotalPass(idx_clp))
%                 aaaaa = 1;
%                 totalCycles
%             end
            
            currentPassTime(idx_clp, currentPass(idx_clp)) = currentPassTime(idx_clp, currentPass(idx_clp)) + intervalCycles;
        end
        
        if all(currentPass == needTotalPass+1)            
            exeCycles(idx_cnt) = totalCycles;
            totalCycles = 0;
            idx_cnt = idx_cnt+1;
            for idxCLP = 1:numCLP
                currentPass(idxCLP) = 1;
            end
            break;
        end
        
    end
end

%}



% accum = 0;
% 
% 
% for idxCLP = 1:numCLP    
%     for idxLayer = 1:2
%         accum = accum + needTotalPass(...
%             idxCLP)*(TE(idxCLP,idxLayer)*TF(idxCLP,idxLayer)*TM(idxCLP,idxLayer)*TC(idxCLP,idxLayer)/C(idxCLP,idxLayer)+ ...
%             (TE(idxCLP,idxLayer)+R(idxCLP,idxLayer)-1)*(TF(idxCLP,idxLayer)+S(idxCLP,idxLayer)-1)*TC(idxCLP,idxLayer)+ ...
%             TM(idxCLP,idxLayer)*TC(idxCLP,idxLayer)*R(idxCLP,idxLayer)*S(idxCLP,idxLayer)  );
%     end
% end

