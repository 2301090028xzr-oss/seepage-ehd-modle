function porous_ehl_final_v1()
%==========================================================================
% porous_ehl_final_v1.m
%
% 功能：
%   多孔点接触弹流润滑（EHL）主程序（A1+A2 + 工况预设整合版，终版）
%
% 当前版本实现：
%   1) 二维界面 Reynolds 方程
%   2) 三维多孔层稳态渗流方程：∇·(kappa ∇Pp)=0
%   3) 标量骨架法向位移 Uz 方程（顶面牵引边界离散）
%   4) 膜厚-压力-渗流-骨架变形外层耦合
%   5) 载荷平衡修正（H00 采用带阻尼 Secant 法）
%   6) Reynolds 几何多重网格 V-cycle 加速
%   7) 输出中心线 Y=0 处的二维压力分布和膜厚分布图
%   8) 支持“无孔基线工况”开关
%   9) A1：Uz 从主弹性变形降级为“多孔附加修正项”
%  10) A2：新增非局部表面弹性卷积模块，作为主 EHL 变形
%  11) 新增工况预设：
%         - baseline_ehl_shape
%         - current_reference
%         - porous_diag_strong
%  12) 新增 Y 方向镜像对称化：
%         - Reynolds 求解后对 P 做镜像平均
%         - 渗流求解后对 Pp 做镜像平均
%         - 骨架求解后对 Uz 做镜像平均
%  13) 新增冻结形态后的最终载荷闭合：
%         - 外层结束后，冻结 deltaEHL/deltaPorous/Pp/Uz 形态
%         - 仅调 H00 并重解 Reynolds，将载荷误差进一步压紧
%
% 统一无量纲约定：
%   X = x / b
%   Y = y / b
%   Z = z / b
%   P = p / pH
%   P* = pp / pH
%   H = hR / b^2
%   Uz = uzR / b^2
%==========================================================================

    clc; clear; close all;

    fprintf('============================================================\n');
    fprintf(' 多孔点接触EHL耦合求解程序（A1+A2 + 工况预设整合版，终版）\n');
    fprintf(' 耦合链：P -> P* -> U -> H -> P\n');
    fprintf('============================================================\n');

    %% 1. 初始化参数
    par = initParameters_local();
    par = applyDefaultSkeletonParams_local(par);

    if par.Nx < 3 || par.Ny < 3 || par.Nz < 3
        error('Nx、Ny、Nz 均必须 >= 3；当前代码的二阶边界/梯度离散至少需要 3 个节点。');
    end

    fprintf('当前工况预设：%s\n', par.casePreset);
    if par.runBaselineNoPorous
        fprintf('当前工况：无孔基线工况（Pi_p=0，关闭渗流与骨架变形）\n');
    else
        fprintf('当前工况：多孔耦合工况\n');
    end

    fprintf('\n闭合得到的无量纲参数：\n');
    fprintf('Lambda      = %.8e\n', par.Lambda);
    fprintf('Pi_p        = %.8e\n', par.Pi_p);
    fprintf('Gamma       = %.8e\n', par.Gamma);
    fprintf('Lambda_s    = %.8e\n', par.Lambda_s);
    fprintf('G_dim       = %.8e Pa\n', par.G_dim);
    fprintf('lambda_dim  = %.8e Pa\n', par.lambda_dim);
    fprintf('------------------------------------------------------------\n');

    %% 2. 建立二维/三维网格
    grid2D = buildGrid2D_local(par);
    grid3D = buildGrid3D_local(par);

    %% 2.1 预计算点接触 EHL 表面弹性影响核（A2）
    par = prepareElasticInfluence_local(par, grid2D);

    %% 3. 初始化二维/三维场变量
    field2D = initField2D_local(par, grid2D);
    field3D = initField3D_local(par);

    %% 4. Hertz 初始压力
    field2D.P = initPressureHertz_local(grid2D);

    %% 5. 初始化三维孔压场
    field3D.Pp(:,:,1) = field2D.P;
    for k = 2:par.Nz
        field3D.Pp(:,:,k) = par.P0_bottom;
    end

    field3D.Ux(:) = 0.0;
    field3D.Uy(:) = 0.0;
    field3D.Uz(:) = 0.0;

    %% 6. 初始更新
    [field2D, field3D] = updateFilmAndProps_local(par, grid2D, grid3D, field2D, field3D);

    %% 7. 外层耦合迭代
    fprintf('\n开始外层耦合迭代...\n');

    for itOuter = 1:par.maxOuter
        P_old        = field2D.P;
        H_old        = field2D.H;
        Pp_old       = field3D.Pp;
        Uz0_old      = field3D.Uz(:,:,1);
        H00_old      = par.H00;
        deltaEHL_old = field2D.deltaEHL;

        %--------------------------------------------------------------
        % Step A: 渗流
        %--------------------------------------------------------------
        if par.enableSeepage
            field3D = solveSeepage3D_local(par, grid3D, field2D, field3D);
        else
            field3D.Pp(:,:,1)     = field2D.P;
            field3D.Pp(:,:,2:end) = par.P0_bottom;
        end

        %--------------------------------------------------------------
        % Step B: 骨架
        %--------------------------------------------------------------
        if par.enableSkeleton
            field3D = solveSkeleton3D_local(par, grid3D, field2D, field3D);
        else
            field3D.Ux(:) = 0.0;
            field3D.Uy(:) = 0.0;
            field3D.Uz(:) = 0.0;
            field3D.Ev(:) = 0.0;
        end

        %--------------------------------------------------------------
        % Step C: 更新膜厚、物性、系数
        %--------------------------------------------------------------
        [field2D, field3D] = updateFilmAndProps_local(par, grid2D, grid3D, field2D, field3D);

        %--------------------------------------------------------------
        % Step D: Reynolds
        %--------------------------------------------------------------
        field2D = solveReynolds2D_local(par, grid2D, field2D, field3D);

        %--------------------------------------------------------------
        % Step E: 载荷平衡修正（Secant 法）
        %--------------------------------------------------------------
        loadNow = computeLoad2D_local(grid2D, field2D.P);
        par     = adjustH00_local(par, loadNow);

        % H00 改变后，先更新膜厚和系数
        [field2D, field3D] = updateFilmAndProps_local(par, grid2D, grid3D, field2D, field3D);

        % 再解一次 Reynolds，使 P 与 H/H00 保持一致
        field2D = solveReynolds2D_local(par, grid2D, field2D, field3D);

        % 再同步一次物性、Xi 和界面梯度
        [field2D, field3D] = updateFilmAndProps_local(par, grid2D, grid3D, field2D, field3D);

        % 最终载荷（本轮用于收敛判断）
        loadNow = computeLoad2D_local(grid2D, field2D.P);

        %--------------------------------------------------------------
        % 收敛判据
        %--------------------------------------------------------------
        errP      = max(abs(field2D.P(:) - P_old(:)));
        errH      = max(abs(field2D.H(:) - H_old(:)));
        errPp     = max(abs(field3D.Pp(:) - Pp_old(:)));
        errUz0    = max(abs(field3D.Uz(:,:,1) - Uz0_old), [], 'all');
        errH00    = abs(par.H00 - H00_old);
        errDelta  = max(abs(field2D.deltaEHL(:) - deltaEHL_old(:)));
        errW      = abs(loadNow - par.Wtarget);

        fprintf(['Outer = %3d | errP = %.3e | errH = %.3e | errPp = %.3e | ' ...
                 'errUz0 = %.3e | errH00 = %.3e | errDelta = %.3e | errW = %.3e | ' ...
                 'cavFrac = %.3e | Load = %.6f | H00 = %.6f\n'], ...
                 itOuter, errP, errH, errPp, errUz0, errH00, errDelta, errW, ...
                 field2D.cavFracRe, loadNow, par.H00);

        if max([errP, errH, errPp, errUz0, errH00, errDelta]) < par.tolOuter && errW < par.tolLoad
            fprintf('>>> 外层耦合迭代已收敛。\n');
            break;
        end
    end

    if itOuter == par.maxOuter
        fprintf('>>> 已达到最大外层迭代次数，程序停止。\n');
    end

    %% 7.1 冻结形态后的最终载荷闭合（仅调 H00，不再更新 deltaEHL/deltaPorous）
    freeze.deltaEHL = field2D.deltaEHL;
    freeze.deltaPorous = field2D.deltaPorous;
    freeze.Pp = field3D.Pp;
    freeze.Uz = field3D.Uz;
    freeze.Ev = field3D.Ev;
    freeze.K  = field3D.K;

    % 重新置空 Secant 历史，避免把外层漂移映射带入冻结闭合阶段
    par.H00_prev_secant = NaN;
    par.F_prev_secant   = NaN;

    for itPost = 1:par.maxPostH00
        finalLoad = computeLoad2D_local(grid2D, field2D.P);
        errW_post = abs(finalLoad - par.Wtarget);

        if errW_post < par.tolLoad
            fprintf('>>> 冻结形态后的最终载荷闭合已满足阈值。\n');
            break;
        end

        par = adjustH00_local(par, finalLoad);

        [field2D, field3D] = updateFilmAndPropsFrozenShape_local( ...
            par, grid2D, grid3D, field2D, field3D, freeze);
        field2D = solveReynolds2D_local(par, grid2D, field2D, field3D);
        [field2D, field3D] = updateFilmAndPropsFrozenShape_local( ...
            par, grid2D, grid3D, field2D, field3D, freeze);

        finalLoad = computeLoad2D_local(grid2D, field2D.P);

        fprintf('Post-FreezeH00 = %2d | Load = %.8f | errW = %.3e | H00 = %.8f\n', ...
                itPost, finalLoad, abs(finalLoad - par.Wtarget), par.H00);
    end

    %% 8. 输出关键数值
    finalLoad = computeLoad2D_local(grid2D, field2D.P);

    fprintf('\n================ 计算结束 ================\n');
    fprintf('最终载荷积分       = %.8f\n', finalLoad);
    fprintf('目标载荷           = %.8f\n', par.Wtarget);
    fprintf('最终载荷误差       = %.8e\n', abs(finalLoad - par.Wtarget));
    fprintf('最终 H00           = %.8f\n', par.H00);
    fprintf('最大压力 Pmax      = %.8f\n', max(field2D.P(:)));
    fprintf('最小膜厚 Hmin      = %.8f\n', min(field2D.H(:)));
    fprintf('最大表面位移 Uzmax = %.8f\n', max(field3D.Uz(:,:,1), [], 'all'));
    fprintf('最小表面位移 Uzmin = %.8f\n', min(field3D.Uz(:,:,1), [], 'all'));
    fprintf('Reynolds 空化截断累计次数 cavCountRe = %d\n', field2D.cavCountRe);
    fprintf('Reynolds 空化截断比例 cavFracRe     = %.8e\n', field2D.cavFracRe);
    fprintf('==========================================\n');

    %% 8.1 位移方向快速自检
    Uzsurf    = field3D.Uz(:,:,1);
    pressMask = field2D.P > 0.10 * max(field2D.P(:));

    if any(pressMask, 'all')
        meanUz_all      = mean(Uzsurf(:), 'omitnan');
        meanUz_loadZone = mean(Uzsurf(pressMask), 'omitnan');
        meanH_loadZone  = mean(field2D.H(pressMask), 'omitnan');

        fprintf('\n位移方向快速自检：\n');
        fprintf('全场平均表面位移 Uz_mean(all)       = %.8e\n', meanUz_all);
        fprintf('受压区平均表面位移 Uz_mean(load)    = %.8e\n', meanUz_loadZone);
        fprintf('受压区平均膜厚     H_mean(load)     = %.8e\n', meanH_loadZone);

        if meanUz_loadZone < 0
            fprintf('诊断：当前符号体系下，受压区 Uz 为负；本版将 -Uz 作为多孔附加修正项的一部分。\n');
        elseif meanUz_loadZone > 0
            fprintf('警告：当前受压区 Uz 为正，需核对顶面牵引符号、孔压梯度符号与膜厚方程符号是否一致。\n');
        else
            fprintf('提示：受压区平均 Uz 接近 0，当前工况下骨架响应较弱或参数过小。\n');
        end
    else
        fprintf('\n位移方向快速自检：未找到明显受压区，跳过。\n');
    end

    %% 8.2 A1+A2 快速诊断
    [~, i0] = min(abs(grid2D.X));
    [~, j0] = min(abs(grid2D.Y));

    [PmaxVal, idxPmax] = max(field2D.P(:));
    [iPmax, jPmax] = ind2sub(size(field2D.P), idxPmax);

    [~, idxHmin] = min(field2D.H(:));
    [iHmin, jHmin] = ind2sub(size(field2D.H), idxHmin);

    centerPressureRatio = field2D.P(i0,j0) / max(PmaxVal, 1e-14);

    fprintf('\nA1+A2 快速诊断：\n');
    fprintf('casePreset                    = %s\n', par.casePreset);
    fprintf('中心点主 EHL 变形 deltaEHL    = %.8e\n', field2D.deltaEHL(i0,j0));
    fprintf('中心点多孔附加项 deltaPorous  = %.8e\n', field2D.deltaPorous(i0,j0));
    fprintf('中心压力比 P(0,0)/Pmax        = %.8e\n', centerPressureRatio);
    fprintf('Pmax 位置                     = (%.6f, %.6f)\n', grid2D.X(iPmax), grid2D.Y(jPmax));
    fprintf('Hmin 位置                     = (%.6f, %.6f)\n', grid2D.X(iHmin), grid2D.Y(jHmin));

    if field2D.deltaEHL(i0,j0) > abs(field2D.deltaPorous(i0,j0))
        fprintf('诊断：主 EHL 非局部弹性变形已超过多孔附加项，A1+A2 的主导权分配是合理的。\n');
    else
        fprintf('警告：多孔附加项已接近或超过主 EHL 变形，建议进一步减小 porousDispWeight。\n');
    end

    if centerPressureRatio > 0.30
        fprintf('诊断：中心区压力已有明显抬升，压力形态开始偏离“单纯上游单峰楔形解”。\n');
    else
        fprintf('提示：中心区压力仍偏低，说明主 EHL 变形还不够强，或出口区仍受简化边界控制。\n');
    end

    %% 8.3 Y方向对称性诊断
    symErrP  = max(abs(field2D.P - field2D.P(:, end:-1:1)), [], 'all');
    symErrH  = max(abs(field2D.H - field2D.H(:, end:-1:1)), [], 'all');
    symErrUz = max(abs(field3D.Uz(:,:,1) - field3D.Uz(:, end:-1:1, 1)), [], 'all');

    fprintf('\nY方向对称性诊断：\n');
    fprintf('symErrP  = %.8e\n', symErrP);
    fprintf('symErrH  = %.8e\n', symErrH);
    fprintf('symErrUz = %.8e\n', symErrUz);

    if max([symErrP, symErrH, symErrUz]) < 1e-8
        fprintf('诊断：Y方向对称性很好，偏侧 Hmin 更可能只是对称双谷中的一个。\n');
    elseif max([symErrP, symErrH, symErrUz]) < 1e-5
        fprintf('提示：Y方向存在轻微数值不对称，但通常仍在可接受范围内。\n');
    else
        fprintf('警告：Y方向不对称较明显，建议检查离散、边界条件与迭代过程中的数值偏置。\n');
    end

    %% 9. 绘图输出
    plotResults_local(grid2D, field2D, field3D);

end

%% =========================================================================
function par = initParameters_local()
% 初始化参数（含 W1/W2/W3 闭合 + A1/A2 + 工况预设 + 稳定化参数）

%======================== 1. 网格与计算域 ========================
par.Nx   = 161;
par.Ny   = 121;
par.Nz   = 15;

par.Xmin = -3.0;
par.Xmax =  2.0;

par.Ymin = -2.0;
par.Ymax =  2.0;

par.Hp   = 1.0;

%======================== 2. 有量纲参考物性与尺度 ========================
par.pH   = 2.83e6;    % Pa
par.eta0 = 0.08;      % Pa.s
par.rho0 = 850;       % kg/m^3

par.Ue_dim = 0.50;    % m/s
par.R_dim  = 1.0e-2;  % m
par.b_dim  = 1.30e-4; % m
par.k0_dim = 7.4e-13; % m^2

par.E_dim  = 3.0e10;  % Pa
par.nu_dim = 0.34;    % -

%======================== 3. 多孔骨架/渗流其他参数 ========================
par.alpha = 0.8;
par.Mhat  = 1.0;

%======================== 4. 由有量纲参数闭合无量纲群 ========================
par = closeDimensionlessGroups_local(par);

%======================== 5. 膜厚参数 ========================
par.H00    = 0.05;
par.Hmin   = 1.0e-8;
par.H00min = -0.20;
par.H00max = 2.0;

%======================== 5.1 主弹性变形（A1/A2）参数 ========================
par.enableEHLDeflection = true;
par.ehlDeflectionWeight = 1.00;
par.porousDispWeight    = 0.20;

par.elasticCoeffMode      = 'physical';   % 'physical' or 'legacy'
par.elasticInfluenceCoeff = NaN;
par.elasticCoreEpsFactor  = 0.5;

par.deflRelax   = 0.15;
par.deltaEHLMax = 2.0;

%======================== 6. 载荷平衡参数（Secant） ========================
par.Wtarget = 2*pi/3;

par.loadGain = 0.10;

par.useSecantH00  = true;
par.secantDamping = 0.80;
par.secantTolDen  = 1.0e-12;
par.maxH00Step    = 0.25;
par.maxPostH00    = 12;

par.H00_prev_secant = NaN;
par.F_prev_secant   = NaN;

%======================== 7. 多孔层底部边界 ========================
par.P0_bottom = 0.0;

%======================== 8. 黏度模型参数 ========================
par.Zeta      = 0.68;
par.Aeta      = log(par.eta0) + 9.67;
par.betaEta   = 30 * (5.1e-9 * par.pH);
par.etaExpMin = -20;
par.etaExpMax =  10;

%======================== 9. 密度模型参数 ========================
par.A3 = 0.059 / (1.0e-9 * par.pH);

%======================== 10. 数值保护参数 ========================
par.Kmin       = 1.0e-6;
par.Kmax       = 1.0e6;

par.kappaFloor = 1.0e-12;
par.kappaMax   = 1.0e3;

par.gradPpMax  = 1.0e4;
par.sourceMax  = 1.0e6;

par.XiMin      = 1.0e-12;
par.XiMax      = 1.0e6;

%======================== 11. 最大迭代次数 ========================
par.maxOuter = 80;
par.maxRe    = 60;
par.maxSeep  = 150;
par.maxSkel  = 120;

%======================== 12. 收敛阈值 ========================
par.tolOuter = 1.0e-5;
par.tolRe    = 1.0e-5;
par.tolSeep  = 1.0e-6;
par.tolSkel  = 1.0e-6;
par.tolLoad  = 1.0e-5;

%======================== 13. 松弛系数 ========================
par.omegaP    = 0.70;
par.omegaSeep = 0.85;

%======================== 14. 界面梯度离散 ========================
par.useSecondOrderDz = true;

%======================== 15. Reynolds 多重网格参数 ========================
par.useMG   = true;
par.mgLevel = 5;
par.mgPre   = 3;
par.mgPost  = 3;

%======================== 16. 骨架边界与骨架松弛参数 ========================
par.skelTopTractionScale = 1.0;
par.skelRelax            = 0.75;

%======================== 17. 渗流内层重算 kappa ========================
par.seepUpdateKappa      = true;
par.seepKappaUpdateEvery = 5;

%======================== 18. 工况开关参数 ========================
par.runBaselineNoPorous = false;
par.enableSeepage       = true;
par.enableSkeleton      = true;

%======================== 18.1 工况预设选择器 ========================
par.casePreset = 'spike_probe_safe';
par = applyCasePreset_local(par);

%======================== 19. 基线模式覆盖逻辑 ========================
if par.runBaselineNoPorous
    par.Pi_p           = 0.0;
    par.enableSeepage  = false;
    par.enableSkeleton = false;
end

end

%% =========================================================================
function par = closeDimensionlessGroups_local(par)
% W1/W2/W3 闭合

    if par.pH <= 0 || par.eta0 <= 0 || par.rho0 <= 0
        error('pH、eta0、rho0 必须为正。');
    end
    if par.Ue_dim <= 0 || par.R_dim <= 0 || par.b_dim <= 0
        error('Ue_dim、R_dim、b_dim 必须为正。');
    end
    if par.k0_dim < 0
        error('k0_dim 不能为负。');
    end
    if par.E_dim <= 0
        error('E_dim 必须为正。');
    end
    if abs(1 - 2*par.nu_dim) < 1e-10 || abs(1 + par.nu_dim) < 1e-10
        error('nu_dim 取值导致 Lamé 常数发散，请检查泊松比。');
    end

    par.G_dim      = par.E_dim / (2 * (1 + par.nu_dim));
    par.lambda_dim = par.E_dim * par.nu_dim / ((1 + par.nu_dim) * (1 - 2*par.nu_dim));

    par.Lambda = 12 * par.eta0 * par.Ue_dim * par.R_dim^2 / ...
                 (par.pH * par.b_dim^3);

    par.Pi_p = par.k0_dim * par.pH * par.R_dim / ...
               (par.eta0 * par.Ue_dim * par.b_dim^2);

    par.Gamma    = par.G_dim      * par.b_dim / (par.R_dim * par.pH);
    par.Lambda_s = par.lambda_dim * par.b_dim / (par.R_dim * par.pH);
end

%% =========================================================================
function par = applyCasePreset_local(par)
% 工况预设：
% 1) baseline_ehl_shape  : 基线 EHL 形态检查
% 2) current_reference   : 当前参考工况
% 3) porous_diag_strong  : 强多孔诊断工况

    tag = lower(strtrim(par.casePreset));

    switch tag
        case 'baseline_ehl_shape'

            par.runBaselineNoPorous = true;
            par.enableSeepage       = false;
            par.enableSkeleton      = false;

            par.enableEHLDeflection = true;
            par.ehlDeflectionWeight = 1.00;
            par.porousDispWeight    = 0.00;

            par.Ue_dim = 0.80;
            par.k0_dim = 1.0e-16;
            par.E_dim  = 3.0e10;
            par.alpha  = 0.80;
            par.Mhat   = 1.0;

            par.skelTopTractionScale = 1.0;

            par.H00 = -0.05;

            par.maxOuter = 100;
            par.maxRe    = 80;
            par.maxSeep  = 120;
            par.maxSkel  = 120;

            par.secantDamping = 0.50;
            par.maxH00Step    = 0.10;

            par.H00min = -0.20;
            par.H00max = 2.0;

            par.deflRelax   = 0.15;
            par.deltaEHLMax = 2.0;
            par.maxPostH00  = 12;

        case 'current_reference'

            par.runBaselineNoPorous = false;
            par.enableSeepage       = true;
            par.enableSkeleton      = true;

            par.enableEHLDeflection = true;
            par.ehlDeflectionWeight = 1.00;
            par.porousDispWeight    = 0.20;

            par.Ue_dim = 0.50;
            par.k0_dim = 1.0e-16;
            par.E_dim  = 3.0e10;
            par.alpha  = 0.80;
            par.Mhat   = 1.0;

            par.skelTopTractionScale = 1.0;

            par.H00 = 0.10;

            par.maxOuter = 100;
            par.maxRe    = 80;
            par.maxSeep  = 150;
            par.maxSkel  = 120;

            par.secantDamping = 0.75;
            par.maxH00Step    = 0.25;

            par.H00min = -0.20;
            par.H00max = 2.0;

            par.deflRelax   = 0.15;
            par.deltaEHLMax = 2.0;
            par.maxPostH00  = 12;

        case 'porous_diag_strong'

            par.runBaselineNoPorous = false;
            par.enableSeepage       = true;
            par.enableSkeleton      = true;

            par.enableEHLDeflection = true;
            par.ehlDeflectionWeight = 1.00;
            par.porousDispWeight    = 0.30;

            par.Ue_dim = 0.20;
            par.k0_dim = 1.0e-14;
            par.E_dim  = 3.0e10;
            par.alpha  = 0.95;
            par.Mhat   = 2.0;

            par.skelTopTractionScale = 5.0;

            par.H00 = 0.05;

            par.maxOuter = 140;
            par.maxRe    = 90;
            par.maxSeep  = 220;
            par.maxSkel  = 180;

            par.secantDamping = 0.60;
            par.maxH00Step    = 0.15;

            par.H00min = -1.0;
            par.H00max = 2.0;

            par.deflRelax   = 0.35;
            par.deltaEHLMax = 4.0;
            par.maxPostH00  = 12;
        
case 'spike_probe_safe'
    %==========================================================
    % 二次压力峰探测（稳妥版）
    % 目标：
    %   在仍可载荷闭合的前提下，增强经典EHL二次峰形成倾向
    %==========================================================

    % A. 近似无孔基线
    par.runBaselineNoPorous = true;
    par.enableSeepage       = false;
    par.enableSkeleton      = false;

    % B. 主EHL弹性开启
    par.enableEHLDeflection = true;
    par.ehlDeflectionWeight = 1.00;
    par.porousDispWeight    = 0.00;

    % C. 网格加密
    par.Nx = 201;
    par.Ny = 121;
    par.Nz = 15;

    par.Xmin = -3.0;
    par.Xmax =  2.0;
    par.Ymin = -2.0;
    par.Ymax =  2.0;

    % D. 提高刚度，但不走极端
    par.E_dim  = 1.0e10;
    par.nu_dim = 0.34;

    % E. 卷吸与基准黏度适度提高
    par.Ue_dim = 0.80;
    par.eta0   = 0.10;

    % F. 中等增强黏压效应
    par.Zeta    = 0.68;
    par.Aeta    = log(par.eta0) + 9.67;
    par.betaEta = 8 * (5.1e-9 * par.pH);

    % G. 中等增强密压效应
    par.A3 = (0.59 / (1.0e-9 * par.pH)) / 5;

    % H. 多孔参数保留但失活
    par.k0_dim = 1.0e-16;
    par.alpha  = 0.80;
    par.Mhat   = 1.0;
    par.skelTopTractionScale = 1.0;

    % I. H00 与边界
    par.H00 = 0.00;

    par.maxOuter = 160;
    par.maxRe    = 120;
    par.maxSeep  = 80;
    par.maxSkel  = 80;

    par.secantDamping = 0.50;
    par.maxH00Step    = 0.08;

    par.H00min = -1.0;
    par.H00max =  4.0;

    % J. 非局部弹性变形参数
    par.deflRelax   = 0.20;
    par.deltaEHLMax = 1.5;

    % K. 收敛阈值
    par.tolOuter = 1.0e-5;
    par.tolRe    = 1.0e-5;
    par.tolLoad  = 1.0e-5;
        otherwise
            error('未知工况预设：%s', par.casePreset);
    end

    % 按当前预设重新闭合无量纲群
    par = closeDimensionlessGroups_local(par);

    % 基线工况覆盖
    if par.runBaselineNoPorous
        par.Pi_p           = 0.0;
        par.enableSeepage  = false;
        par.enableSkeleton = false;
    end
end

%% =========================================================================
function par = applyDefaultSkeletonParams_local(par)

    if ~isfield(par,'tolSkel'), par.tolSkel = 1.0e-6; end
    if ~isfield(par,'maxSkel'), par.maxSkel = 120; end
    if ~isfield(par,'loadGain'), par.loadGain = 0.10; end

    if ~isfield(par,'H00min'), par.H00min = -0.20; end
    if ~isfield(par,'H00max'), par.H00max = 5.0;   end

    if ~isfield(par,'omegaP'), par.omegaP = 0.70; end
    if ~isfield(par,'tolRe'), par.tolRe = 1.0e-5; end
    if ~isfield(par,'maxRe'), par.maxRe = 60; end
    if ~isfield(par,'tolOuter'), par.tolOuter = 1.0e-5; end
    if ~isfield(par,'maxOuter'), par.maxOuter = 80; end
    if ~isfield(par,'Hp'), par.Hp = 1.0; end
    if ~isfield(par,'tolLoad'), par.tolLoad = 1.0e-5; end

    if ~isfield(par,'useSecantH00'), par.useSecantH00 = true; end
    if ~isfield(par,'secantDamping'), par.secantDamping = 0.80; end
    if ~isfield(par,'secantTolDen'), par.secantTolDen = 1.0e-12; end
    if ~isfield(par,'maxH00Step'), par.maxH00Step = 0.25; end
    if ~isfield(par,'maxPostH00'), par.maxPostH00 = 12; end
    if ~isfield(par,'H00_prev_secant'), par.H00_prev_secant = NaN; end
    if ~isfield(par,'F_prev_secant'), par.F_prev_secant = NaN; end

    if ~isfield(par,'seepUpdateKappa'), par.seepUpdateKappa = true; end
    if ~isfield(par,'seepKappaUpdateEvery'), par.seepKappaUpdateEvery = 5; end

    if ~isfield(par,'enableEHLDeflection'), par.enableEHLDeflection = true; end
    if ~isfield(par,'ehlDeflectionWeight'), par.ehlDeflectionWeight = 1.00; end
    if ~isfield(par,'porousDispWeight'),    par.porousDispWeight    = 0.20; end

    if ~isfield(par,'elasticCoeffMode'),      par.elasticCoeffMode = 'physical'; end
    if ~isfield(par,'elasticInfluenceCoeff'), par.elasticInfluenceCoeff = NaN;   end
    if ~isfield(par,'elasticCoreEpsFactor'),  par.elasticCoreEpsFactor  = 0.5;   end

    if ~isfield(par,'deflRelax'),   par.deflRelax   = 0.15; end
    if ~isfield(par,'deltaEHLMax'), par.deltaEHLMax = 2.0;  end
end

%% =========================================================================
function grid2D = buildGrid2D_local(par)
    grid2D.X  = linspace(par.Xmin, par.Xmax, par.Nx).';
    grid2D.Y  = linspace(par.Ymin, par.Ymax, par.Ny).';
    grid2D.dX = grid2D.X(2) - grid2D.X(1);
    grid2D.dY = grid2D.Y(2) - grid2D.Y(1);
    [grid2D.XX, grid2D.YY] = ndgrid(grid2D.X, grid2D.Y);
end

%% =========================================================================
function grid3D = buildGrid3D_local(par)
    grid3D.X   = linspace(par.Xmin, par.Xmax, par.Nx).';
    grid3D.Y   = linspace(par.Ymin, par.Ymax, par.Ny).';
    grid3D.Z   = linspace(0, -par.Hp, par.Nz).';
    grid3D.dX  = grid3D.X(2) - grid3D.X(1);
    grid3D.dY  = grid3D.Y(2) - grid3D.Y(1);
    grid3D.dZ  = abs(grid3D.Z(2) - grid3D.Z(1));
    grid3D.dZs = grid3D.Z(2) - grid3D.Z(1);
end

%% =========================================================================
function par = prepareElasticInfluence_local(par, grid2D)
% A2: 预计算点接触 EHL 表面弹性影响核

    Nx = par.Nx;
    Ny = par.Ny;

    dX = grid2D.dX;
    dY = grid2D.dY;

    switch lower(strtrim(par.elasticCoeffMode))
        case 'physical'
            Eprime = par.E_dim / (1 - par.nu_dim^2);
            if ~isfinite(Eprime) || Eprime <= 0
                error('Eprime 非法，请检查 E_dim 和 nu_dim。');
            end
            par.elasticInfluenceCoeff = 2 * par.pH * par.R_dim / (pi * Eprime * par.b_dim);

        case 'legacy'
            par.elasticInfluenceCoeff = 2 / pi^2;

        otherwise
            error('未知 elasticCoeffMode：%s', par.elasticCoeffMode);
    end

    fprintf('A2 弹性影响系数 C_e = %.8e (mode = %s)\n', ...
            par.elasticInfluenceCoeff, par.elasticCoeffMode);

    xk = (-(Nx-1):(Nx-1)) * dX;
    yk = (-(Ny-1):(Ny-1)) * dY;
    [XXk, YYk] = ndgrid(xk, yk);

    epsReg = par.elasticCoreEpsFactor * sqrt(dX^2 + dY^2);
    Rk = sqrt(XXk.^2 + YYk.^2 + epsReg^2);

    K = par.elasticInfluenceCoeff * dX * dY ./ Rk;

    padNx = 2^nextpow2(3*Nx - 2);
    padNy = 2^nextpow2(3*Ny - 2);

    par.elasticPadNx = padNx;
    par.elasticPadNy = padNy;

    par.elasticKernelFFT = fft2(K, padNx, padNy);

    par.elasticIx1 = Nx;
    par.elasticIx2 = 2*Nx - 1;
    par.elasticIy1 = Ny;
    par.elasticIy2 = 2*Ny - 1;

    par.elasticKernelCenter = K(Nx, Ny);
end

%% =========================================================================
function field2D = initField2D_local(par, grid2D)
    field2D.P           = zeros(par.Nx, par.Ny);
    field2D.H           = zeros(par.Nx, par.Ny);

    field2D.rho         = ones(par.Nx, par.Ny);
    field2D.eta         = ones(par.Nx, par.Ny);
    field2D.Xi          = ones(par.Nx, par.Ny);
    field2D.S           = zeros(par.Nx, par.Ny);

    field2D.Uz_surf     = zeros(par.Nx, par.Ny);
    field2D.deltaEHL    = zeros(par.Nx, par.Ny);
    field2D.deltaPorous = zeros(par.Nx, par.Ny);

    field2D.dPp_dZ      = zeros(par.Nx, par.Ny);
    field2D.Hgeom       = 0.5 * (grid2D.XX.^2 + grid2D.YY.^2);

    field2D.cavCountRe  = 0;
    field2D.cavFracRe   = 0.0;
end

%% =========================================================================
function field3D = initField3D_local(par)
    field3D.Pp    = zeros(par.Nx, par.Ny, par.Nz);
    field3D.Ux    = zeros(par.Nx, par.Ny, par.Nz);
    field3D.Uy    = zeros(par.Nx, par.Ny, par.Nz);
    field3D.Uz    = zeros(par.Nx, par.Ny, par.Nz);
    field3D.Ev    = zeros(par.Nx, par.Ny, par.Nz);
    field3D.K     = ones(par.Nx, par.Ny, par.Nz);
    field3D.etaP  = ones(par.Nx, par.Ny, par.Nz);
    field3D.kappa = ones(par.Nx, par.Ny, par.Nz);
end

%% =========================================================================
function P = initPressureHertz_local(grid2D)
    R2 = 1 - grid2D.XX.^2 - grid2D.YY.^2;
    P  = zeros(size(R2));
    idx = (R2 > 0);
    P(idx) = sqrt(R2(idx));
end

%% =========================================================================
function field2D = solveReynolds2D_local(par, grid2D, field2D, field3D)

    Nx = par.Nx;
    Ny = par.Ny;

    cavCountTotal = 0;

    for it = 1:par.maxRe

        Pold = field2D.P;

        field2D.S = buildReynoldsSource_local(par, grid2D, field2D, field3D);

        if par.useMG
            field2D.P = reynoldsVcycle_local( ...
                field2D.P, field2D.Xi, field2D.S, ...
                grid2D.dX, grid2D.dY, par.mgLevel, par, false);
        else
            field2D.P = reynoldsSmoothPressureGS_local( ...
                field2D.P, field2D.Xi, field2D.S, ...
                grid2D.dX, grid2D.dY, 4, par);
        end

        field2D.P(~isfinite(field2D.P)) = 0.0;

        negMask = (field2D.P < 0);
        cavCountTotal = cavCountTotal + nnz(negMask);
        field2D.P(negMask) = 0.0;

        field2D.P(1,:)  = 0.0;
        field2D.P(Nx,:) = 0.0;
        field2D.P(:,1)  = 0.0;
        field2D.P(:,Ny) = 0.0;

        err = max(abs(field2D.P(:) - Pold(:)));
        if err < par.tolRe
            break;
        end
    end

    %--------------------------------------------------------------
    % Y方向镜像对称化（抑制顺序扫描导致的横向偏置）
    %--------------------------------------------------------------
    field2D.P = enforceYSymmetry2D_local(field2D.P);

    field2D.P(~isfinite(field2D.P)) = 0.0;
    field2D.P(field2D.P < 0) = 0.0;

    field2D.P(1,:)  = 0.0;
    field2D.P(Nx,:) = 0.0;
    field2D.P(:,1)  = 0.0;
    field2D.P(:,Ny) = 0.0;

    field2D.cavCountRe = cavCountTotal;
    field2D.cavFracRe  = cavCountTotal / max(1, par.maxRe * par.Nx * par.Ny);
end

%% =========================================================================
function S = buildReynoldsSource_local(par, grid2D, field2D, field3D)
% Reynolds 右端源项：
%   S = d(rho_bar*H)/dX - Pi_p * rho_bar * kappa * dP*/dZ

    Nx = par.Nx;
    Ny = par.Ny;
    S  = zeros(Nx, Ny);

    RH = field2D.rho .* field2D.H;

    for j = 2:Ny-1
        for i = 2:Nx-1
            convTerm = (RH(i,j) - RH(i-1,j)) / grid2D.dX;
            seepTerm = par.Pi_p * field2D.rho(i,j) * field3D.kappa(i,j,1) * field2D.dPp_dZ(i,j);

            val = convTerm - seepTerm;
            if ~isfinite(val)
                val = 0.0;
            end
            S(i,j) = min(max(val, -par.sourceMax), par.sourceMax);
        end
    end
end

%% =========================================================================
function U = reynoldsVcycle_local(U, Xi, S, dX, dY, level, par, isErrorProblem)

    if isErrorProblem
        U = reynoldsSmoothErrorGS_local(U, Xi, S, dX, dY, par.mgPre, par);
    else
        U = reynoldsSmoothPressureGS_local(U, Xi, S, dX, dY, par.mgPre, par);
    end

    if level <= 1 || min(size(U)) < 7
        if isErrorProblem
            U = reynoldsSmoothErrorGS_local(U, Xi, S, dX, dY, 12, par);
        else
            U = reynoldsSmoothPressureGS_local(U, Xi, S, dX, dY, 12, par);
        end
        return;
    end

    R = S - reynoldsApplyA_local(U, Xi, dX, dY);

    Rc  = restrictFW_local(R);
    Xic = max(restrictFW_local(Xi), par.XiMin);

    Ec = zeros(size(Rc));
    Ec = reynoldsVcycle_local(Ec, Xic, Rc, 2*dX, 2*dY, level-1, par, true);

    E = prolongBilinear_local(Ec, size(U));
    U = U + E;

    if isErrorProblem
        U = reynoldsSmoothErrorGS_local(U, Xi, S, dX, dY, par.mgPost, par);
    else
        U = reynoldsSmoothPressureGS_local(U, Xi, S, dX, dY, par.mgPost, par);
    end
end

%% =========================================================================
function P = reynoldsSmoothPressureGS_local(P, Xi, S, dX, dY, nSweep, par)

    [Nx, Ny] = size(P);

    for it = 1:nSweep
        for j = 2:Ny-1
            for i = 2:Nx-1

                Xi_e = 0.5 * (Xi(i+1,j) + Xi(i,j));
                Xi_w = 0.5 * (Xi(i-1,j) + Xi(i,j));
                Xi_n = 0.5 * (Xi(i,j+1) + Xi(i,j));
                Xi_s = 0.5 * (Xi(i,j-1) + Xi(i,j));

                Ae = Xi_e / dX^2;
                Aw = Xi_w / dX^2;
                An = Xi_n / dY^2;
                As = Xi_s / dY^2;
                Ap = Ae + Aw + An + As;

                if ~isfinite(Ap) || Ap <= 0
                    Ap = 1.0;
                end

                rhs = Ae * P(i+1,j) + Aw * P(i-1,j) + ...
                      An * P(i,j+1) + As * P(i,j-1) - S(i,j);

                pNew   = rhs / Ap;
                P(i,j) = (1 - par.omegaP) * P(i,j) + par.omegaP * pNew;
            end
        end

        P(~isfinite(P)) = 0.0;
        P(P < 0) = 0.0;
        P(1,:)   = 0.0;
        P(end,:) = 0.0;
        P(:,1)   = 0.0;
        P(:,end) = 0.0;
    end
end

%% =========================================================================
function E = reynoldsSmoothErrorGS_local(E, Xi, S, dX, dY, nSweep, par)

    [Nx, Ny] = size(E);

    for it = 1:nSweep
        for j = 2:Ny-1
            for i = 2:Nx-1

                Xi_e = 0.5 * (Xi(i+1,j) + Xi(i,j));
                Xi_w = 0.5 * (Xi(i-1,j) + Xi(i,j));
                Xi_n = 0.5 * (Xi(i,j+1) + Xi(i,j));
                Xi_s = 0.5 * (Xi(i,j-1) + Xi(i,j));

                Ae = Xi_e / dX^2;
                Aw = Xi_w / dX^2;
                An = Xi_n / dY^2;
                As = Xi_s / dY^2;
                Ap = Ae + Aw + An + As;

                if ~isfinite(Ap) || Ap <= 0
                    Ap = 1.0;
                end

                rhs = Ae * E(i+1,j) + Aw * E(i-1,j) + ...
                      An * E(i,j+1) + As * E(i,j-1) - S(i,j);

                eNew   = rhs / Ap;
                E(i,j) = (1 - par.omegaP) * E(i,j) + par.omegaP * eNew;
            end
        end

        E(~isfinite(E)) = 0.0;
        E(1,:)   = 0.0;
        E(end,:) = 0.0;
        E(:,1)   = 0.0;
        E(:,end) = 0.0;
    end
end

%% =========================================================================
function AP = reynoldsApplyA_local(P, Xi, dX, dY)

    [Nx, Ny] = size(P);
    AP = zeros(Nx, Ny);

    for j = 2:Ny-1
        for i = 2:Nx-1

            Xi_e = 0.5 * (Xi(i+1,j) + Xi(i,j));
            Xi_w = 0.5 * (Xi(i-1,j) + Xi(i,j));
            Xi_n = 0.5 * (Xi(i,j+1) + Xi(i,j));
            Xi_s = 0.5 * (Xi(i,j-1) + Xi(i,j));

            Ae = Xi_e / dX^2;
            Aw = Xi_w / dX^2;
            An = Xi_n / dY^2;
            As = Xi_s / dY^2;
            Ap = Ae + Aw + An + As;

            AP(i,j) = Ae * P(i+1,j) + Aw * P(i-1,j) + ...
                      An * P(i,j+1) + As * P(i,j-1) - Ap * P(i,j);
        end
    end
end

%% =========================================================================
function C = restrictFW_local(F)

    [Nx, Ny] = size(F);
    Nxc = floor((Nx + 1) / 2);
    Nyc = floor((Ny + 1) / 2);

    C = zeros(Nxc, Nyc);

    for jc = 2:Nyc-1
        for ic = 2:Nxc-1
            i = 2*ic - 1;
            j = 2*jc - 1;

            C(ic,jc) = 0.25   * F(i,j) + ...
                       0.125  * (F(i-1,j) + F(i+1,j) + F(i,j-1) + F(i,j+1)) + ...
                       0.0625 * (F(i-1,j-1) + F(i-1,j+1) + F(i+1,j-1) + F(i+1,j+1));
        end
    end
end

%% =========================================================================
function F = prolongBilinear_local(C, targetSize)

    Nxf = targetSize(1);
    Nyf = targetSize(2);
    F = zeros(Nxf, Nyf);

    [Nxc, Nyc] = size(C);

    for jc = 1:Nyc
        for ic = 1:Nxc
            i = 2*ic - 1;
            j = 2*jc - 1;
            if i <= Nxf && j <= Nyf
                F(i,j) = C(ic,jc);
            end
        end
    end

    for j = 1:2:Nyf
        for i = 2:2:Nxf-1
            F(i,j) = 0.5 * (F(i-1,j) + F(i+1,j));
        end
    end

    for j = 2:2:Nyf-1
        for i = 1:2:Nxf
            F(i,j) = 0.5 * (F(i,j-1) + F(i,j+1));
        end
    end

    for j = 2:2:Nyf-1
        for i = 2:2:Nxf-1
            F(i,j) = 0.25 * (F(i-1,j-1) + F(i-1,j+1) + F(i+1,j-1) + F(i+1,j+1));
        end
    end

    F(1,:)   = 0.0;
    F(end,:) = 0.0;
    F(:,1)   = 0.0;
    F(:,end) = 0.0;
end

%% =========================================================================
function loadNow = computeLoad2D_local(grid2D, P)
    tmp = trapz(grid2D.Y.', P, 2);
    loadNow = trapz(grid2D.X, tmp);
end

%% =========================================================================
function par = adjustH00_local(par, loadNow)
% 带阻尼 Secant 法更新 H00

    Hcur = par.H00;
    Fcur = loadNow - par.Wtarget;

    useSecant = false;
    if par.useSecantH00
        if isfinite(par.H00_prev_secant) && isfinite(par.F_prev_secant)
            denom = (Fcur - par.F_prev_secant);
            if abs(denom) > par.secantTolDen && abs(Hcur - par.H00_prev_secant) > eps
                useSecant = true;
            end
        end
    end

    if useSecant
        dH = -Fcur * (Hcur - par.H00_prev_secant) / (Fcur - par.F_prev_secant);
        dH = par.secantDamping * dH;
    else
        dH = par.loadGain * Fcur;
    end

    if ~isfinite(dH)
        dH = par.loadGain * Fcur;
    end

    dH = min(max(dH, -par.maxH00Step), par.maxH00Step);

    Hnew = Hcur + dH;
    Hnew = min(max(Hnew, par.H00min), par.H00max);

    par.H00_prev_secant = Hcur;
    par.F_prev_secant   = Fcur;
    par.H00 = Hnew;
end

%% =========================================================================
function plotResults_local(grid2D, field2D, field3D)

    XX = grid2D.XX;
    YY = grid2D.YY;
    X  = grid2D.X;
    Y  = grid2D.Y;

    [~, j0] = min(abs(Y));
    xLine   = X;
    pLine   = field2D.P(:, j0);
    hLine   = field2D.H(:, j0);

    figure('Name','Pressure Surface');
    surf(XX, YY, field2D.P);
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('P');
    title('无量纲压力分布 P(X,Y)');

    figure('Name','Pressure Contour');
    safeContour_local(XX, YY, field2D.P); colorbar;
    xlabel('X'); ylabel('Y');
    title('无量纲压力等值线');

    figure('Name','Film Thickness Surface');
    surf(XX, YY, field2D.H);
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('H');
    title('无量纲膜厚分布 H(X,Y)');

    figure('Name','Film Thickness Contour');
    safeContour_local(XX, YY, field2D.H); colorbar;
    xlabel('X'); ylabel('Y');
    title('无量纲膜厚等值线');

    figure('Name','EHL Deflection Surface');
    surf(XX, YY, field2D.deltaEHL);
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('\\delta_{EHL}');
    title('主 EHL 非局部表面弹性变形');

    figure('Name','Porous Extra Deflection Surface');
    surf(XX, YY, field2D.deltaPorous);
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('\\delta_{porous}');
    title('多孔附加修正项');

    figure('Name','Surface Uz');
    surf(XX, YY, field3D.Uz(:,:,1));
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('U_z');
    title('多孔表面法向位移 U_z(X,Y,0)');

    figure('Name','Surface Uz Contour');
    safeContour_local(XX, YY, field3D.Uz(:,:,1)); colorbar;
    xlabel('X'); ylabel('Y');
    title('多孔表面法向位移等值线');

    figure('Name','Interface dPp_dZ');
    safeContour_local(XX, YY, field2D.dPp_dZ); colorbar;
    xlabel('X'); ylabel('Y');
    title('界面孔压法向梯度 \\partialP^*/\\partialZ');

    figure('Name','Interface Pp');
    surf(XX, YY, field3D.Pp(:,:,1));
    shading interp; colorbar;
    xlabel('X'); ylabel('Y'); zlabel('P^*_{surf}');
    title('界面孔隙压力 P^*(X,Y,0)');

    figure('Name','Centerline Pressure');
    plot(xLine, pLine, 'LineWidth', 1.8);
    grid on;
    xlabel('X');
    ylabel('Pressure P');
    title(sprintf('中心线压力分布（Y = %.4f）', Y(j0)));

    figure('Name','Centerline Film Thickness');
    plot(xLine, hLine, 'LineWidth', 1.8);
    grid on;
    xlabel('X');
    ylabel('Film Thickness H');
    title(sprintf('中心线膜厚分布（Y = %.4f）', Y(j0)));

    figure('Name','Centerline Profiles Summary');
    subplot(2,1,1);
    plot(xLine, pLine, 'LineWidth', 1.8);
    grid on;
    xlabel('X');
    ylabel('Pressure P');
    title(sprintf('中心线压力分布（Y = %.4f）', Y(j0)));

    subplot(2,1,2);
    plot(xLine, hLine, 'LineWidth', 1.8);
    grid on;
    xlabel('X');
    ylabel('Film Thickness H');
    title(sprintf('中心线膜厚分布（Y = %.4f）', Y(j0)));

    [~, i0] = min(abs(X));
    fprintf('\n中心截面输出信息：\n');
    fprintf('中心线取 Y = %.8f（最接近 0 的网格点）\n', Y(j0));
    fprintf('中心点取 X = %.8f, Y = %.8f\n', X(i0), Y(j0));
    fprintf('中心点压力 P(0,0) ≈ %.8e\n', field2D.P(i0, j0));
    fprintf('中心点膜厚 H(0,0) ≈ %.8e\n', field2D.H(i0, j0));
    fprintf('中心点主 EHL 变形 deltaEHL(0,0) ≈ %.8e\n', field2D.deltaEHL(i0, j0));
    fprintf('中心点多孔附加项 deltaPorous(0,0) ≈ %.8e\n', field2D.deltaPorous(i0, j0));
end

%% =========================================================================
function [field2D, field3D] = updateFilmAndProps_local(par, grid2D, grid3D, field2D, field3D)

    if ~isfield(field2D, 'deltaEHL') || isempty(field2D.deltaEHL)
        field2D.deltaEHL = zeros(par.Nx, par.Ny);
    end
    if ~isfield(field2D, 'deltaPorous') || isempty(field2D.deltaPorous)
        field2D.deltaPorous = zeros(par.Nx, par.Ny);
    end
    if ~isfield(field2D, 'Uz_surf') || isempty(field2D.Uz_surf)
        field2D.Uz_surf = zeros(par.Nx, par.Ny);
    end

    P2 = max(field2D.P, 0.0);

    if par.enableEHLDeflection
        deltaEHL_raw = par.ehlDeflectionWeight * ...
            solveElasticSurfaceDeflection_local(par, grid2D, P2);

        deltaEHL_raw(~isfinite(deltaEHL_raw)) = 0.0;
        deltaEHL_raw = min(max(deltaEHL_raw, 0.0), par.deltaEHLMax);

        field2D.deltaEHL = (1 - par.deflRelax) * field2D.deltaEHL + ...
                            par.deflRelax  * deltaEHL_raw;
    else
        field2D.deltaEHL = zeros(par.Nx, par.Ny);
    end

    if par.enableSkeleton
        field2D.Uz_surf = field3D.Uz(:,:,1);
        deltaPorousRaw = -field2D.Uz_surf;
        field2D.deltaPorous = par.porousDispWeight * deltaPorousRaw;
    else
        field2D.Uz_surf     = zeros(par.Nx, par.Ny);
        field2D.deltaPorous = zeros(par.Nx, par.Ny);

        field3D.Ux(:) = 0.0;
        field3D.Uy(:) = 0.0;
        field3D.Uz(:) = 0.0;
        field3D.Ev(:) = 0.0;
    end

    field2D.deltaPorous(~isfinite(field2D.deltaPorous)) = 0.0;

    field2D.H = par.H00 + field2D.Hgeom + field2D.deltaEHL + field2D.deltaPorous;
    field2D.H = max(field2D.H, par.Hmin);

    field2D.eta = etaBarFromP_local(par, P2);
    field2D.rho = rhoBarFromP_local(par, P2);

    Pp = max(field3D.Pp, 0.0);
    field3D.etaP = etaBarFromP_local(par, Pp);

    if isempty(field3D.Ev)
        field3D.Ev = zeros(par.Nx, par.Ny, par.Nz);
    end

    field3D.K = exp(par.Mhat * field3D.Ev);
    field3D.K = min(max(field3D.K, par.Kmin), par.Kmax);

    field3D.kappa = field3D.K ./ max(field3D.etaP, par.kappaFloor);
    field3D.kappa = min(max(field3D.kappa, par.kappaFloor), par.kappaMax);

    field2D.Xi = field2D.rho .* (field2D.H.^3) ./ (par.Lambda * field2D.eta);
    field2D.Xi(~isfinite(field2D.Xi)) = par.XiMin;
    field2D.Xi = min(max(field2D.Xi, par.XiMin), par.XiMax);

    if (~par.enableSeepage) || (par.Pi_p == 0)
        field2D.dPp_dZ = zeros(par.Nx, par.Ny);
    else
        if par.useSecondOrderDz && par.Nz >= 3
            field2D.dPp_dZ = (-3*field3D.Pp(:,:,1) + 4*field3D.Pp(:,:,2) - field3D.Pp(:,:,3)) / (2*grid3D.dZs);
        else
            field2D.dPp_dZ = (field3D.Pp(:,:,2) - field3D.Pp(:,:,1)) / grid3D.dZs;
        end
    end

    field2D.dPp_dZ(~isfinite(field2D.dPp_dZ)) = 0.0;
    field2D.dPp_dZ = min(max(field2D.dPp_dZ, -par.gradPpMax), par.gradPpMax);
end


%% =========================================================================
function [field2D, field3D] = updateFilmAndPropsFrozenShape_local(par, ~, grid3D, field2D, field3D, freeze)
% 冻结形态后的更新：
%   - 固定 deltaEHL、deltaPorous、Pp、Uz、Ev、K
%   - 仅随当前 H00 与当前 P 更新 H、rho、eta、Xi、etaP、kappa、dPp_dZ
%
% 用于外层耦合完成后的最终载荷闭合，不再允许形态继续漂移。

    field2D.deltaEHL    = freeze.deltaEHL;
    field2D.deltaPorous = freeze.deltaPorous;
    field2D.Uz_surf     = freeze.Uz(:,:,1);

    field3D.Pp = freeze.Pp;
    field3D.Uz = freeze.Uz;
    field3D.Ev = freeze.Ev;
    field3D.K  = freeze.K;

    % 1) 当前表面压力（非负）
    P2 = max(field2D.P, 0.0);

    % 2) 冻结形态下的膜厚
    field2D.H = par.H00 + field2D.Hgeom + field2D.deltaEHL + field2D.deltaPorous;
    field2D.H = max(field2D.H, par.Hmin);

    % 3) 二维压力相关物性
    field2D.eta = etaBarFromP_local(par, P2);
    field2D.rho = rhoBarFromP_local(par, P2);

    % 4) 三维孔隙流体黏度与渗流系数（Pp、K 固定）
    Pp = max(field3D.Pp, 0.0);
    field3D.etaP = etaBarFromP_local(par, Pp);

    field3D.K = min(max(field3D.K, par.Kmin), par.Kmax);
    field3D.kappa = field3D.K ./ max(field3D.etaP, par.kappaFloor);
    field3D.kappa = min(max(field3D.kappa, par.kappaFloor), par.kappaMax);

    % 5) Reynolds 系数
    field2D.Xi = field2D.rho .* (field2D.H.^3) ./ (par.Lambda * field2D.eta);
    field2D.Xi(~isfinite(field2D.Xi)) = par.XiMin;
    field2D.Xi = min(max(field2D.Xi, par.XiMin), par.XiMax);

    % 6) 固定 Pp 下的界面法向梯度
    if (~par.enableSeepage) || (par.Pi_p == 0)
        field2D.dPp_dZ = zeros(par.Nx, par.Ny);
    else
        if par.useSecondOrderDz && par.Nz >= 3
            field2D.dPp_dZ = (-3*field3D.Pp(:,:,1) + 4*field3D.Pp(:,:,2) - field3D.Pp(:,:,3)) / (2*grid3D.dZs);
        else
            field2D.dPp_dZ = (field3D.Pp(:,:,2) - field3D.Pp(:,:,1)) / grid3D.dZs;
        end
    end

    field2D.dPp_dZ(~isfinite(field2D.dPp_dZ)) = 0.0;
    field2D.dPp_dZ = min(max(field2D.dPp_dZ, -par.gradPpMax), par.gradPpMax);
end


%% =========================================================================
function deltaEHL = solveElasticSurfaceDeflection_local(par, ~, P)

    if ~par.enableEHLDeflection
        deltaEHL = zeros(size(P));
        return;
    end

    P = max(P, 0.0);
    P(~isfinite(P)) = 0.0;

    Ufull = real(ifft2( ...
        fft2(P, par.elasticPadNx, par.elasticPadNy) .* par.elasticKernelFFT ));

    deltaEHL = Ufull(par.elasticIx1:par.elasticIx2, par.elasticIy1:par.elasticIy2);

    deltaEHL(~isfinite(deltaEHL)) = 0.0;
    deltaEHL(deltaEHL < 0) = 0.0;
end

%% =========================================================================
function field3D = solveSeepage3D_local(par, grid3D, field2D, field3D)

    if (~par.enableSeepage) || (par.Pi_p == 0)
        field3D.Pp(:,:,1) = field2D.P;
        field3D.Pp(:,:,2:end) = par.P0_bottom;
        return;
    end

    Nx = par.Nx;
    Ny = par.Ny;
    Nz = par.Nz;

    dX = grid3D.dX;
    dY = grid3D.dY;
    dZ = grid3D.dZ;

    omega = par.omegaSeep;

    Pp = field3D.Pp;

    if par.seepUpdateKappa
        [~, kappa_now] = rebuildSeepageKappa_local(par, field3D.K, Pp);
    else
        
        kappa_now = max(field3D.kappa, par.kappaFloor);
    end

    for it = 1:par.maxSeep
        Pold = Pp;

        if par.seepUpdateKappa
            if it == 1 || mod(it-1, par.seepKappaUpdateEvery) == 0
                [~, kappa_now] = rebuildSeepageKappa_local(par, field3D.K, Pp);
            end
        end

        kappa = max(kappa_now, par.kappaFloor);

        Pp(1,:,2:Nz-1)   = Pp(2,:,2:Nz-1);
        Pp(Nx,:,2:Nz-1)  = Pp(Nx-1,:,2:Nz-1);
        Pp(:,1,2:Nz-1)   = Pp(:,2,2:Nz-1);
        Pp(:,Ny,2:Nz-1)  = Pp(:,Ny-1,2:Nz-1);

        Pp(:,:,1)  = field2D.P;
        Pp(:,:,Nz) = par.P0_bottom;

        for k = 2:Nz-1
            for j = 2:Ny-1
                for i = 2:Nx-1

                    ke = harmonicMean_local(kappa(i+1,j,k), kappa(i,j,k));
                    kw = harmonicMean_local(kappa(i-1,j,k), kappa(i,j,k));
                    kn = harmonicMean_local(kappa(i,j+1,k), kappa(i,j,k));
                    ks = harmonicMean_local(kappa(i,j-1,k), kappa(i,j,k));
                    kt = harmonicMean_local(kappa(i,j,k-1), kappa(i,j,k));
                    kb = harmonicMean_local(kappa(i,j,k+1), kappa(i,j,k));

                    aE = ke / dX^2;
                    aW = kw / dX^2;
                    aN = kn / dY^2;
                    aS = ks / dY^2;
                    aT = kt / dZ^2;
                    aB = kb / dZ^2;
                    aP = aE + aW + aN + aS + aT + aB;

                    rhs = aE * Pp(i+1,j,k) + aW * Pp(i-1,j,k) + ...
                          aN * Pp(i,j+1,k) + aS * Pp(i,j-1,k) + ...
                          aT * Pp(i,j,k-1) + aB * Pp(i,j,k+1);

                    pNew = rhs / max(aP, par.kappaFloor);

                    Pp(i,j,k) = (1 - omega) * Pp(i,j,k) + omega * pNew;
                end
            end
        end

        Pp(1,:,2:Nz-1)   = Pp(2,:,2:Nz-1);
        Pp(Nx,:,2:Nz-1)  = Pp(Nx-1,:,2:Nz-1);
        Pp(:,1,2:Nz-1)   = Pp(:,2,2:Nz-1);
        Pp(:,Ny,2:Nz-1)  = Pp(:,Ny-1,2:Nz-1);

        Pp(:,:,1)  = field2D.P;
        Pp(:,:,Nz) = par.P0_bottom;

        Pp(~isfinite(Pp)) = 0.0;

        err = max(abs(Pp(:) - Pold(:)));
        if err < par.tolSeep
            break;
        end
    end

    %--------------------------------------------------------------
    % Y方向镜像对称化（抑制顺序扫描导致的横向偏置）
    %--------------------------------------------------------------
    Pp = enforceYSymmetry3D_local(Pp);

    % 重新强制边界
    Pp(1,:,2:Nz-1)   = Pp(2,:,2:Nz-1);
    Pp(Nx,:,2:Nz-1)  = Pp(Nx-1,:,2:Nz-1);
    Pp(:,1,2:Nz-1)   = Pp(:,2,2:Nz-1);
    Pp(:,Ny,2:Nz-1)  = Pp(:,Ny-1,2:Nz-1);

    Pp(:,:,1)  = field2D.P;
    Pp(:,:,Nz) = par.P0_bottom;

    [field3D.etaP, field3D.kappa] = rebuildSeepageKappa_local(par, field3D.K, Pp);
    field3D.Pp = Pp;
end

%% =========================================================================
function field3D = solveSkeleton3D_local(par, grid3D, field2D, field3D)

    [Nx, Ny, Nz] = size(field3D.Pp);

    dx  = grid3D.dX;
    dy  = grid3D.dY;
    dz  = grid3D.dZ;
    dZs = grid3D.dZs;

    G  = par.Gamma;
    L2 = par.Lambda_s + 2 * par.Gamma;
    al = par.alpha;

    Uz = field3D.Uz;

    field3D.Ux(:) = 0.0;
    field3D.Uy(:) = 0.0;

    for it = 1:par.maxSkel
        UzOld = Uz;

        Uz(1,:,2:Nz-1)   = Uz(2,:,2:Nz-1);
        Uz(Nx,:,2:Nz-1)  = Uz(Nx-1,:,2:Nz-1);
        Uz(:,1,2:Nz-1)   = Uz(:,2,2:Nz-1);
        Uz(:,Ny,2:Nz-1)  = Uz(:,Ny-1,2:Nz-1);

        Tn = -par.skelTopTractionScale * (1 - al) * field2D.P;
        Uz(:,:,1) = (4*Uz(:,:,2) - Uz(:,:,3) + 2*dz * Tn / max(L2,1e-12)) / 3.0;

        Uz(:,:,Nz) = 0.0;

        for k = 2:Nz-1
            for j = 2:Ny-1
                for i = 2:Nx-1

                    dPdz = (field3D.Pp(i,j,k+1) - field3D.Pp(i,j,k-1)) / (2*dZs);

                    aE = G  / dx^2;
                    aW = G  / dx^2;
                    aN = G  / dy^2;
                    aS = G  / dy^2;
                    aT = L2 / dz^2;
                    aB = L2 / dz^2;
                    aP = aE + aW + aN + aS + aT + aB;

                    rhs = aE * Uz(i+1,j,k) + aW * Uz(i-1,j,k) + ...
                          aN * Uz(i,j+1,k) + aS * Uz(i,j-1,k) + ...
                          aT * Uz(i,j,k-1) + aB * Uz(i,j,k+1) - ...
                          al * dPdz;

                    uNew = rhs / max(aP, 1e-12);

                    Uz(i,j,k) = (1 - par.skelRelax) * Uz(i,j,k) + par.skelRelax * uNew;
                end
            end
        end

        Uz(1,:,2:Nz-1)   = Uz(2,:,2:Nz-1);
        Uz(Nx,:,2:Nz-1)  = Uz(Nx-1,:,2:Nz-1);
        Uz(:,1,2:Nz-1)   = Uz(:,2,2:Nz-1);
        Uz(:,Ny,2:Nz-1)  = Uz(:,Ny-1,2:Nz-1);

        Uz(:,:,1)  = (4*Uz(:,:,2) - Uz(:,:,3) + 2*dz * Tn / max(L2,1e-12)) / 3.0;
        Uz(:,:,Nz) = 0.0;

        Uz(~isfinite(Uz)) = 0.0;

        if max(abs(Uz(:) - UzOld(:))) < par.tolSkel
            break;
        end
    end

    %--------------------------------------------------------------
    % Y方向镜像对称化（抑制顺序扫描导致的横向偏置）
    %--------------------------------------------------------------
    Uz = enforceYSymmetry3D_local(Uz);

    % 再强制一次顶/底边界
    Tn = -par.skelTopTractionScale * (1 - al) * field2D.P;
    Uz(:,:,1)  = (4*Uz(:,:,2) - Uz(:,:,3) + 2*dz * Tn / max(L2,1e-12)) / 3.0;
    Uz(:,:,Nz) = 0.0;

    field3D.Uz = Uz;

    %--------------------------------------------------------------
    % 由对称化后的 Uz 重新计算 Ev
    %--------------------------------------------------------------
    field3D.Ev = computeEvFromUz_local(Uz, dZs);
end

%% =========================================================================
function etaBar = etaBarFromP_local(par, Pbar)

    Pbar = max(Pbar, 0.0);

    etaExp = par.Aeta * ((1 + par.betaEta * Pbar).^par.Zeta - 1);
    etaExp = min(max(etaExp, par.etaExpMin), par.etaExpMax);
    etaBar = exp(etaExp);

    etaBar(~isfinite(etaBar)) = 1.0;
end

%% =========================================================================
function rhoBar = rhoBarFromP_local(par, Pbar)

    Pbar = max(Pbar, 0.0);

    rhoBar = (par.A3 + 1.34 * Pbar) ./ (par.A3 + Pbar);
    rhoBar(~isfinite(rhoBar)) = 1.0;
end

%% =========================================================================
function [etaP, kappa] = rebuildSeepageKappa_local(par, K, Pp)

    Pp = max(Pp, 0.0);
    etaP = etaBarFromP_local(par, Pp);

    kappa = K ./ max(etaP, par.kappaFloor);
    kappa = min(max(kappa, par.kappaFloor), par.kappaMax);

    etaP(~isfinite(etaP))   = 1.0;
    kappa(~isfinite(kappa)) = par.kappaFloor;
end

%% =========================================================================
function val = harmonicMean_local(a, b)
    denom = a + b;
    if abs(denom) < 1e-30
        val = 0.0;
    else
        val = 2 * a * b / denom;
    end
end

%% =========================================================================
function A = enforceYSymmetry2D_local(A)
% 对二维场 A(i,j) 做关于 Y=0 的镜像平均
    A = 0.5 * (A + A(:, end:-1:1));
end

%% =========================================================================
function A = enforceYSymmetry3D_local(A)
% 对三维场 A(i,j,k) 做关于 Y=0 的镜像平均
    A = 0.5 * (A + A(:, end:-1:1, :));
end

%% =========================================================================
function Ev = computeEvFromUz_local(Uz, dZs)
% 由 Uz 重新计算 Ev ≈ dUz/dZ

    Nz = size(Uz, 3);
    Ev = zeros(size(Uz));

    Ev(:,:,2:Nz-1) = (Uz(:,:,3:Nz) - Uz(:,:,1:Nz-2)) / (2*dZs);
    Ev(:,:,1)      = (Uz(:,:,2)    - Uz(:,:,1))      / dZs;
    Ev(:,:,Nz)     = (Uz(:,:,Nz)   - Uz(:,:,Nz-1))   / dZs;

    Ev(~isfinite(Ev)) = 0.0;
end

%% =========================================================================
function safeContour_local(XX, YY, Z)
    Zplot = Z;
    Zplot(~isfinite(Zplot)) = 0.0;
    zmin = min(Zplot(:));
    zmax = max(Zplot(:));

    if abs(zmax - zmin) < 1e-14
        surf(XX, YY, Zplot);
        view(2);
        shading flat;
    else
        contourf(XX, YY, Zplot, 30, 'LineColor', 'none');
    end
end
