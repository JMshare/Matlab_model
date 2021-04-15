function [] = Aero_static_fits(filenamestatic, filenamecc)

    %% Static
    CL_static_load(filenamestatic);
    CD_static_load(filenamestatic);
    Cm_static_load();
    
    
    %% Static Custer
    CL_Custer_stall_load(filenamecc, filenamestatic);
    CD_Custer_stall_load(filenamecc, filenamestatic);
    Cm_Custer_stall_load();


    %% Plots
    plot_coefs(filenamestatic);


    %%-------------------------------------------------------------------------
    function [CL, Alphas, CLs] = Extend_and_interpolate_CL(Alphas, CLs)
        % theory Stengel Flight Dynamics 2.4.2 p68 LIFT; data Flat Plate https://doi.org/10.3390/en8042438 
        global CLaoa
        %flat plate
        CLaoa = @(A) (2*sind(A).^2.*cosd(A)).*(A>=0) - (2*sind(-A).^2.*cosd(-A)).*(A<0);

        % data merge
        Alphas = [Alphas, 2*Alphas(end)-Alphas(end-1)];
        CLs = [CLs, CLs(end-1)];
        Aoas_prepend = -90:-45;
        Aoas_add = 60:90;
        CLs_prepend = CLaoa(Aoas_prepend);
        CLs_add = CLaoa(Aoas_add);
        Alphas = [Aoas_prepend, Alphas, Aoas_add];
        CLs = [CLs_prepend, CLs, CLs_add];
        CL = @(A) interp1([Alphas], [CLs], A, 'spline', 'extrap');
    end

    function [CD, Alphas, CDs] = Extend_and_interpolate_CD(Alphas, CDs)
        % theory Stengel Flight Dynamics 2.4.2 p77 DRAG; data Flat Plate https://doi.org/10.3390/en8042438 
        global CDaoa
        % flat plate
        CDaoa = @(A) (2*sind(A).^3).*(A>=0) + (2*sind(-A).^3).*(A<0);

        % data merge
        Aoas_add = 60:90;
        CDs_add = CDaoa(Aoas_add);
        Alphas = [Alphas, Aoas_add];
        CDs = [CDs, CDs_add];
        CDplus = @(A) interp1([Alphas], [CDs], A, 'spline', 'extrap');
        CD = @(A) CDplus(A).*(A>=0) + CDplus(-A).*(A<0);
    end

    function [Cm] = Extend_and_interpolate_Cm(CLstatic, CDstatic)
        % Just combines CL and CD acting XCP from the CG. Remove the
        % friction drag though.
        % the CD_static_load and CL_static_load have to be executed already
        % to have them defined in global
        global XCP
        Cmtemp = @(A) XCP*sqrt(CLstatic(A).^2 + (CDstatic(A)-CDstatic(0)).^2);
        Cmtempneg = @(A) -XCP*sqrt(CLstatic(A).^2 + (CDstatic(A)-CDstatic(0)).^2);
        [Across] = fzero(CLstatic, 0);
        Cm = @(A) Cmtemp(A).*(A>=Across) + Cmtempneg(A).*(A<Across);
%         [A0, fval] = fminbnd(@(A) -Cmtemp(A), -20, 20);
%         Cm = @(A) Cmtemp(A).*(A>=A0) - (Cmtemp(A)+2*fval).*(A<A0);
    end




    function [] = CL_static_load(filenamestatic)
        global CLstatic
        A = importdata(filenamestatic.CL);
        A = A.data;
        Alphas = A(:,1)';
        CLs = A(:,2)';
        [CLstatic] = Extend_and_interpolate_CL(Alphas, CLs);
    end

    function [] = CD_static_load(filenamestatic)
        global CDstatic
        A = importdata(filenamestatic.CD);
        A = A.data;
        Alphas = A(:,1)';
        CDs = A(:,2)';
        [CDstatic] = Extend_and_interpolate_CD(Alphas, CDs);
    end

    function [] = Cm_static_load()
        global Cmstatic CLstatic CDstatic
        [Cmstatic] = Extend_and_interpolate_Cm(CLstatic, CDstatic); % the CD_static_load and CL_static_load have to be executed already
    end




    function [] = CL_Custer_stall_load(filenamecc, filenamestatic) 
        global CL_Custer_stall
        global CLstatic

        AS = importdata(filenamecc.stall.CL);
        AS = AS.data;
        Vs = AS(:,1);

        A = importdata(filenamestatic.CL);
        A = A.data;
        Alphas = A(1:end-2,1)';
        Alphas = linspace(Alphas(1),Alphas(end), 100);
        CLs = CLstatic(Alphas);
        
        AoaS = @(v) interp1(Vs, AS(:,2), v, 'linear', 'extrap');
        CLS = @(v) interp1(Vs, AS(:,3), v, 'linear', 'extrap');
        AoaM = @(v) interp1(Vs, AS(:,4), v, 'linear', 'extrap');
        CLM = @(v) interp1(Vs, AS(:,5), v, 'linear', 'extrap');

        Alphas = @(v) [Alphas, AoaS(v), AoaM(v)];
        CLs = @(v) [CLs, CLS(v), CLM(v)];
        CL_Custer = @(v,aoa) feval(Extend_and_interpolate_CL(Alphas(v), CLs(v)), aoa);

        CL_Custer_stall = @(v,aoa) CL_Custer(v,aoa) - CLstatic(aoa);
    end

    function [] = CD_Custer_stall_load(filenamecc, filenamestatic) 
        global CD_Custer_stall
        global CDstatic

        AS = importdata(filenamecc.stall.CD);
        AS = AS.data;
        Vs = AS(:,1);

        A = importdata(filenamestatic.CD);
        A = A.data;
        Alphas = A(1:end-1,1)';
        Alphas = linspace(Alphas(1),Alphas(end), 100);
        CDs = CDstatic(Alphas);

        AoaS = @(v) interp1(Vs, AS(:,2), v, 'linear', 'extrap');
        CDS = @(v) interp1(Vs, AS(:,3), v, 'linear', 'extrap');

        Alphas = @(v) [Alphas, AoaS(v)];
        CDs = @(v) [CDs, CDS(v)];
        CD_Custer = @(v,aoa) feval(Extend_and_interpolate_CD(Alphas(v), CDs(v)), aoa);

        CD_Custer_stall = @(v,aoa) CD_Custer(v,aoa) - CDstatic(aoa);
    end

    function [] = Cm_Custer_stall_load() 
        global Cm_Custer_stall 
        
        % This would be ok in general but we can assume custer changes the
        % XCP such that it doesn't really have any cm effect
%         Cm_Custer_stall = @(v, aoa) feval(get_fun_v(v), aoa);
%         function [Cm_Custer_stall_v] = get_fun_v(v)
%             global CL_Custer_stall CD_Custer_stall
%             CL_Custer_stall_v = @(aoa) CL_Custer_stall(v, aoa);
%             CD_Custer_stall_v = @(aoa) CD_Custer_stall(v, aoa);
%             Cm_Custer_v = @(aoa) feval(Extend_and_interpolate_Cm(CL_Custer_stall_v, CD_Custer_stall_v), aoa);
%             Cm_Custer_stall_v = @(aoa) Cm_Custer_v(aoa);
%         end
        
        global Cmstatic
        Cm_Custer_stall = @(v, aoa) 0*Cmstatic(aoa);
    end



    function [] = plot_coefs(filenamestatic)

        %% CL
        global CLstatic CL_Custer_stall CLaoa
        A = importdata(filenamestatic.CL);
        A = A.data;
        Alphas = A(:,1)';
        CLs = A(:,2)';
        % Alphas = [Alphas, 2*Alphas(end)-Alphas(end-1)];
        % CLs = [CLs, CLs(end-1)];
        figure(200);
        clf
        subplot(1,3,1);
        plot(Alphas, CLs, '*r', 'MarkerSize', 10, 'DisplayName', 'measured');
        hold on;
        A0 = -10;
        As = linspace(A0, 90, 500);
        plot(As, CLstatic(As), 'b-', 'DisplayName', 'extrapolated');
        xlabel('aoa [deg]');
        ylabel('CL');
        plot(As, CLaoa(As), 'b--', 'DisplayName', 'flat plate');
        title('Lift coef model');
        grid minor;
        
        hold on;
        Asp = As(1:10:end);
        plot(As, CL_Custer_stall(5, As) + CLstatic(As), 'k--', 'DisplayName', sprintf('stall-delayed V=%dm/s', 5));
        hold on;
        plot(As, CL_Custer_stall(10, As) + CLstatic(As), 'k', 'DisplayName', sprintf('stall-delayed V=%dm/s', 10));
        hold on;
        plot(Asp, CL_Custer_stall(20, Asp) + CLstatic(Asp), 'k.', 'DisplayName', sprintf('stall-delayed V=%dm/s', 20));
        plot(As, CLstatic(As), 'b-', 'HandleVisibility', 'Off');
        
        legend('location', 'best');
        xlim([A0-1,91]);
        
        %% CD
        global CDstatic CD_Custer_stall CDaoa
        A = importdata(filenamestatic.CD);
        A = A.data;
        Alphas = A(:,1)';
        CDs = A(:,2)';
        figure(200);
        subplot(1,3,2);
        plot(Alphas, CDs, '*r', 'MarkerSize', 10);
        hold on;
        plot(As, CDstatic(As), 'b-');
        xlabel('aoa [deg]');
        ylabel('CD');
        plot(As, CDaoa(As), 'b--');
        title('Drag coef model');
        grid minor;
        xlim([A0-1,91]);

        for v = [5,10,20,30]
            hold on;
            plot(As, CD_Custer_stall(v, As) + CDstatic(As), 'k');
        end
        plot(As, CDstatic(As), 'b-');

        
        %% CM
        global Cmstatic Cm_Custer_stall
        figure(200);
        subplot(1,3,3);
        hold on;
        plot(As, Cmstatic(As), 'b-');
        xlabel('aoa [deg]');
        ylabel('Cm');
        title('Pitch moment coef model');
        grid minor;
        for v = [5,10,20,30]
            hold on;
            plot(As, Cm_Custer_stall(v, As) + Cmstatic(As), 'k');
        end
        plot(As, Cmstatic(As), 'b-');
        xlim([A0-1,91]);
        
        %% for the legend to be on top of all subplots (not in the background of subplot 2..)
        subplot(1,3,1); 
        ax = get(gcf,'children');
        ind = find(isgraphics(ax,'Legend'));
        set(gcf,'children',ax([ind:end,1:ind-1]))
    end

end
