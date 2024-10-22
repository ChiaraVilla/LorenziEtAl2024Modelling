%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%      Plot stripes and cells at final time (macroscopic model)       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_stripes(rho0,rhoA,rhoB,par,x1,x2)

set(gca,'TickLabelInterpreter','latex')
hold on
set(0,'DefaultAxesFontSize',14)
rhomax = max([max(max(rhoA)),max(max(rhoB))]);
%rhomax = 1.0729;

clf
subplot(2,3,1)
Plot_StripeS_Timet(rho0,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,rhomax,'n')
title(['Initial conditions (0h)'])

subplot(2,3,2)
Plot_StripeS_Timet(rhoA,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,rhomax,'n')
title(['A) $\bar{y}_K^L=\bar{y}_K^F=0$ (48h)'])
subplot(2,3,5)
% Plot_StripeS_Timet(rhoA,x1,x2,par,11./4,21./4,-5./4,5./4,rhomax,'z')
% Plot_StripeS_Timet(rhoA,x1,x2,par,11./4,21./4,-2-5./4,-2+5./4,rhomax,'z')
Plot_StripeS_Timet(rhoA,x1,x2,par,11./4,21./4,0-5./4,0+5./4,rhomax,'z')

subplot(2,3,3)
Plot_StripeS_Timet(rhoB,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,rhomax,'n')
title(['B) $\bar{y}_K^L=1$, $\bar{y}_K^F=0$ (48h)'])
subplot(2,3,6)
Plot_StripeS_Timet(rhoB,x1,x2,par,11./4,21./4,2-5./4,2+5./4,rhomax,'z')

end

function Plot_StripeS_Timet(rho,x1,x2,par,x1min,x1max,x2min,x2max,rhomax,zoomin)

    % Laminin and Fibronectin stripes
    L_f=@(s1,s2) par.L0*((s1<par.str1) | (s1 >par.str2)) + 0*((s1>par.str1) & (s1 <par.str2));
    F_f=@(s1,s2) 0*((s1<par.str1) | (s1 >par.str2)) + par.F0*((s1>par.str1) & (s1 <par.str2));
    L = repmat(L_f(x1',x2),1,length(x2));
    F = repmat(F_f(x1',x2),1,length(x2));

    % To better correct numerical approximation of zero
    % rho(rho<=1e-5) = NaN;
    L(rho>0) = NaN;
    F(rho>0) = NaN;

    surf(x1,x2,L','FaceColor','k','FaceAlpha',0.4);
    hold on 
    surf(x1,x2,F','FaceColor','k','FaceAlpha',0.25)
    hold on
    surf(x1,x2,rho')
    view(0,90)
    xlim([x1min,x1max])
    ylim([x2min,x2max])
    grid off
    shading interp
    axis square
    xlabel('$x_1$')
    ylabel('$x_2$')
    colorbar
    caxis([0 rhomax])
    daspect([1 1 1])
    box on
    if zoomin == 'n' % Normal plot
        xticks([])
        colorln = 'black';
        colorln = 'white';
        text(0.35,4.75,20,'LN','Color',colorln,'FontSize',18)  
        text(2.35,4.75,20,'FN','Color',colorln,'FontSize',18)
        text(4.35,4.75,20,'LN','Color',colorln,'FontSize',18)
    elseif zoomin == 'z' % Zoom in
        xticks([x1min,x1max])
        colorln = 'black';
        colorln = 'white';
        x2_pos = (x2min+x2max)/2.0;
        text(3.25,x2_pos+0.75,20,'FN','Color',colorln,'FontSize',18)
        text(4.4,x2_pos+0.75,20,'LN','Color',colorln,'FontSize',18)
    end
    colorstripeb = 'black';
    colorstripeb = [.7 .7 .7];
    lineconnect = line('XData',[x1min,x1max], 'YData',[0,0], 'ZData',[20,20], ...
        'Color',colorstripeb,'LineWidth',2,'LineStyle','--');
    lineconnect = line('XData',[par.str1,par.str1], 'YData',[0,x2max], 'ZData',[20,20], ...
        'Color',colorstripeb,'LineWidth',2,'LineStyle','--');
    lineconnect = line('XData',[par.str2,par.str2], 'YData',[0,x2max], 'ZData',[20,20], ...
        'Color',colorstripeb,'LineWidth',2,'LineStyle','--');
    lineconnect = line('XData',[par.str1,par.str1], 'YData',[x2min,0], 'ZData',[20,20], ...
        'Color',[.7 .7 .7],'LineWidth',2,'LineStyle','--');
    lineconnect = line('XData',[par.str2,par.str2], 'YData',[x2min,0], 'ZData',[20,20], ...
        'Color',[.7 .7 .7],'LineWidth',2,'LineStyle','--');
    drawnow


end