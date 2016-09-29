% Plot Simcase Disp for all three methods for different Window Size
addpath Y:\Projects\TurbulentFlame\analysis\src\readimx_v2.0_win64\readimx;
% addpath Z:\Planar_Uncertainty_work\Results\final_plots\applyhatch_pluscolor_bundle\;
addpath Z:\Planar_Uncertainty_work\codes\export_fig\;
addpath Z:\Planar_Uncertainty_work\codes\color_code\;
clear;clc;

simcase='Dia';
printfig=0;

if ispc
    basedir='Z:\Planar_Uncertainty_work\';
else
    basedir='/home/shannon/a/bhattac3/Planar_Uncertainty_work/';
end

outputdir=fullfile(basedir,'Simulated_cases','Tests_01_24_2016','plots','combined');
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

winsize={'32','64'};

dia1=0.5:0.5:5;
utrue=0.3;
vtrue=0.6;

pranabase='PIV_scc_deform_pass4_';

for j=1:length(winsize)
    
    
casename=['final_test_ws_',winsize{j}];
pranadir=fullfile(basedir,'Simulated_cases','Tests_01_24_2016','results',['win',winsize{j}],casename,simcase);

if strcmp(winsize{j},'32')
    davisdir=fullfile(basedir,'uncertainty_tests_Davis','Davis_synthetic_image_test','Diameter','CreateMultiframe','PIV_MP(4x32x32_50%ov)');
else
    davisdir=fullfile(basedir,'uncertainty_tests_Davis','Davis_synthetic_image_test','Diameter','CreateMultiframe','PIV_MP(4x64x64_75%ov)');
end

    for i=1:length(dia1);
        
         %% Load prana solution
        pranasol=load(fullfile(pranadir,[pranabase, num2str(2*i-1,'%05.0f'),'.mat']));
        
       
        Up=pranasol.U(2:end-1,2:end-1);
        Vp=-pranasol.V(2:end-1,2:end-1);
        
        if i==1
            Xp=pranasol.X(2:end-1,2:end-1);
            Yp=pranasol.Y(2:end-1,2:end-1);
        end
        %% Load Davis Solution
        davissol=readimx(fullfile(davisdir,['B', num2str(i,'%05.0f'),'.vc7']));
        
        Utemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        Ud=Utemp(2:end-1,2:end-1);
        Vd=Vtemp(2:end-1,2:end-1);
        
        %Davis X and Y grid
        if i==1
            [D] = create2DVec(davissol.Frames{1,1});
            Xd=D.X(2:end-1,2:end-1);
            Yd=D.Y(2:end-1,2:end-1);
            clear D;             
        end
        
        %% Get True Solution
        
        Utrue=utrue;
        Vtrue=vtrue;
        
        %Interpolated true solution onto Davis Grid
        
        Utrued=utrue;
        Vtrued=vtrue;
        %% Calculate Error
        
        err_prana=[abs(Utrue-Up(:));abs(Vtrue-Vp(:))];
        err_davis=[abs(Utrued-Ud(:));abs(Vtrued-Vd(:))];
        

        %% Get MC uncertainty estimates
        [syp,sxp]=size(pranasol.X);
        Ixx=pranasol.ixx(2:end-1,2:end-1);
        Iyy=pranasol.iyy(2:end-1,2:end-1);
        scaling=pranasol.mi(2:end-1,2:end-1);
        biasx=reshape(pranasol.uncertainty(:,15),sxp,syp);
        biasy=reshape(pranasol.uncertainty(:,16),sxp,syp);
        Autod=reshape(pranasol.uncertainty(:,6),sxp,syp);
        biasx=biasx(2:end-1,2:end-1);
        biasy=biasy(2:end-1,2:end-1);
        Autod=Autod(2:end-1,2:end-1);
        
        %Gradient correction
        Udiff=socdiff(Up,abs(Xp(2,1)-Xp(1,1)),1);
        Vdiff=socdiff(Vp,abs(Yp(1,1)-Yp(2,1)),2);
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        
        %scaling and bias correction
        UMCx=sqrt(biasx.^2+(Ixxt.^2)./scaling);
        UMCy=sqrt(biasy.^2+(Iyyt.^2)./scaling);
        MC=[UMCx;UMCy];
        
        %% Get IM uncertainty estimates
        UIMx=pranasol.imx(2:end-1,2:end-1);
        UIMy=pranasol.imy(2:end-1,2:end-1);
        IM=[UIMx;UIMy];
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
        
        Unrx=Unrxtemp(2:end-1,2:end-1);
        Unry=Unrytemp(2:end-1,2:end-1);
        Unbx=Unbxtemp(2:end-1,2:end-1);
        Unby=Unbytemp(2:end-1,2:end-1);
        
        UCSx=sqrt(Unrx.^2+Unbx.^2);
        UCSy=sqrt(Unry.^2+Unby.^2);
        CS=[UCSx;UCSy];
        %% RMS values
        
        rmsprana(j,i)=rms(err_prana(:));
        rmsdavis(j,i)=rms(err_davis(:));
        rmsMC(j,i)=rms(MC(:));
        rmsIM(j,i)=rms(IM(:));
        rmsCS(j,i)=rms(CS(:));

    end
end

%% Plot
lw=1.2;fs=20;
set(0,'DefaultAxesFontName', 'Times New Roman');
crimson=rgb('crimson');
dblue=rgb('deepskyblue');
blueviolet=rgb('blueviolet');
gray=rgb('gray');
mshape={'+','o'};
msize=[5 5];
lstyle={'-','-'};
figure;set(gcf,'DefaultLineLineWidth',lw);set(gca,'FontSize',fs);hold on;
for j=1:2
    sh=mshape{j};
    ms=msize(j);
    lst=lstyle{j};
%     plot(utrue,rmsprana(j,:),'k-','Marker',sh,'MarkerSize',ms);
%     plot(utrue,rmsdavis(j,:),'k--','Marker',sh,'MarkerSize',ms);
%     plot(utrue,rmsMC(j,:),'r','Marker',sh,'MarkerSize',ms,'color',crimson);
%     plot(utrue,rmsIM(j,:),'b','Marker',sh,'MarkerSize',ms,'color',dblue);
%     plot(utrue,rmsCS(j,:),'c','Marker',sh,'MarkerSize',ms,'color',blueviolet);
    plot(dia1,rmsprana(j,:),'k-','Marker',sh,'Linestyle',lst);
    plot(dia1,rmsdavis(j,:),'k-','Marker',sh,'Linestyle',lst,'color',gray);
    plot(dia1,rmsMC(j,:),'r','Marker',sh,'Linestyle',lst,'color',crimson);
    plot(dia1,rmsIM(j,:),'b','Marker',sh,'Linestyle',lst,'color',dblue);
    plot(dia1,rmsCS(j,:),'c','Marker',sh,'Linestyle',lst','color',blueviolet);
end
hold off;
% [LH1,LH2,LH3,LH4]=legend();
xlabel('Diameter (pix)');
ylabel('Rms Error, Uncertainty (pix)');
title('Variation With Diameter');
lgndstr1={'\sigma^{32}_{prana}','\sigma^{32}_{davis}','\sigma^{32}_{MC}','\sigma^{32}_{PD}','\sigma^{32}_{CS}'};
lgndstr2={'\sigma^{64}_{prana}','\sigma^{64}_{davis}','\sigma^{64}_{MC}','\sigma^{64}_{PD}','\sigma^{64}_{CS}'};
lgndstr={lgndstr1{:} lgndstr2{:}};
lbox=legend(lgndstr,'location','eastoutside','FontSize',12);
axis([0 6 0 0.15]);
axis square;
% end

if printfig==1
    grid on;box on;
    set(gcf,'color','white');
    %             print(gcf,'-dpng',[outputdir,M,'.png'],'-r300');
    
%     export_fig(gcf,fullfile(outputdir,[simcase,'.png']),'-painters');
%     export_fig(gcf,fullfile(outputdir,simcase),'-png','-eps','-painters','-m3','-r300');
%     export_fig(gcf,fullfile(outputdir,'Disp2'),'-png','-pdf','-painters','-m1.5','-r300');
    export_fig(gcf,fullfile(outputdir,simcase),'-png','-pdf','-painters','-m1.5','-r300');
end
