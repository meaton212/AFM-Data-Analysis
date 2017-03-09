clear all; close all
% filenames to analye (should be in the 'data' directory')
[filenames1] = dir(fullfile('Matlab Example Files','*.mat'));
source_directory2 = pwd;
TestFileDirectory2 = strcat(source_directory2,'\Matlab Example Files');
AdhesionDirectory=strcat(source_directory2,'\Matlab Example Files\Adhesion\files\');
[filenames_adh] = dir(fullfile('Matlab Example Files\Adhesion\files','*.mat'));
[m,n] = size(filenames1);
[o,p] = size(filenames_adh);
fileID = fopen('VEData.txt','w');
formatSpec='%6s \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \n';
fprintf(fileID,formatSpec,'Temperature','Name','Loss Tangent','Std Loss Tan','E''','Standard Deviation','E"','StD','contact radius (nm)','E*a')

set(0,'defaultaxesfontsize', 16)
set(0,'defaultlinelinewidth',1)

for filenum=1:m
    clear TanDelta; clear StorageMod; clear LossMod; clear contact_radius; clear Fadh
    clear pp; clear tt; clear dd; clear sined; clear sinep
    filenames=filenames1(filenum).name;
    filenames=strrep(filenames,'.mat','');
    load([TestFileDirectory2 '\' filenames '.mat']) % read in all of the data.
    area=strsplit(filenames,{'_','-'});
    temp=strsplit(char(area(4)),{'pt','C'},'CollapseDelimiters',true);
    temp(1)=strrep(temp(1), 'm', '-');
    temp2=strjoin({char(temp(1)),char(temp(2))},'.');
    
    for i=1:o
        adh_area=strsplit(filenames_adh(i).name,{'_','-'});
        if strcmp(adh_area(1),area(1))==1  & strcmp(adh_area(4),area(4))==1 & strcmp(adh_area(3),area(3))==1
            adh_file=strrep(filenames_adh(i).name,'.mat','');
            load([AdhesionDirectory '\' adh_file '.mat']); 
        end
    end    
     
%     Fadh1=mean(Fadh);
    
    clear NSMU % NSMatUtilities - not generally useful here
    period=(freq^-1)*10^6;
    delt=period/length(t_loadcurve.p1);
    delete(h)
    % define a generalized sine fit that we'll use to fit things
    sinefit=@(C,t) C(1)+C(2)*sin(2*pi.*t/period + C(3));
    
    pixelstocalc=1:numel( fieldnames(t_loadcurve)); % chage this if we dont want to calculate everything
    
    hwait = waitbar(0,['Calculating for file ' num2str(filenum) ' of ' num2str(length(filenames1))]);

    pixiter=0;
    for pixel=pixelstocalc  % pixels to include in analysis
        pixiter=pixiter+1;
        waitbar(pixiter / length(pixelstocalc))
        pstring=['p' num2str(pixel)];
        t=t_loadcurve.(pstring); % time from raw load-time curve
        p=p_loadcurve_v.(pstring)*Defl_Sens*spring_constant;
        p_trace=p_trace_v.(pstring)*Defl_Sens*spring_constant;
        p_retrace=p_retrace_v.(pstring)*Defl_Sens*spring_constant;
        dpiezo=vertcat(flipud(d_trace.(pstring)), d_retrace.(pstring));  % combine trace and retrace data into single curve
        ppiezo=vertcat(flipud(p_trace), p_retrace);
        d=dpiezo-ppiezo/spring_constant;
        
        
        % now generate initial guesses for C for our fits to dpiezo, d and p
        
        Cdpiezoguess=[0.5*(max(dpiezo)+min(dpiezo)), 0.5*(max(dpiezo)-min(dpiezo)), -pi/2];
        Cdguess=[0.5*(max(d)+min(d)), 0.5*(max(d)-min(d)), -pi/2];
        Cpguess=[0.5*(max(p)+min(p)), 0.5*(max(p)-min(p)), -pi/2];
        
        % fit the displacement to a sign wave and plot the fit
        opts = optimset('Display','off');
        dfit{pixel}=lsqcurvefit(sinefit,Cdguess,t,d,[],[], opts);
        pfit{pixel}=lsqcurvefit(sinefit,Cpguess,t,p,[],[], opts);
        phi(pixel)=180*(pfit{pixel}(3)-dfit{pixel}(3))/pi; % phase angle from sinusoidal fits
        phi2(pixel)=asind(trapz(d,p)/(pi*pfit{pixel}(2)*dfit{pixel}(2))); % phase angle from area of hysteresis loop
        Era(pixel)=(3/8)*pfit{pixel}(2)/dfit{pixel}(2);
        TanDelta(pixel)=tand(phi(pixel));
        
        P1=max(p)+2*Fadh(pixel)+2*(Fadh(pixel)*(Fadh(pixel)+max(p)))^0.5;
%         P1=max(p)+2*Fadh1+2*(Fadh1*(Fadh1+max(p)))^0.5;
        hertz=(tip_radius*max(d))^0.5;
        factor=1/((1-(4/3)*(Fadh(pixel)/P1)^0.5)^0.5);
%         factor=1/((1-(4/3)*(Fadh1/P1)^0.5)^0.5);
        contact_radius(pixel)=hertz*factor;
        
        StorageMod(pixel)=(Era(pixel)/contact_radius(pixel))*cosd(phi(pixel))*1000;
        LossMod(pixel)=(Era(pixel)/contact_radius(pixel))*sind(phi(pixel))*1000;
%   
        if  isreal(StorageMod(pixel))==0
             StorageMod(pixel)=0;
        end
        
        if  isreal(LossMod(pixel))==0
            LossMod(pixel)=0;
        end

        tt.(pstring)=t;
        dd.(pstring)=d;
        pp.(pstring)=p;
        sinep.(pstring)=sinefit(pfit{pixel},t);
        sined.(pstring)=sinefit(dfit{pixel},t);
    end
    delete(hwait)
    close all
    
    % now we make a plot of the data from the last pixel considered
    combinedfig=figure('PaperPosition', [0 0 12 5], 'PaperSize', [12 5]);
    
    pdtplot=subplot(1,2,1);
    yyaxis left
    plot(t, p, 'b+')   
    xlabel('t (\mu s)')
    hold on
    ylabel('P (nN)')
    
    yyaxis right
    plot(t, dpiezo, 'r-')
    
    plot(t, d, 'r+')
    ylabel('\delta (nm)')
    legend('P', '\delta_{piezo}', '\delta')
    
    yyaxis left
    plot(pdtplot,t,sinefit(pfit{pixel},t),'b-')   
    yyaxis right
    plot(pdtplot,t,sinefit(dfit{pixel},t),'r-')
    
    % put the phase information as a title on the left plot
    title(pdtplot, ['\phi = ' num2str(phi(pixel),3) '^\circ (' num2str(phi2(pixel),3) '^\circ)'])
    
    pdplot=subplot(1,2,2);
    trace_d=d_trace.(pstring);
    trace_p=p_trace_v.(pstring)*Defl_Sens*spring_constant;
    trace_d_reversed=flipud(trace_d);
    trace_p=flipud(trace_p);
    
    retrace_d=d_retrace.(pstring);
    retrace_p=p_retrace_v.(pstring)*Defl_Sens*spring_constant;
    
    %plot(trace_d-trace_p/spring_constant, trace_p, 'k+')
    hold on
    plot(trace_d_reversed-trace_p/spring_constant, trace_p, 'b+')  % reverse the trace curve
    plot(retrace_d-retrace_p/spring_constant, retrace_p, 'r+')
    plot(d, p, 'k-') % these are the data from the left plot (pdtplot) 
    
    pdplot.YAxisLocation='right';
    plot(pdplot,sinefit(dfit{pixel},t),sinefit(pfit{pixel},t),'r-')
    xlim([min(d) max(d)])
    
    xlabel('\delta (nm)')
    ylabel('P (nN)')
    legend('trace (reversed)', 'retrace', 'P-\delta data', 'P-\delta fit', 'location', 'best')
       
    % put the contact stiffness on the right plot
    title(['|E_r^*|a = ' num2str(0.5*pfit{pixel}(2)/dfit{pixel}(2),3) ' N/m'])

    % use the suplabel command to label the overal figure according to the
    % data filename (from MATLAB file exchange)
    suplabel([strrep(filenames, '_', '\_') '\_' pstring], 't');
    %tightfig  % from MATLAB file exchange
    
    % maximize the figure on secondary monitor if we have one
    pos = get(0,'MonitorPositions');
    sz = size(pos);
    if (sz(1) > 1)
        combinedfig.Position=pos(2,:);
    else
        combinedfig.Position=pos;
    end
    
    
    print(gcf,[source_directory2 filenames '_' pstring '.svg'], '-dsvg')
    
    figure
    tdelt=reshape(TanDelta,[forVolPixel,forVolPixel])'; 
    x=linspace(0,xlength,forVolPixel);
    contourf(x,x,tdelt)
    c=colorbar;
    title('Loss Tangent')
    
    
    ylabel(c,'tan(\delta)')
    ylabel('Distance (nm)')
    xlabel('Distance (nm)')
    savefig(strcat(filenames,'_LossTan.fig'))
    figure
    smod=reshape(StorageMod,[forVolPixel,forVolPixel])'; 
    
    contourf(x,x,smod)
    c=colorbar;
    title('Storage Modulus')
    
    ylabel(c,'E'' (MPa)')
    ylabel('Distance (nm)')
    xlabel('Distance (nm)')
     savefig(strcat(filenames,'_StorMod.fig'))
    figure
    lmod=reshape(LossMod,[forVolPixel,forVolPixel])'; 
  
    contourf(x,x,lmod)
    c=colorbar;
    ylabel(c,'E" (MPa)')
    ylabel('Distance (nm)')
    xlabel('Distance (nm)')
    title('Loss Modulus')
    savefig(strcat(filenames,'_LossMod.fig'))
    name=FilenameSet{2};
    % now save the extracted data as a separate .mat file
    outputfile=strcat(source_directory2, filenames, '_solved.mat');
    save(outputfile,'phi', 'phi2', 'Era', 'dfit', 'pfit','TanDelta','StorageMod','LossMod','tdelt','smod','lmod','name','contact_radius','tt','pp','dd','sinep','sined','x')
    
    formatSpec = '%6s \t %6s \t %8.8f \t %8.8f \t %8.8f \t %8.8f \t %8.8f \t %8.8f \t %8.8f \t %8.8f\n';
 fprintf(fileID,formatSpec,temp2,FilenameSet{2},mean(TanDelta),std(TanDelta),mean(StorageMod),std(StorageMod),mean(LossMod),std(LossMod),mean(contact_radius),mean(Era))
   
end


fclose(fileID);



