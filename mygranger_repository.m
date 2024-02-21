close all;
beta_range=[5 30];
figure_path0='/Volumes/hd5/Data/vmpfc/';
cd(figure_path0);
load('granger0_sakura.mat');
granger0_sakura=granger0;
granger11_sakura=granger11;
asym0_sakura=asym0;

addsakura=1;

clear granger0 granger11;
effective_sessions0={
           {'ABLCnet_effective_all','2014-01-17_14-47-15','2014-01-08_16-14-34','2013-09-06_15-11-27','2013-09-13_14-38-48','2014-09-09_16-32-19'};%,'2011-06-22_16-33-16'} %}, ...
   };

clear selectedsessions;
for sessionnum0=1:numel(effective_sessions0)
    foldername0=effective_sessions0{sessionnum0}{1};
    mkdir(foldername0);
    for sessionnum1=2:numel(effective_sessions0{sessionnum0})
        selectedsessions{sessionnum1-1}=effective_sessions0{sessionnum0}{sessionnum1};
    end
    effective1non_effective2=1;
    effective_sessions=selectedsessions;
    
    resample_num=500;
    try
        cd ('C:\Dropbox (MIT)\Data\vmpfc');
        load aloc8
    catch
        cd ('/Volumes/hd6/Data/');
        load aloc8
    end
    fsum
    cd ('/Volumes/hd5/Data/vmpfc/granger/compactdata/');
    a=dir;anum1=0;
    clear popg0;
    for anum0=1:numel(a)
        if strfind(a(anum0).name,'cmpt_granger(')==1
                for xnum0=1:numel(effective_sessions)
                    if strfind(a(anum0).name,effective_sessions{xnum0})==14
                        anum1=anum1+1;
                        load(a(anum0).name);
                        popg0(anum1).granger=granger1;
                        popg0(anum1).name=granger1.filename;
                    end
                end
        end
    end
    cd ('/Volumes/hd5/Data/vmpfc/');
    clear periodmat0 granger0;
    a=datetime('2012/9/1');
    b=datetime(popg0(1).name(9:18));
    if b>a
        monkey0=3;
    else
        monkey0=2;
    end
    cue1precue2all3_loop=1:1;
    for cue1precue2all3=cue1precue2all3_loop
        switch cue1precue2all3
            case 1
                cuename='cue';
            case 2
                cuename='precue';
            case 3
                cuename='allperiod';
        end
        ap1av2all3_loop=3:3;
        for ap1av2all3=ap1av2all3_loop
            switch ap1av2all3
                case 1
                   apavname='ap';
                case 2
                   apavname='av';
                case 3
                   apavname='apav';
            end
            if cue1precue2all3==1
                period00=[1,3];
            else
                period00=[1,2,3];
            end
            clear asym0 granger11 asym2 grangerout1 grangerout asymout;
            for period0=period00
                switch period0
                    case 1
                        periodname='0stim0';
                    case 2
                        periodname='1stim1';
                    case 3
                        periodname='2resid';
                    case 4
                        periodname='3all';
                end
                fsum=region_assignment_lav14_vmpfc_combined(monkey0,fsum);
                
                subplot(2,2,1);
                if fsum.monkey==2
                    mindepth=4;
                    maxdepth=14;
                else
                    mindepth=0;
                    maxdepth=12;
                end
                if fsum.monkey==2
                    xmin=15;xmax=28;ymin=5;ymax=40;zmin=-maxdepth;zmax=-mindepth;cmin=-maxdepth;cmax=-mindepth;
                else
                    xmin=-20;xmax=8;ymin=5;ymax=40;zmin=-maxdepth;zmax=-mindepth;cmin=-maxdepth;cmax=-mindepth;
                end
                A=colormap(jet);
                [a1,a2]=size(A);
                
                caxis([-10 -2]);
                
                clear connectivity0alpha connectivity0beta connectivity0gamma channel0;
                clear granger12alpha granger12beta granger12gamma;
                clear granger21alpha granger21beta granger21gamma;
                clear connectivity0alpha2 connectivity0beta2 connectivity0gamma2;
                clear connectivity0alpha3 connectivity0beta3 connectivity0gamma3;
                clear granger12alpha2 granger12beta2 granger12gamma2;
                clear granger21alpha2 granger21beta2 granger21gamma2;
                clear grangertunebeta12 grangertunebeta21;
                clear granger21alpha3 granger21beta3 granger21gamma3;
                clear grangerdiffalpha3 grangerdiffbeta3 grangerdiffgamma3;
                clear grangertunealpha3 grangertunebeta3 grangertunegamma3;
                
                xnum0=0;
                for nnum1=1:numel(fsum.freq_channel)
                    for nnum2=1:numel(fsum.freq_channel(nnum1).channel)
                        xnum0=xnum0+1;
                        channel0(xnum0).num=fsum.freq_channel(nnum1).channel(nnum2);
                        channel0(xnum0).name=[fsum.freq_channel(nnum1).group_name,':',num2str(channel0(xnum0).num)];
                        channel0(xnum0).group=nnum1;
                    end
                end
                close;
                figure;
                snum0=0;
                clear Spectra00 Spectra12 Spectradiff freq1 freq00 Group0;
                
                
                numel(popg0)

                for xnum9=1:numel(popg0)
                    channel1=popg0(xnum9).granger.channel1;
                    channel2=popg0(xnum9).granger.channel2;
                    if isempty(channel1)||isempty(channel2)
                        continue
                    end
            
                    nnum1=find([channel0(:).num]==channel1,1);
                    nnum2=find([channel0(:).num]==channel2,1);
            
                    if isempty(nnum1)
                        continue
                    end
                    if isempty(nnum2)
                        continue
                    end
                    connectivity0alpha2(nnum1,nnum2).num=0;
                    connectivity0beta2(nnum1,nnum2).num=0;
                    connectivity0gamma2(nnum1,nnum2).num=0;
                    granger12alpha2(nnum1,nnum2).num=0;
                    granger12beta2(nnum1,nnum2).num=0;
                    granger12gamma2(nnum1,nnum2).num=0;
                    granger21alpha2(nnum1,nnum2).num=0;
                    granger21beta2(nnum1,nnum2).num=0;
                    granger21gamma2(nnum1,nnum2).num=0;
                    grangerdiffalpha2(nnum1,nnum2).num=0;
                    grangerdiffbeta2(nnum1,nnum2).num=0;
                    grangerdiffgamma2(nnum1,nnum2).num=0;
                    grangertunealpha2(nnum1,nnum2).num=0;
                    grangertunebeta2(nnum1,nnum2).num=0;
                    grangertunegamma2(nnum1,nnum2).num=0;
            
                    granger12tunebeta(nnum1,nnum2).num=0;
                    granger21tunebeta(nnum1,nnum2).num=0;
                    granger12tunealpha(nnum1,nnum2).num=0;
                    granger21tunealpha(nnum1,nnum2).num=0;
                    granger12tunegamma(nnum1,nnum2).num=0;
                    granger21tunegamma(nnum1,nnum2).num=0;
                end
                for xnum12=1:7
                    for xnum22=1:7
                        connectivity0alpha3(xnum12,xnum22).num=0;
                        connectivity0beta3(xnum12,xnum22).num=0;
                        connectivity0gamma3(xnum12,xnum22).num=0;
                        granger12alpha3(xnum12,xnum22).num=0;
                        granger12beta3(xnum12,xnum22).num=0;
                        granger12gamma3(xnum12,xnum22).num=0;
                        granger21alpha3(xnum12,xnum22).num=0;
                        granger21beta3(xnum12,xnum22).num=0;
                        granger21gamma3(xnum12,xnum22).num=0;
                        grangerdiffalpha3(xnum12,xnum22).num=0;
                        grangerdiffbeta3(xnum12,xnum22).num=0;
                        granger0(period0).grangerdiffbeta3(xnum12,xnum22).num=0;
                        granger0(period0).granger12tunebeta3(xnum12,xnum22).num=0;
                        granger0(period0).granger21tunebeta3(xnum12,xnum22).num=0;
                        granger0(period0).granger12beta3(xnum12,xnum22).num=0;
                        granger0(period0).granger21beta3(xnum12,xnum22).num=0;
                        grangerdiffgamma3(xnum12,xnum22).num=0;
                    end
                end
            
                clear channel01 grangeralpha grangerbeta grangergamma cohalpha cohbeta cohgamma tuningalpha tuningbeta tuninggamma ;
            
                xnum1=0;gb0=0;ca0=0;cb0=0;cg0=0;ga0=0;gg0=0;ta0=0;tb0=0;tg0=0;
                for xnum2=1:numel(popg0)

                    channel1=popg0(xnum2).granger.channel1;
                    channel2=popg0(xnum2).granger.channel2;
                    
                    if isempty(popg0(xnum2).granger.channel1)==0 && isempty(popg0(xnum2).granger.channel2)==0
                        try
                            freq0=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21;
                        catch
                            continue
                        end
                        freq0=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.freq;
                        nnum1=find([channel0(:).num]==channel1,1);
                        nnum2=find([channel0(:).num]==channel2,1);
                        if isempty(nnum1)
                            continue
                        end
                        if isempty(nnum2)
                            continue
                        end
                        
                        minfidx = find(freq0 >= 5, 1);
                        maxfidx = find(freq0 <= 13, 1, 'last');
                        fpts_alpha = minfidx:maxfidx;
                        minfidx = find(freq0 >= beta_range(1), 1);
                        maxfidx = find(freq0 <= beta_range(2), 1, 'last');
                        fpts_beta = minfidx:maxfidx;
                        minfidx = find(freq0 >= 30, 1);
                        maxfidx = find(freq0 <= 55, 1, 'last');
                        fpts_gamma0 = minfidx:maxfidx;
                        minfidx = find(freq0 >= 65, 1);
                        maxfidx = find(freq0 <= 150, 1, 'last');
                        fpts_gamma = [fpts_gamma0 minfidx:maxfidx];
                        
                        freq00=freq0;
                        fpts_alpha0=fpts_alpha;
                        fpts_beta0=fpts_beta;
                        fpts_gamma00=fpts_gamma;
                                            
                        subplot(3,3,1);
                        connectivity0alpha2(nnum1,nnum2).num=connectivity0alpha2(nnum1,nnum2).num+1;
                        cnum=connectivity0alpha2(nnum1,nnum2).num;
                        connectivity0alpha2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_alpha));
                        plot(freq0(fpts_alpha),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_alpha));hold on;       
                        ca0=ca0+1;
                        cohalpha(ca0).freq=freq0(fpts_alpha);
                        cohalpha(ca0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_alpha);
                        title('coh alpha');
            
                        subplot(3,3,2);
                        connectivity0beta2(nnum1,nnum2).num=connectivity0beta2(nnum1,nnum2).num+1;
                        cnum=connectivity0beta2(nnum1,nnum2).num;
                        connectivity0beta2(nnum1,nnum2).data(cnum).val=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_beta));
                        connectivity0beta2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_beta);
                        plot(freq0(fpts_beta),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_beta));hold on;
                        title('coh beta');
            
%                        cohbeta                        
                        cb0=cb0+1;
                        cohbeta(cb0).freq=freq0(fpts_beta);
                        cohbeta(cb0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_beta);

                        subplot(3,3,3);
                        connectivity0gamma2(nnum1,nnum2).num=connectivity0gamma2(nnum1,nnum2).num+1;
                        cnum=connectivity0gamma2(nnum1,nnum2).num;
                        connectivity0gamma2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_gamma));
                        plot(freq0(fpts_gamma),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_gamma));hold on;                    
                        title('coh gamma');

                        %    cohgamma                        
                        cg0=cg0+1;
                        cohgamma(cg0).fpts=fpts_gamma;
                        cohgamma(cg0).freq=freq0(fpts_gamma);
                        cohgamma(cg0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_gamma);
            
                        freq0=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.freq;
                        minfidx = find(freq0 >= 5, 1);
                        maxfidx = find(freq0 <= 13, 1, 'last');
                        fpts_alpha = minfidx:maxfidx;
                        minfidx = find(freq0 >= beta_range(1), 1);
                        maxfidx = find(freq0 <= beta_range(2), 1, 'last');
                        fpts_beta = minfidx:maxfidx;
                        minfidx = find(freq0 >= 30, 1);
                        maxfidx = find(freq0 <= 50, 1, 'last');
                        fpts_gamma0 = minfidx:maxfidx;
                        minfidx = find(freq0 >= 80, 1);
                        maxfidx = find(freq0 <= 150, 1, 'last');
                        fpts_gamma = [fpts_gamma0 minfidx:maxfidx];
                        
                        freq1=freq0;
                        fpts_alpha1=fpts_alpha;
                        fpts_beta1=fpts_beta;
                        fpts_gamma1=fpts_gamma;
            
                        subplot(3,3,4);
                        granger12alpha2(nnum1,nnum2).num=granger12alpha2(nnum1,nnum2).num+1;
                        cnum=granger12alpha2(nnum1,nnum2).num;
                        granger12alpha2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha));
                        
                        granger21alpha2(nnum1,nnum2).num=granger21alpha2(nnum1,nnum2).num+1;
                        cnum=granger21alpha2(nnum1,nnum2).num;
        
                        granger12alpha2(nnum1,nnum2).data(cnum).val=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha));
                        granger12alpha2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha);
                        granger12alpha2(nnum1,nnum2).data(cnum).freq0=freq0(fpts_alpha);
                        granger21alpha2(nnum1,nnum2).data(cnum).val=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_alpha));
                        granger21alpha2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_alpha);
                        granger21alpha2(nnum1,nnum2).data(cnum).freq0=freq0(fpts_alpha);
                        
                        grangerdiffalpha2(nnum1,nnum2).num=grangerdiffalpha2(nnum1,nnum2).num+1;
                        cnum=grangerdiffalpha2(nnum1,nnum2).num;
                        grangerdiffalpha2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha) - ...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_alpha));            
                        plot(freq0(fpts_alpha),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha));hold on;   
                        plot(freq0(fpts_alpha),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_alpha));hold on;  
                        
                        granger12tunealpha(nnum1,nnum2).num=granger12tunealpha(nnum1,nnum2).num+1;
                        cnum=granger12tunealpha(nnum1,nnum2).num;
                        granger12tunealpha(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum12(fpts_alpha)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(1).granger.spectrum12(fpts_alpha);
                        granger21tunealpha(nnum1,nnum2).num=granger21tunealpha(nnum1,nnum2).num+1;
                        cnum=granger21tunealpha(nnum1,nnum2).num;
                        granger21tunealpha(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum21(fpts_alpha)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(1).granger.spectrum21(fpts_alpha);
                        
                        ga0=ga0+1;
                        grangeralpha(ga0).freq=freq0(fpts_alpha);
                        grangeralpha(ga0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_alpha);
                        ga0=ga0+1;
                        grangeralpha(ga0).freq=freq0(fpts_alpha);
                        grangeralpha(ga0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_alpha);

                        title('granger alpha');
            
                        subplot(3,3,5);
                        granger12beta2(nnum1,nnum2).num=granger12beta2(nnum1,nnum2).num+1;
                        granger21beta2(nnum1,nnum2).num=granger21beta2(nnum1,nnum2).num+1;
                        cnum=granger12beta2(nnum1,nnum2).num;
                        granger12beta2(nnum1,nnum2).data(cnum).val=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta));
                        granger12beta2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta);
                        granger12beta2(nnum1,nnum2).data(cnum).freq0=freq0(fpts_beta);
                        granger21beta2(nnum1,nnum2).data(cnum).val=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta));
                        granger21beta2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta);
                        granger21beta2(nnum1,nnum2).data(cnum).freq0=freq0(fpts_beta);
            
                        grangerdiffbeta2(nnum1,nnum2).num=grangerdiffbeta2(nnum1,nnum2).num+1;
                        cnum=grangerdiffbeta2(nnum1,nnum2).num;
                        grangerdiffbeta2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta));
                        grangerdiffbeta2(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta) - ...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta);
                        granger12tunebeta(nnum1,nnum2).num=granger12tunebeta(nnum1,nnum2).num+1;
                        cnum=granger12tunebeta(nnum1,nnum2).num;
                        granger12tunebeta(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum12(fpts_beta)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(1).granger.spectrum12(fpts_beta);
                        granger21tunebeta(nnum1,nnum2).num=granger21tunebeta(nnum1,nnum2).num+1;
                        cnum=granger21tunebeta(nnum1,nnum2).num;
                        granger21tunebeta(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum21(fpts_beta)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(1).granger.spectrum21(fpts_beta);

                        plot(freq0(fpts_beta),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta));hold on;   
                        plot(freq0(fpts_beta),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta));hold on;  
        
                        gb0=gb0+1;
                        grangerbeta(gb0).freq=freq0(fpts_beta);
                        grangerbeta(gb0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_beta);
                        gb0=gb0+1;
                        grangerbeta(gb0).freq=freq0(fpts_beta);
                        grangerbeta(gb0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_beta);

                        title('granger beta');
            
                        subplot(3,3,6);
                        granger12gamma2(nnum1,nnum2).num=granger12gamma2(nnum1,nnum2).num+1;
                        cnum=granger12gamma2(nnum1,nnum2).num;
                        granger12gamma2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_gamma));
                        granger21gamma2(nnum1,nnum2).num=granger21gamma2(nnum1,nnum2).num+1;
                        cnum=granger21gamma2(nnum1,nnum2).num;
                        granger21gamma2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_gamma));
                        grangerdiffgamma2(nnum1,nnum2).num=grangerdiffgamma2(nnum1,nnum2).num+1;
                        cnum=grangerdiffgamma2(nnum1,nnum2).num;
                        grangerdiffgamma2(nnum1,nnum2).val(cnum)=max(popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_gamma)-popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_gamma));
                                      
                        granger12tunegamma(nnum1,nnum2).num=granger12tunegamma(nnum1,nnum2).num+1;
                        cnum=granger12tunegamma(nnum1,nnum2).num;
                        granger12tunegamma(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum12(fpts_gamma)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_gamma);
                        granger21tunegamma(nnum1,nnum2).num=granger21tunegamma(nnum1,nnum2).num+1;
                        cnum=granger21tunegamma(nnum1,nnum2).num;
                        granger21tunegamma(nnum1,nnum2).data(cnum).spect=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(2).granger.spectrum21(fpts_gamma)-...
                            popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_gamma);
                        
                        plot(freq0(fpts_gamma),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_gamma));hold on;   
                        plot(freq0(fpts_gamma),popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_gamma));hold on;   
                        
                        gg0=gg0+1;
                        grangergamma(gg0).freq=fpts_gamma;
                        grangergamma(gg0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum12(fpts_gamma);
                        gg0=gg0+1;
                        grangergamma(gg0).freq=fpts_gamma;
                        grangergamma(gg0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).granger.spectrum21(fpts_gamma);
                        title('granger gamma');
            
                        subplot(3,3,7);
                        plot(freq0(fpts_alpha),granger12tunealpha(nnum1,nnum2).data(cnum).spect);hold on;   
                        plot(freq0(fpts_alpha),granger21tunealpha(nnum1,nnum2).data(cnum).spect);hold on;                           
                        ta0=ta0+1;
                        tuningalpha(ta0).freq=fpts_alpha;
                        tuningalpha(ta0).val=granger12tunealpha(nnum1,nnum2).data(cnum).spect;                       
                        ta0=ta0+1;
                        tuningalpha(ta0).freq=fpts_alpha;
                        tuningalpha(ta0).val=granger21tunealpha(nnum1,nnum2).data(cnum).spect;
            
                        subplot(3,3,8);
                        cnum=granger12tunegamma(nnum1,nnum2).num;
                        plot(freq0(fpts_beta),granger12tunebeta(nnum1,nnum2).data(cnum).spect);hold on;   
                        cnum=granger21tunegamma(nnum1,nnum2).num;
                        plot(freq0(fpts_beta),granger21tunebeta(nnum1,nnum2).data(cnum).spect);hold on; 
                        tb0=tb0+1;
                        tuningbeta(tb0).freq=fpts_beta;
                        tuningbeta(tb0).val=granger12tunebeta(nnum1,nnum2).data(cnum).spect;                       
                        tb0=tb0+1;
                        tuningbeta(tb0).freq=fpts_beta;
                        tuningbeta(tb0).val=granger21tunebeta(nnum1,nnum2).data(cnum).spect;

                        titname=['av-ap',cuename,periodname,apavname];
                        title(titname);
                        subplot(3,3,9);
                        plot(freq0(fpts_gamma),granger12tunegamma(nnum1,nnum2).data(cnum).spect);hold on;   
                        plot(freq0(fpts_gamma),granger21tunegamma(nnum1,nnum2).data(cnum).spect);hold on;   
                        tg0=tg0+1;
                        tuninggamma(tg0).freq=fpts_gamma;
                        tuninggamma(tg0).val=granger12tunegamma(nnum1,nnum2).data(cnum).spect;                       
                        tg0=tg0+1;
                        tuninggamma(tg0).freq=fpts_gamma;
                        tuninggamma(tg0).val=granger21tunegamma(nnum1,nnum2).data(cnum).spect;

                        xnum1=xnum1+1;
                        popg0(xnum2).name
            
                        channel01(xnum1).nnum1=nnum1;
                        channel01(xnum1).nnum2=nnum2;
                        channel01(xnum1).origin_channel=channel1;
                        channel01(xnum1).target_channel=channel2;
                        
                        for xnum5=1:numel(fsum.freq_channel)
                            if isempty(find(fsum.freq_channel(xnum5).channel==channel1,1))~=1
                                channel01(xnum1).origin_channel=channel1;
                                channel01(xnum1).origin_group=xnum5;
                                channel01(xnum1).origin_name=[fsum.freq_channel(xnum5).group_name,':',num2str(channel0(xnum0).num)];
                                channel01(xnum1).origin_group_name=fsum.freq_channel(xnum5).group_name;
                            end
                            if isempty(find(fsum.freq_channel(xnum5).channel==channel2,1))~=1
                                channel01(xnum1).target_channel=channel2;
                                channel01(xnum1).target_group=xnum5;
                                channel01(xnum1).target_name=[fsum.freq_channel(xnum5).group_name,':',num2str(channel0(xnum0).num)];
                                channel01(xnum1).target_group_name=fsum.freq_channel(xnum5).group_name;
                            end
                        end
                        origin_group_num=channel01(xnum1).origin_group;
                        target_group_num=channel01(xnum1).target_group;
            
                        % coh alpha
                        connectivity0alpha3(origin_group_num,target_group_num).num=connectivity0alpha3(origin_group_num,target_group_num).num+1;
                        cnum=connectivity0alpha3(origin_group_num,target_group_num).num;
                        cnum0=connectivity0alpha2(nnum1,nnum2).num;
                        connectivity0alpha3(origin_group_num,target_group_num).val(cnum)=connectivity0alpha2(nnum1,nnum2).val(cnum0);
                        connectivity0alpha3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % coh beta
                        connectivity0beta3(origin_group_num,target_group_num).num=connectivity0beta3(origin_group_num,target_group_num).num+1;
                        cnum=connectivity0beta3(origin_group_num,target_group_num).num;
                        cnum0=connectivity0beta2(nnum1,nnum2).num;
                        connectivity0beta3(origin_group_num,target_group_num).val(cnum)=connectivity0beta2(nnum1,nnum2).data(cnum0).val;
                        connectivity0beta3(origin_group_num,target_group_num).data(cnum).spect=connectivity0beta2(nnum1,nnum2).data(cnum0).spect;
                        connectivity0beta3(origin_group_num,target_group_num).channel1(cnum)=channel1;
                        connectivity0beta3(origin_group_num,target_group_num).channel2(cnum)=channel2;
                        connectivity0beta3(origin_group_num,target_group_num).origin_group_name=channel01(xnum1).origin_group_name;
                        connectivity0beta3(origin_group_num,target_group_num).target_group_name=channel01(xnum1).target_group_name;
            %            connectivity0beta3(origin_group_num,target_group_num).time=time0;
                        % coh gamma
                        connectivity0gamma3(origin_group_num,target_group_num).num=connectivity0gamma3(origin_group_num,target_group_num).num+1;
                        cnum=connectivity0gamma3(origin_group_num,target_group_num).num;
                        cnum0=connectivity0gamma2(nnum1,nnum2).num;
                        connectivity0gamma3(origin_group_num,target_group_num).val(cnum)=connectivity0gamma2(nnum1,nnum2).val(cnum0);
                        connectivity0gamma3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        
                        % granger12 alpha
                        granger12alpha3(origin_group_num,target_group_num).num=granger12alpha3(origin_group_num,target_group_num).num+1;
                        cnum=granger12alpha3(origin_group_num,target_group_num).num;
                        cnum0=granger12alpha2(nnum1,nnum2).num;
                        granger12alpha3(origin_group_num,target_group_num).val(cnum)=granger12alpha2(nnum1,nnum2).val(cnum0);
                        granger12alpha3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % granger12 beta  %granger12tunebeta(nnum1,nnum2).num
                        granger12beta3(origin_group_num,target_group_num).num=granger12beta3(origin_group_num,target_group_num).num+1;
                        cnum=granger12beta3(origin_group_num,target_group_num).num;
                        cnum0=granger12beta2(nnum1,nnum2).num;
                        granger12beta3(origin_group_num,target_group_num).data(cnum).val=granger12beta2(nnum1,nnum2).data(cnum0).val;
                        granger0(period0).granger12beta3(origin_group_num,target_group_num).num=granger12beta3(origin_group_num,target_group_num).num;
                        granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect=granger12beta2(nnum1,nnum2).data(cnum0).spect;
                        granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0=granger12beta2(nnum1,nnum2).data(cnum0).freq0;
                        granger12beta3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % granger12 gamma
                        granger12gamma3(origin_group_num,target_group_num).num=granger12gamma3(origin_group_num,target_group_num).num+1;
                        cnum=granger12gamma3(origin_group_num,target_group_num).num;
                        cnum0=granger12gamma2(nnum1,nnum2).num;
                        granger12gamma3(origin_group_num,target_group_num).val(cnum)=granger12gamma2(nnum1,nnum2).val(cnum0);
                        granger12gamma3(origin_group_num,target_group_num).channel(cnum)=channel1;
            
                        % granger21 alpha
                        granger21alpha3(origin_group_num,target_group_num).num=granger21alpha3(origin_group_num,target_group_num).num+1;
                        cnum=granger21alpha3(origin_group_num,target_group_num).num;
                        cnum0=granger21alpha2(nnum1,nnum2).num;
                        granger21alpha3(origin_group_num,target_group_num).data(cnum).val=granger21alpha2(nnum1,nnum2).data(cnum0).val;
                        granger21alpha3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % granger21 beta
                        granger21beta3(origin_group_num,target_group_num).num=granger21beta3(origin_group_num,target_group_num).num+1;
                        cnum=granger21beta3(origin_group_num,target_group_num).num;
                        cnum0=granger21beta2(nnum1,nnum2).num;
            %            granger21beta3(origin_group_num,target_group_num).data(cnum).val=granger21beta2(nnum1,nnum2).data(cnum0).val;
                        granger0(period0).granger21beta3(origin_group_num,target_group_num).num=granger21beta3(origin_group_num,target_group_num).num;
                        granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect=granger21beta2(nnum1,nnum2).data(cnum0).spect;
                        granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).freq0=granger21beta2(nnum1,nnum2).data(cnum0).freq0;
                        granger21beta3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % granger21 gamma
                        granger21gamma3(origin_group_num,target_group_num).num=granger21gamma3(origin_group_num,target_group_num).num+1;
                        cnum=granger21gamma3(origin_group_num,target_group_num).num;
                        cnum0=granger21gamma2(nnum1,nnum2).num;
                        granger21gamma3(origin_group_num,target_group_num).val(cnum)=granger21gamma2(nnum1,nnum2).val(cnum0);
                        granger21gamma3(origin_group_num,target_group_num).channel(cnum)=channel1;
            
                        % grangerdiff alpha
                        grangerdiffalpha3(origin_group_num,target_group_num).num=grangerdiffalpha3(origin_group_num,target_group_num).num+1;
                        cnum=grangerdiffalpha3(origin_group_num,target_group_num).num;
                        cnum0=grangerdiffalpha2(nnum1,nnum2).num;
                        grangerdiffalpha3(origin_group_num,target_group_num).val(cnum)=grangerdiffalpha2(nnum1,nnum2).val(cnum0);
                        grangerdiffalpha3(origin_group_num,target_group_num).channel(cnum)=channel1;
                        % grangerdiff beta
                        grangerdiffbeta3(origin_group_num,target_group_num).num=grangerdiffbeta3(origin_group_num,target_group_num).num+1;
                        cnum=grangerdiffbeta3(origin_group_num,target_group_num).num;
                        cnum0=grangerdiffbeta2(nnum1,nnum2).num;
                        grangerdiffbeta3(origin_group_num,target_group_num).val(cnum)=grangerdiffbeta2(nnum1,nnum2).val(cnum0);
                        grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).spect=grangerdiffbeta2(nnum1,nnum2).data(cnum0).spect;
                        grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).freq0=granger12beta2(nnum1,nnum2).data(cnum0).freq0;
                        granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).num=grangerdiffbeta3(origin_group_num,target_group_num).num;
                        granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).spect=grangerdiffbeta2(nnum1,nnum2).data(cnum0).spect;
                        granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).freq0=granger12beta2(nnum1,nnum2).data(cnum0).freq0;
                        grangerdiffbeta3(origin_group_num,target_group_num).channel(cnum)=channel1;

                        % granger tune beta %granger12tunebeta(nnum1,nnum2).num
                        granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).num=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).num+1;
                        cnum=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).num;
                        cnum0=granger12tunebeta(nnum1,nnum2).num;%granger0(period0).granger12tunebeta3(
                        granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).spect=granger12tunebeta(nnum1,nnum2).data(cnum0).spect;
                        granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).freq0=granger12beta2(nnum1,nnum2).data(cnum0).freq0;
                        granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).num=granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).num+1;
                        cnum=granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).num;
                        cnum0=granger21tunebeta(nnum1,nnum2).num;
                        granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum).spect=granger21tunebeta(nnum1,nnum2).data(cnum0).spect;
                        granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum).freq0=granger21beta2(nnum1,nnum2).data(cnum0).freq0;

                        % grangerdiff gamma
                        grangerdiffgamma3(origin_group_num,target_group_num).num=grangerdiffgamma3(origin_group_num,target_group_num).num+1;
                        cnum=grangerdiffgamma3(origin_group_num,target_group_num).num;
                        cnum0=grangerdiffgamma2(nnum1,nnum2).num;
                        grangerdiffgamma3(origin_group_num,target_group_num).val(cnum)=grangerdiffgamma2(nnum1,nnum2).val(cnum0);
                        grangerdiffgamma3(origin_group_num,target_group_num).channel(cnum)=channel1;
                    end
                end
                % connectivity beta
                nnum1
                nnum2
                nnum00=0;
                              
                clear channel1val myLabel1 matrix0 matrix3 matrix11 matrix4 matrix2 matrix_sig matrix00 matrix_asym;
                onum0=0;
                for origin_group_num=1:numel(fsum.freq_channel)
                    for target_group_num=1:numel(fsum.freq_channel)
                        totnum=connectivity0beta3(origin_group_num,target_group_num).num;
                        clear val;
                        for tnum0=1:totnum
                            val(tnum0)=connectivity0beta3(origin_group_num,target_group_num).val(tnum0);
                        end
                        if totnum>0
                            matrix0{origin_group_num,target_group_num}=val;
                        else
                            matrix0{origin_group_num,target_group_num}=nan;
                        end
                        matrix00(origin_group_num,target_group_num)=0;
                        matrix_asym(origin_group_num,target_group_num)=0;
                        matrix2{origin_group_num,target_group_num}=[];
                        matrix_sig{origin_group_num,target_group_num}=[];
            %             matrix3{origin_group_num,target_group_num}=nan;
            %             matrix11{origin_group_num,target_group_num}=nan;
            %             matrix4{origin_group_num,target_group_num}=nan;
                    end
                end
                matrix0
                clear matrix1
                clear myLabel myColorMap
                for origin_group_num=1:numel(fsum.freq_channel)
                    for target_group_num=1:numel(fsum.freq_channel)
                        val=matrix0{origin_group_num,target_group_num};
                        if monkey0==3
                            if target_group_num==1 && (origin_group_num==4 || origin_group_num==5)
                                val=0;
                            end
                        end
                        if ~isnan(val)
                            matrix1(origin_group_num,target_group_num)=nan;
                        else
                            matrix1(origin_group_num,target_group_num)=mean(val);
                        end
                    end
                    myLabel{origin_group_num} = fsum.freq_channel(origin_group_num).group_name;
                    switch origin_group_num
                        case 1
                          myColorMap(origin_group_num,:) =[1 0 0];
                        case 2
                          myColorMap(origin_group_num,:) =[0 1 0];          
                        case 3
                          myColorMap(origin_group_num,:) =[0 0 1];  
                        case 4
                          myColorMap(origin_group_num,:) =[0.5 0.8 1];  
                        case 5
                          myColorMap(origin_group_num,:) =[1 0 1];  
                        case 6
                          myColorMap(origin_group_num,:) =[0 0.5 0.5];  
                        case 7
                          myColorMap(origin_group_num,:) =[0.8 0.8 0.1];  
                        otherwise
                          myColorMap(origin_group_num,:) =[0 0 0];  
                    end     
                end
                matrix1
                myLabel{4} = [myLabel{4},'coherence']; 
                try
                    myLabel{5} = [myLabel{5},[cuename,'-',periodname,'-',apavname]];
                catch
                    myLabel{2} = [myLabel{2},[cuename,'-',periodname,'-',apavname]];
                end
        
                % [figure_path0,foldername0]
                % figname=['Coh_net1',cuename,periodname,apavname];
                % print_file10([figure_path0,foldername0],figname,1,1);
                close;
        %        exportgraphics(gcf,['Figure-',cuename,'-',periodname,'2.eps'],'ContentType','vector');
                
                figure;
                clear val22 mean0 sem0;
                numel1=numel(grangeralpha(1).freq);
                for gb1=1:ga0
                    for xnum10=1:grangeralpha(gb1).freq(numel1)-grangeralpha(gb1).freq(1)
                        val22(xnum10,gb1)=grangeralpha(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:grangeralpha(gb1).freq(numel1)-grangeralpha(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(ga0);
                end
                subplot(3,3,4);
                plot(grangeralpha(gb1).freq(1):grangeralpha(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(grangeralpha(gb1).freq(1):grangeralpha(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(grangeralpha(gb1).freq(1):grangeralpha(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('granger alpha');

                clear val22 mean0 sem0;
                numel1=numel(grangerbeta(1).freq);
                for gb1=1:gb0
                    for xnum10=1:grangerbeta(gb1).freq(numel1)-grangerbeta(gb1).freq(1)
                        val22(xnum10,gb1)=grangerbeta(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:grangerbeta(gb1).freq(numel1)-grangerbeta(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(gb0);
                end
                subplot(3,3,5);
                plot(grangerbeta(gb1).freq(1):grangerbeta(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(grangerbeta(gb1).freq(1):grangerbeta(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(grangerbeta(gb1).freq(1):grangerbeta(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('granger beta');
%                         cb0=cb0+1;
%                         cohbeta(cb0).freq=freq0(fpts_beta);
%                         cohbeta(cb0).val=popg0(xnum2).granger.stim01stim12resid3all4(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).coh.spectrum(fpts_beta);
                clear val22 mean0 sem0;
                numel1=numel(cohbeta(1).freq);
                for gb1=1:cb0
                    for xnum10=1:cohbeta(gb1).freq(numel1)-cohbeta(gb1).freq(1)
                        val22(xnum10,gb1)=cohbeta(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:cohbeta(gb1).freq(numel1)-cohbeta(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(cb0);
                end
                subplot(3,3,2);
                plot(cohbeta(gb1).freq(1):cohbeta(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(cohbeta(gb1).freq(1):cohbeta(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(cohbeta(gb1).freq(1):cohbeta(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('coh beta');

                clear val22 mean0 sem0;
                numel1=numel(cohalpha(1).freq);
                for gb1=1:ca0
                    for xnum10=1:cohalpha(gb1).freq(numel1)-cohalpha(gb1).freq(1)
                        val22(xnum10,gb1)=cohalpha(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:cohalpha(gb1).freq(numel1)-cohalpha(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(ca0);
                end
                subplot(3,3,1);
                plot(cohalpha(gb1).freq(1):cohalpha(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(cohalpha(gb1).freq(1):cohalpha(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(cohalpha(gb1).freq(1):cohalpha(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('coh alpha');

                clear val22;
                for gb1=1:cg0
                    numel1(gb1)=numel(cohgamma(gb1).freq);
                    for xnum10=1:numel1(gb1)
                        val22(xnum10,gb1)=cohgamma(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:numel1(gb1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(cg0);
                end
                subplot(3,3,3);
                plot(cohgamma(gb1).freq(1):cohgamma(gb1).freq(26),mean0(1:26)','-r');hold on;
                plot(cohgamma(gb1).freq(27):cohgamma(gb1).freq(numel1),mean0(27:numel1)','-r');hold on;
                plot(cohgamma(gb1).freq(1):cohgamma(gb1).freq(26),mean0(1:26)'+sem0(1:26)','-r');hold on;
                plot(cohgamma(gb1).freq(27):cohgamma(gb1).freq(numel1),mean0(27:numel1)'+sem0(27:numel1)','-r');hold on;
                plot(cohgamma(gb1).freq(1):cohgamma(gb1).freq(26),mean0(1:26)'-sem0(1:26)','-r');hold on;
                plot(cohgamma(gb1).freq(27):cohgamma(gb1).freq(numel1),mean0(27:numel1)'-sem0(27:numel1)','-r');hold on;
                title('coh gamma');
                                
                clear val22;
                for gb1=1:gg0
                    numel1(gb1)=numel(grangergamma(gb1).freq);
                    for xnum10=1:numel1(gb1)
                        val22(xnum10,gb1)=grangergamma(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:numel1(gb1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(gg0);
                end
                subplot(3,3,6);
                plot(grangergamma(gb1).freq(1):grangergamma(gb1).freq(21),mean0(1:21)','-r');hold on;
                plot(grangergamma(gb1).freq(22):grangergamma(gb1).freq(numel1),mean0(22:numel1)','-r');hold on;
                plot(grangergamma(gb1).freq(1):grangergamma(gb1).freq(21),mean0(1:21)'+sem0(1:21)','-r');hold on;
                plot(grangergamma(gb1).freq(22):grangergamma(gb1).freq(numel1),mean0(22:numel1)'+sem0(22:numel1)','-r');hold on;
                plot(grangergamma(gb1).freq(1):grangergamma(gb1).freq(21),mean0(1:21)'-sem0(1:21)','-r');hold on;
                plot(grangergamma(gb1).freq(22):grangergamma(gb1).freq(numel1),mean0(22:numel1)'-sem0(22:numel1)','-r');hold on;
                title('granger gamma');

                clear val22 mean0 sem0;
                numel1=numel(tuningalpha(1).freq);
                for gb1=1:ta0
                    for xnum10=1:tuningalpha(gb1).freq(numel1)-tuningalpha(gb1).freq(1)
                        val22(xnum10,gb1)=tuningalpha(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:tuningalpha(gb1).freq(numel1)-tuningalpha(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(ta0);
                end
                subplot(3,3,7);
                plot(tuningalpha(gb1).freq(1):tuningalpha(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(tuningalpha(gb1).freq(1):tuningalpha(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(tuningalpha(gb1).freq(1):tuningalpha(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('av-ap alpha');

                clear val22 mean0 sem0;
                numel1=numel(tuningbeta(1).freq);
                for gb1=1:tb0
                    for xnum10=1:tuningbeta(gb1).freq(numel1)-tuningbeta(gb1).freq(1)
                        val22(xnum10,gb1)=tuningbeta(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:tuningbeta(gb1).freq(numel1)-tuningbeta(gb1).freq(1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(tb0);
                end
                subplot(3,3,8);
                plot(tuningbeta(gb1).freq(1):tuningbeta(gb1).freq(numel1)-1,mean0','-r');hold on;
                plot(tuningbeta(gb1).freq(1):tuningbeta(gb1).freq(numel1)-1,mean0'+sem0','-r');hold on;
                plot(tuningbeta(gb1).freq(1):tuningbeta(gb1).freq(numel1)-1,mean0'-sem0','-r');hold on;
                title('av-ap beta');

                clear val22;
                for gb1=1:tg0
                    numel1(gb1)=numel(tuninggamma(gb1).freq);
                    for xnum10=1:numel1(gb1)
                        val22(xnum10,gb1)=tuninggamma(gb1).val(xnum10);
                    end                            
                end                
                for x0=1:numel1(gb1)
                   mean0(x0)=mean(val22(x0,:));                   
                   sem0(x0)=std(val22(x0,:))/sqrt(tg0);
                end
                subplot(3,3,9);
                plot(tuninggamma(gb1).freq(1):tuninggamma(gb1).freq(21),mean0(1:21)','-r');hold on;
                plot(tuninggamma(gb1).freq(22):tuninggamma(gb1).freq(numel1),mean0(22:numel1)','-r');hold on;
                plot(tuninggamma(gb1).freq(1):tuninggamma(gb1).freq(21),mean0(1:21)'+sem0(1:21)','-r');hold on;
                plot(tuninggamma(gb1).freq(22):tuninggamma(gb1).freq(numel1),mean0(22:numel1)'+sem0(22:numel1)','-r');hold on;
                plot(tuninggamma(gb1).freq(1):tuninggamma(gb1).freq(21),mean0(1:21)'-sem0(1:21)','-r');hold on;
                plot(tuninggamma(gb1).freq(22):tuninggamma(gb1).freq(numel1),mean0(22:numel1)'-sem0(22:numel1)','-r');hold on;
                title(['av-ap gamma',' ',num2str(numel(tuninggamma))]);

                [figure_path0,foldername0]
                figname=['Coh_net2',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);        
           
                % granger beta
                nnum1
                nnum2
                nnum00=0;
                clear channel1val myLabel1;
                onum0=0;

                if addsakura==1
                    for origin_group_num=1:7%numel(fsum.freq_channel)
                        for target_group_num=1:7%numel(fsum.freq_channel)
                            if isempty(granger0_sakura(period0).grangerdiffbeta3)
                                continue
                            end
                            last00=numel(granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).data);
                            for cnum=last00+1:last00+numel(granger0_sakura(period0).grangerdiffbeta3(origin_group_num,target_group_num).data)
                                if isempty(granger0_sakura(period0).grangerdiffbeta3(origin_group_num,target_group_num).data)==1
                                    continue
                                end
                                granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).spect=...
                                     granger0_sakura(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum-last00).spect;
                                granger0(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum).freq0=...
                                     granger0_sakura(period0).grangerdiffbeta3(origin_group_num,target_group_num).data(cnum-last00).freq0;
                            end
                        end
                    end
                    for origin_group_num=1:7%numel(fsum.freq_channel)
                        for target_group_num=1:7%numel(fsum.freq_channel)
                            if isempty(granger0_sakura(period0).granger12beta3)
                                continue
                            end
                            last00=numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data);
                            for cnum=last00+1:last00+numel(granger0_sakura(period0).granger12beta3(origin_group_num,target_group_num).data)
                                if isempty(granger0_sakura(period0).granger12beta3(origin_group_num,target_group_num).data)==1
                                    continue
                                end
                                granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect=...
                                     granger0_sakura(period0).granger12beta3(origin_group_num,target_group_num).data(cnum-last00).spect;
                                granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0=...
                                     granger0_sakura(period0).granger12beta3(origin_group_num,target_group_num).data(cnum-last00).freq0;
                                granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect=...
                                     granger0_sakura(period0).granger21beta3(origin_group_num,target_group_num).data(cnum-last00).spect;
                                granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).freq0=...
                                     granger0_sakura(period0).granger21beta3(origin_group_num,target_group_num).data(cnum-last00).freq0;
                            end
                        end
                    end
                    for origin_group_num=1:7%numel(fsum.freq_channel)
                        for target_group_num=1:7%numel(fsum.freq_channel)
                            if isempty(granger0_sakura(period0).granger12tunebeta3)
                                continue
                            end
                            last00=numel(granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data);
                            for cnum=last00+1:last00+numel(granger0_sakura(period0).granger12tunebeta3(origin_group_num,target_group_num).data)
                                if isempty(granger0_sakura(period0).granger12tunebeta3(origin_group_num,target_group_num).data)==1
                                    continue
                                end
                                granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).spect=...
                                     granger0_sakura(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum-last00).spect;
                                granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).freq0=...
                                     granger0_sakura(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum-last00).freq0;
                                granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum).spect=...
                                     granger0_sakura(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum-last00).spect;
                                granger0(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum).freq0=...
                                     granger0_sakura(period0).granger21tunebeta3(origin_group_num,target_group_num).data(cnum-last00).freq0;
                            end
                        end
                    end
                end
%%%
%%%
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];
                        if monkey0==3
                            if origin_group_num==3 && target_group_num==1 
                                for cnum=1:numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data)
                                    granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect;%./2;
                                end
                                for cnum=1:numel(granger0(period0).granger21beta3(target_group_num,origin_group_num).data)
                                    granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect;%./2;
                                end
                            end
                        end
                        for cnum=1:numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data)
                            val0(cnum).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect; 
                            freq0=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum).val);hold on;
                        end
                        for cnum=1:numel(granger0(period0).granger21beta3(target_group_num,origin_group_num).data)
                            val0(cnum).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect; 
                            freq0=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum).val);hold on;
                        end
                        title([cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(numel(val0))]);                      
                        grangerout(period0).region(origin_group_num,target_group_num).val0=val0;
                        grangerout(period0).region(origin_group_num,target_group_num).freq0=freq0;
                    end
                end
                [figure_path0,foldername0]
                figname=['Granger_all2',cuename,periodname];
                print_file10([figure_path0,foldername0],figname,1,1);
%                close;
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];
                        for cnum=1:numel(granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data)
                            val0(cnum).val=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).spect; 
                            freq0=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum).val);hold on;
                        end
                        for cnum=1:numel(granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data)
                            val0(cnum).val=granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data(cnum).spect; 
                            freq0=granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum).val);hold on;
                        end
                        title([cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(numel(val0))]);
                    end
                end
                [figure_path0,foldername0]
                figname=['Granger_tune1',cuename,periodname];
                print_file10([figure_path0,foldername0],figname,1,1);
                %close;
        %        exportgraphics(gcf,['Figure-',cuename,'-',periodname,'4.eps'],'ContentType','vector');
                % assymmetry index
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];val0o=[];val0t=[];cnum1=0;val1=[];
                        for cnum=1:numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data)
                            cnum1=cnum1+1;
                            val0(cnum1).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect-granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect; 
                            val0(cnum1).valsum=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect+granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect; 
                            %val0(cnum1).valsum=granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect; 
                            val0o(cnum1).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect;
                            val0t(cnum1).val=granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect;
                            freq0=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum1).val./val0(cnum1).valsum);hold on;
                            val1(cnum1).val=val0(cnum1).val./val0(cnum1).valsum;
%                            plot(freq0,val0(cnum1).val);hold on;
                        end
                        for cnum=1:numel(granger0(period0).granger21beta3(target_group_num,origin_group_num).data)
                            cnum1=cnum1+1;
                            val0(cnum1).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect-granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect; 
                            val0(cnum1).valsum=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect+granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect; 
                            %val0(cnum1).valsum=granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect; 
                            val0o(cnum1).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect;
                            val0t(cnum1).val=granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect;
                            freq0=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).freq0;
                            plot(freq0,val0(cnum1).val./val0(cnum1).valsum);hold on;
                            val1(cnum1).val=val0(cnum1).val./val0(cnum1).valsum;
                        end
                        title(['asym:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(numel(val0))]);
                        asymout(period0).region(origin_group_num,target_group_num).val1=val1;
                        asymout(period0).region(origin_group_num,target_group_num).freq0=freq0;
                    end
                end
%                close;
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                clear val122
                for origin_group_num=1:numel(fsum.freq_channel)
                    for target_group_num=1:numel(fsum.freq_channel)
                        matrix_sig{origin_group_num,target_group_num}=[];
                        val122(origin_group_num,target_group_num).val2=[];
                    end
                end
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];
                        val1=[];
                        cnum1=0;
                        for cnum=1:numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data)
                            val0(cnum).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val1(cnum1).val=val0(cnum).val;
                        end
                        for cnum=1:numel(granger0(period0).granger21beta3(target_group_num,origin_group_num).data)
                            val0(cnum).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val1(cnum1).val=val0(cnum).val;
                        end
%                         try
%                             val11=[val11,val1]
%                         end
                        clear mean1 mean0 sem0 
                        val2=[];mean1=[];
                        h00=0;
                        minnum0=nan;
                        if isempty(val1)==0
                            for vnum1=1:length(val1)
                                val2(:,vnum1)=val1(vnum1).val;
                            end
                            [x1,y1]=size(val2);
                            for x0=1:x1
                                mean0(x0)=mean(val2(x0,:));
                                sem0(x0)=std(val2(x0,:))/sqrt(y1);
                            end
                            for y0=1:y1
                                mean1(y0)=mean(val2(:,y0));
%                                plot(20,mean1(y0),'or','markersize',4);hold on;                   
                            end
                            val122(origin_group_num,target_group_num).val2=val2;%x: freq, y: channel
                            plot(freq0,mean0','-r');hold on;
                            plot(freq0,mean0'-sem0','-r');hold on;
                            plot(freq0,mean0'+sem0','-r');hold on;
                            plot([freq0(1) freq0(numel(freq0))],[0 0],'-k');
                            plot(18,mean(mean1),'+r','markersize',10);
                            plot(18,mean(mean1)+std(mean1),'+r','markersize',10);
                            plot(18,mean(mean1)-std(mean1),'+r','markersize',10);
                            h00=0;h01=0;
                            try
                                [h00,p00]=ttest(mean1,0);
                            catch
                                h00=0;p00=1;
                            end
                            matrix_sig{target_group_num,origin_group_num}=mean1;
                            h01=0;p01=1;
                            stats.tstat=0;
                            if isempty(matrix_sig{origin_group_num,target_group_num})~=1
                                mean2=matrix_sig{origin_group_num,target_group_num};
                                [h01,p01,ci,stats]=ttest(mean1,mean2);
%                                stats;
                            end
                            clear h012 p012 stats2 val2val val1222val;
                            alpha0=0.01;
                            % if target_group_num==2 && origin_group_num==3
                            %     alpha0=0.05;
                            % end
                            if isempty(val122(target_group_num,origin_group_num).val2)~=1 && target_group_num~=origin_group_num
                                val1222=val122(target_group_num,origin_group_num).val2;
                                [x111,y111]=size(val2);
                                for xnum1222=1:x111
                                    [h012(xnum1222),p012(xnum1222),ci,stats2(xnum1222).stats]=ttest(val2(xnum1222,:),val1222(xnum1222,:),'Alpha',alpha0/x111);
                                    val2val(xnum1222)=mean(val2(xnum1222,:));
                                    val1222val(xnum1222)=mean(val1222(xnum1222,:));
                                end
                            else
                                try
                                    x111;
                                catch
                                    continue
                                end
                                for xnum1222=1:x111
                                    h012(xnum1222)=0;
                                    val2val(xnum1222)=nan;
                                    val1222val(xnum1222)=nan;
                                end
                            end
                            matrix_sig1(origin_group_num,target_group_num)=0;
                            matrix0{origin_group_num,target_group_num}=0;
                            matrix00(origin_group_num,target_group_num)=0;
                            if isnan(h01)==0
                                if h01
                                    text(25,mean(mean1)+std(mean1),['p=',num2str(p01)]);
                                    matrix_sig1(origin_group_num,target_group_num)=1;
                                    if stats.tstat>0
                                        matrix0{target_group_num,origin_group_num}=stats.tstat;%1000*(max(mean1-mean2));
                                    else
                                        matrix0{origin_group_num,target_group_num}=-stats.tstat;%1000*(max(mean2-mean1));
                                    end
                                    matrix00(origin_group_num,target_group_num)=stats.tstat;
                                    matrix00(target_group_num,origin_group_num)=-stats.tstat;
                                end
                            end
                            title([cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(y1)]);
            
                            subplot(4,4,origin_group_num+4*(target_group_num-1));
                            clear diff000
                            for y0=1:numel(h012)
                                if h012(y0)
                                    line([freq0(y0) freq0(y0)],[val2val(y0) val1222val(y0)],'Color',[1,0.8,0.8],'LineWidth',8);hold on;
%                                    text(freq0(y0),max([val2val(y0) val1222val(y0)]),num2str(abs(val2val(y0)-val1222val(y0)),4));
                                    diff000(y0)=(abs(val2val(y0)-val1222val(y0)));
                                else
                                    diff000(y0)=nan;
                                end
                            end
                            [x12,y12]=max(diff000);
                            text(freq0(y12),max([val2val(y12) val1222val(y12)]),num2str(x12,4));
                            for y0=1:y1
                                mean1(y0)=mean(val2(:,y0));
                            end
                            plot(freq0,mean0','-b');hold on;
                            plot(freq0,mean0'-sem0','-b');hold on;
                            plot(freq0,mean0'+sem0','-b');hold on;
                            plot(25,mean(mean1),'+b','markersize',10);
                            plot(25,mean(mean1)+std(mean1),'+b','markersize',10);
                            plot(25,mean(mean1)-std(mean1),'+b','markersize',10);

                            subplot(4,4,target_group_num+4*(origin_group_num-1));
                            clear diff000
                            diff000=[];x12=nan;
                            for y0=1:numel(h012)
                                if h012(y0)
                                    line([freq0(y0) freq0(y0)],[val2val(y0) val1222val(y0)],'Color',[1,0.8,0.8],'LineWidth',8);hold on;
%                                    text(freq0(y0),max([val2val(y0) val1222val(y0)]),num2str(abs(val2val(y0)-val1222val(y0)),4));
                                    diff000(y0)=(abs(val2val(y0)-val1222val(y0)));
                                else
                                    diff000(y0)=nan;
                                end
                            end
                            [x12,y12]=max(diff000);
                            text(freq0(y12),max([val2val(y12) val1222val(y12)]),num2str(x12,4));
                            minnum0=min([val2val(y12) val1222val(y12)]);
                            for y0=1:y1
                                mean1(y0)=mean(val2(:,y0));
%                                plot(23,mean1(y0),'ob','markersize',4);hold on;                   
                            end
                            plot(freq0,mean0','-r');hold on;
                            plot(freq0,mean0'-sem0','-r');hold on;
                            plot(freq0,mean0'+sem0','-r');hold on;
                            plot([freq0(1) freq0(numel(freq0))],[0 0],'-k');

                            try
                                maxy(origin_group_num,target_group_num)
                            catch
                                maxy(origin_group_num,target_group_num)=0;
                            end
                            if maxy(origin_group_num,target_group_num)==0
                                maxy(origin_group_num,target_group_num)=max(mean0+sem0);
                            else
                                if maxy(origin_group_num,target_group_num)<=max(mean0+sem0)
                                    maxy(origin_group_num,target_group_num)=max(mean0+sem0);
                                end
                            end
                            axis([freq0(1) freq0(numel(freq0)) 0 3.1*max(mean0+sem0)]);%maxy(origin_group_num,target_group_num)]);
                            if period0==2
                                if (origin_group_num==4 && target_group_num==3) || (origin_group_num==3 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.025]);
                                elseif (origin_group_num==4 && target_group_num==2) || (origin_group_num==2 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.020]);
                                elseif (origin_group_num==4 && target_group_num==1) || (origin_group_num==1 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.05]);
                                elseif (origin_group_num==3 && target_group_num==1) || (origin_group_num==1 && target_group_num==3)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.03]);
                                elseif (origin_group_num==3 && target_group_num==2) || (origin_group_num==2 && target_group_num==3)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.15]);
                                elseif (origin_group_num==2 && target_group_num==1) || (origin_group_num==1 && target_group_num==2)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.045]);
                                end
                            else
                                if (origin_group_num==4 && target_group_num==3) || (origin_group_num==3 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.017]);
                                elseif (origin_group_num==4 && target_group_num==2) || (origin_group_num==2 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.007]);
                                elseif (origin_group_num==4 && target_group_num==1) || (origin_group_num==1 && target_group_num==4)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.009]);
                                elseif (origin_group_num==3 && target_group_num==1) || (origin_group_num==1 && target_group_num==3)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.006]);
                                elseif (origin_group_num==3 && target_group_num==2) || (origin_group_num==2 && target_group_num==3)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.02]);
                                elseif (origin_group_num==2 && target_group_num==1) || (origin_group_num==1 && target_group_num==2)
                                    axis([freq0(1) freq0(numel(freq0)) 0 0.008]);
                                end
                            end
                        end
                        val3=mean(mean1);
                        try 
                            x12;
                        catch
                            x12=nan;
                        end
%                        matrix0{origin_group_num,target_group_num}=mean1;
                        matrix2{origin_group_num,target_group_num}=mean1;
                        matrix0{origin_group_num,target_group_num}=mean1;%x12/minnum0;
                        periodmat0(period0).matrix_final(origin_group_num,target_group_num)=x12/minnum0;
                        periodmat0(period0).matrix_final1(origin_group_num,target_group_num)=x12;


                        grangerout1(period0).val122=val122;
                    end
                end
                [figure_path0,foldername0]
                figname=['Granger_direction12',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
%                close;
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val4=(val122(origin_group_num,target_group_num).val2)';
                        granger11(period0).direction(origin_group_num,target_group_num).val=val4;
                        val=val4;
                        val333=sorty(val4);
                        [x11,y11]=size(val4);
                        imagesc(val333);title(['gr:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(x11)]);   
                        caxis([0 0.15]);
                        load('customcolormap1.mat');
                        colormap(CustomColormap);
                    end
                end
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        last00=numel(granger11(period0).direction(origin_group_num,target_group_num).val);
                        if isempty(granger11_sakura(period0).direction)
                            continue
                        end
                        for cnum=last00+1:last00+numel(granger11_sakura(period0).direction(origin_group_num,target_group_num).val)
                            if isempty(granger11_sakura(period0).direction(origin_group_num,target_group_num).val)==1
                                continue
                            end
                            granger11(period0).direction(origin_group_num,target_group_num).data(cnum).val(cnum)=...
                                 granger11_sakura(period0).direction(origin_group_num,target_group_num).val(cnum-last00);
                        end
                    end
                end
                %colormap redblue;
                colorbar;
                [figure_path0,foldername0]
                figname=['Granger_sort1',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
%                close
                % tune
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        matrix_sig{target_group_num,origin_group_num}=[];
                    end
                end
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];val1=[];
                        cnum1=0;
                        for cnum=1:numel(granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data)
                            val0(cnum).val=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger12tunebeta3(origin_group_num,target_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val1(cnum1).val=val0(cnum).val;
                        end
                        for cnum=1:numel(granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data)
                            val0(cnum).val=granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger21tunebeta3(target_group_num,origin_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val1(cnum1).val=val0(cnum).val;
                        end
                        clear mean1 mean0 sem0 
                        val2=[];mean1=[];
                        h00=0;
                        if isempty(val1)==0
                            for vnum1=1:length(val1)
                                val2(:,vnum1)=val1(vnum1).val;
                            end
                            [x1,y1]=size(val2);
                            for x0=1:x1
                                mean0(x0)=mean(val2(x0,:));
                                sem0(x0)=std(val2(x0,:))/sqrt(y1);
                            end
                            for y0=1:y1
                                mean1(y0)=mean(val2(:,y0));
                                plot(mean(freq0)-2,mean1(y0),'ok','markersize',4);hold on;                   
                            end
                            plot(freq0,mean0','-r');hold on;
                            plot(freq0,mean0'-sem0','-r');hold on;
                            plot(freq0,mean0'+sem0','-r');hold on;
                            plot([freq0(1) freq0(numel(freq0))],[0 0],'-k');
                            plot(mean(freq0)-1,mean(mean1),'+r','markersize',10);
                            plot(mean(freq0)-1,mean(mean1)+std(mean1),'+r','markersize',10);
                            plot(mean(freq0)-1,mean(mean1)-std(mean1),'+r','markersize',10);
                            h00=0;h01=0;
                            try
                                [h00,p00]=ttest(mean1,0);%ここは差の検定ではない
                            catch
                                h00=0;p00=1;
                            end
                            if isnan(h00)==0
                                if h00
                                    if mean(mean1)>0
                                        text(mean(freq0),mean(mean1)+std(mean1),['av p=',num2str(p00)]);
                                    else
                                        text(mean(freq0),mean(mean1)+std(mean1),['ap p=',num2str(p00)]);
                                    end
                                end
                            end
                            title(['av-ap:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(numel(mean1))]);            
                        end
                    end
                end
                [figure_path0,foldername0]
                figname=['Granger_tune12',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
                close;
                matrix0              
                clear myLabel myColorMap
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        val=matrix0{origin_group_num,target_group_num};
                        val2=matrix0{target_group_num,origin_group_num};
                        if origin_group_num>target_group_num
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix1(origin_group_num,target_group_num)=nan;
                            else
                                matrix1(origin_group_num,target_group_num)=mean(val);
                            end
                            
                            if sum(isnan(val2))>0 || mean(val2)<0
                                matrix3(origin_group_num,target_group_num)=nan;
                            else
                                matrix3(origin_group_num,target_group_num)=mean(val2);
                            end
                        else
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix1(origin_group_num,target_group_num)=nan;
                            else
                                matrix1(origin_group_num,target_group_num)=mean(val);
                            end
                            if sum(isnan(val2))>0 || mean(val2)<0
                                matrix3(origin_group_num,target_group_num)=nan;
                            else
                                matrix3(origin_group_num,target_group_num)=mean(val2);
                            end
                        end
                    end
                    myLabel{origin_group_num} = fsum.freq_channel(origin_group_num).group_name;
                    switch origin_group_num
                        case 1
                          myColorMap(origin_group_num,:) =[1 0 0];
                        case 2
                          myColorMap(origin_group_num,:) =[0 1 0];          
                        case 3
                          myColorMap(origin_group_num,:) =[0 0 1];  
                        case 4
                          myColorMap(origin_group_num,:) =[0.5 0.8 1];  
                        case 5
                          myColorMap(origin_group_num,:) =[1 0 1];  
                        case 6
                          myColorMap(origin_group_num,:) =[0 0.5 0.5];  
                        case 7
                          myColorMap(origin_group_num,:) =[0.8 0.8 0.1];  
                        otherwise
                          myColorMap(origin_group_num,:) =[0 0 0];  
                    end     
                end
                myLabel{4} = [myLabel{4},'granger 12']; 
                try
                    myLabel{5} = [myLabel{5},cuename,'-',periodname,'-',apavname];
                catch
                    myLabel{2} = [myLabel{2},cuename,'-',periodname,'-',apavname];
                end
                % matrix1 matrix3 matrix11 matrix4
                myLabel3=myLabel;
                try
                    myLabel3{6}  = [myLabel{6},'only sig'];
                catch
                    myLabel3{3}  = [myLabel{3},'only sig'];
                end

               
                figure;
                imagesc(matrix00);colorbar;colormap redblue;
                title(['Matrix',cuename,periodname,apavname]);
                periodmat0(period0).matrix00=matrix00;
                periodmat0(period0).matrix00max=max(matrix00);
                if period0==3
                    m1=max(periodmat0(1).matrix00max);
                    caxis([-m1 m1]);
                end
                [figure_path0,foldername0]
                figname=['Matrix',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
                close;

                periodmat0(period0).matrix2=matrix1;
                close;
                               
                figure;
                clear matrix0 matrix00
                % assym
                clear val222
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        val222(origin_group_num,target_group_num).val2=[];
                    end
                end
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val0=[];
                        val1=[];
                        cnum1=0;val0o=[];val0t=[];
                        for cnum=1:numel(granger0(period0).granger12beta3(origin_group_num,target_group_num).data)
                            val0(cnum).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect-granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect;    
                            val0(cnum).valsum=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect+granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val0o(cnum1).val=granger0(period0).granger12beta3(origin_group_num,target_group_num).data(cnum).spect;
                            val0t(cnum1).val=granger0(period0).granger21beta3(origin_group_num,target_group_num).data(cnum).spect;
                            val1(cnum1).val=val0(cnum).val./val0(cnum).valsum;
                        end
                        for cnum=1:numel(granger0(period0).granger21beta3(target_group_num,origin_group_num).data)
                            val0(cnum).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect-granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect;    
                            val0(cnum).valsum=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect+granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect;    
                            freq0=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).freq0;
                            [max00,i00]=max(val0(cnum).val);
                            cnum1=cnum1+1;
                            val0o(cnum1).val=granger0(period0).granger21beta3(target_group_num,origin_group_num).data(cnum).spect;
                            val0t(cnum1).val=granger0(period0).granger12beta3(target_group_num,origin_group_num).data(cnum).spect;
                            val1(cnum1).val=val0(cnum).val./val0(cnum).valsum;
                        end
                        cnum1
                        clear mean1 mean0 sem0 
                        val2=[];mean1=[];val2o=[];val2t=[];mean1o=[];mean0o=[];mean1t=[];mean0t=[];
                        h00=0;
                        matrix0{origin_group_num,target_group_num}=nan;
%                        matrix00{origin_group_num,target_group_num}=nan;
                        matrix0diff{origin_group_num,target_group_num}=nan;
                        if isempty(val1)==0
                            %val2=[val1(:).val];
                            for vnum1=1:length(val1)
                                val2(:,vnum1)=val1(vnum1).val;
                            end
                            [x1,y1]=size(val2);
                            for x0=1:x1
                                mean0(x0)=mean(val2(x0,:));
                                sem0(x0)=std(val2(x0,:))/sqrt(y1);
                            end
                            asym2(period0).direction(origin_group_num,target_group_num).val2=val2;
                            for y0=1:y1
                                mean1(y0)=mean(val2(:,y0));
                                plot(20,mean1(y0),'ok','markersize',4);hold on;                   
                            end
                            val222(origin_group_num,target_group_num).val2=val2;
                            plot(freq0,mean0','-b');hold on;
                            plot(freq0,mean0'-sem0','-r');hold on;
                            plot(freq0,mean0'+sem0','-r');hold on;
                            plot([freq0(1) freq0(numel(freq0))],[0 0],'-k');
                            plot(25,mean(mean1),'+','markersize',10);
                            plot(25,mean(mean1)+std(mean1),'+','markersize',10);
                            plot(25,mean(mean1)-std(mean1),'+','markersize',10);
                            matrix_asym(origin_group_num,target_group_num)=0;
                            try
                                [h00,p00]=ttest(mean1,0);
                            catch
                                h00=0;p00=1;
                            end
                            if ~isnan(h00)
                                if h00
                                    plot([beta_range(1) beta_range(2)],[mean(mean1) mean(mean1)],'-k');
                                    text(25,mean(mean1)+std(mean1),['p=',num2str(p00)]);
                                    matrix_asym(origin_group_num,target_group_num)=mean(mean(val2));
                                end
                            end
%                            matrix11{origin_group_num,target_group_num}.val2=val2;
                            title(['asym:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(y1)]);

                            for vnum1=1:length(val0o)
                                val2o(:,vnum1)=val0o(vnum1).val;
                            end
                            [x1o,y1o]=size(val2o);
                            clear mean0o sem0o;
                            for x0=1:x1o
                                mean0o(x0)=mean(val2o(x0,:));
                                sem0o(x0)=std(val2o(x0,:))/sqrt(y1o);
                            end
                            for y0=1:y1o
                                mean1o(y0)=mean(val2o(:,y0));
                                plot(20,mean1o(y0),'ok','markersize',4);hold on;                   
                            end
                            plot(freq0,mean0o','-b');hold on;
                            plot(freq0,mean0o'-sem0o','-r');hold on;
                            plot(freq0,mean0o'+sem0o','-r');hold on;
                            clear mean0t sem0t;
                            for vnum1=1:length(val0t)
                                val2t(:,vnum1)=val0t(vnum1).val;
                            end
                            [x1t,y1t]=size(val2t);
                            for x0=1:x1t
                                mean0t(x0)=mean(val2t(x0,:));
                                sem0t(x0)=std(val2t(x0,:))/sqrt(y1t);
                            end
                            for y0=1:y1t
                                mean1t(y0)=mean(val2t(:,y0));
                                plot(20,mean1t(y0),'ok','markersize',4);hold on;                   
                            end
                            plot(freq0,mean0t','-b');hold on;
                            plot(freq0,mean0t'-sem0t','-r');hold on;
                            plot(freq0,mean0t'+sem0t','-r');hold on;

                            try
                                [hot,pot]=ttest(mean1o,mean1t);
                            catch
                                hot=0;pot=1;
                            end         
                        end
                        val3=mean(mean1);
                        if ~isnan(h00)
                            if h00
                                matrix0{origin_group_num,target_group_num}=mean1;
                                matrix00{origin_group_num,target_group_num}=mean1;
                            else
                                matrix00{origin_group_num,target_group_num}=mean1;
                            end
                        end
                        if ~isnan(hot)
                            if hot
                                matrix0diff{origin_group_num,target_group_num}=mean1;
                            end
                        end
                    end
                end
                matrix0
                [figure_path0,foldername0]
                figname=['Asym_net4',cuename,periodname];
                print_file10([figure_path0,foldername0],figname,1,1);
                close;
                h=figure;h.Position(3:4) = [4*280 4*210];
                imagesc(matrix_asym);colorbar;colormap redblue;
                title(['Matrix-asym',cuename,periodname,apavname]);
                periodmat0(period0).matrix_asym=matrix_asym;
                periodmat0(period0).matrix_asym_max=max(matrix_asym);
                if period0==3
                    m1=max(periodmat0(1).matrix_asym_max);
                    caxis([-m1 m1]);
                end
                [figure_path0,foldername0]
                figname=['Matrix_asym',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
                close;
                %figure;
                h=figure;
                h.Position(3:4) = [4*280 4*210];
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        subplot(4,4,target_group_num+4*(origin_group_num-1));
                        origin_region=fsum.freq_channel(origin_group_num).group_name;
                        target_region=fsum.freq_channel(target_group_num).group_name;
                        val4=(val222(origin_group_num,target_group_num).val2)';
                        asym0(period0).direction(origin_group_num,target_group_num).val=val4;
                        val33=sorty(val4);
                        [x12,y12]=size(val33);
                        imagesc(val33);title(['asym:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(x12)]);     
                        caxis([-1 1]);
                    end
                end

                colormap redblue;colorbar;
                [figure_path0,foldername0]
                figname=['Asym_graph1',cuename,periodname,apavname];
                print_file10([figure_path0,foldername0],figname,1,1);
                colorbar;
                close
                clear matrix1 matrix2
                clear myLabel myColorMap
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        val=matrix0{origin_group_num,target_group_num};
                        if origin_group_num>target_group_num               
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix1(origin_group_num,target_group_num)=0;
                            else
                                matrix1(origin_group_num,target_group_num)=mean(val);
                            end
                        else
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix2(origin_group_num,target_group_num)=0;
                            else
                                matrix2(origin_group_num,target_group_num)=mean(val);
                            end
                        end
                    end
                    myLabel{origin_group_num} = fsum.freq_channel(origin_group_num).group_name;
                    switch origin_group_num
                        case 1
                          myColorMap(origin_group_num,:) =[1 0 0];
                        case 2
                          myColorMap(origin_group_num,:) =[0 1 0];          
                        case 3
                          myColorMap(origin_group_num,:) =[0 0 1];  
                        case 4
                          myColorMap(origin_group_num,:) =[0.5 0.8 1];  
                        case 5
                          myColorMap(origin_group_num,:) =[1 0 1];  
                        case 6
                          myColorMap(origin_group_num,:) =[0 0.5 0.5];  
                        case 7
                          myColorMap(origin_group_num,:) =[0.8 0.8 0.1];  
                        otherwise
                          myColorMap(origin_group_num,:) =[0 0 0];  
                    end     
                end
                myLabel{4} = [myLabel{4},'granger asym sig']; 
                try
                    myLabel{5} = [myLabel{5},cuename,'-',periodname,'-',apavname];
                catch
                    myLabel{2} = [myLabel{2},cuename,'-',periodname,'-',apavname];
                end
                matrix1(isnan(matrix1))=0;
                matrix2(isnan(matrix2))=0;
                [figure_path0,foldername0]

%                 figname=['Asym_diff1',cuename,periodname,apavname];
%                 print_file10([figure_path0,foldername0],figname,1,1);

                clear matrix1diff matrix2diff
                clear myLabeldiff myColorMapdiff
                for origin_group_num=1:4%numel(fsum.freq_channel)
                    for target_group_num=1:4%numel(fsum.freq_channel)
                        val=matrix0diff{origin_group_num,target_group_num};
                        if origin_group_num>target_group_num               
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix1diff(origin_group_num,target_group_num)=0;
                            else
                                matrix1diff(origin_group_num,target_group_num)=mean(val);
                            end
                        else
                            if sum(isnan(val))>0 || mean(val)<0
                                matrix2diff(origin_group_num,target_group_num)=0;
                            else
                                matrix2diff(origin_group_num,target_group_num)=mean(val);
                            end
                        end
                    end
                    myLabeldiff{origin_group_num} = fsum.freq_channel(origin_group_num).group_name;
                    switch origin_group_num
                        case 1
                          myColorMapdiff(origin_group_num,:) =[1 0 0];
                        case 2
                          myColorMapdiff(origin_group_num,:) =[0 1 0];          
                        case 3
                          myColorMapdiff(origin_group_num,:) =[0 0 1];  
                        case 4
                          myColorMapdiff(origin_group_num,:) =[0.5 0.8 1];  
                        case 5
                          myColorMapdiff(origin_group_num,:) =[1 0 1];  
                        case 6
                          myColorMapdiff(origin_group_num,:) =[0 0.5 0.5];  
                        case 7
                          myColorMapdiff(origin_group_num,:) =[0.8 0.8 0.1];  
                        otherwise
                          myColorMapdiff(origin_group_num,:) =[0 0 0];  
                    end     
                end
                myLabeldiff{4} = [myLabeldiff{4},'granger asym sig']; 
                try
                    myLabeldiff{5} = [myLabeldiff{5},cuename,'-',periodname,'-',apavname];
                catch
                    myLabeldiff{2} = [myLabeldiff{2},cuename,'-',periodname,'-',apavname];
                end
                matrix1diff(isnan(matrix1diff))=0;
                matrix2diff(isnan(matrix2diff))=0;
                [figure_path0,foldername0]

                figure;

                periodmat0(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).matrix=matrix0;
                periodmat0(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).matrix00=matrix00;
            end
%%%%    
            asym0 
            close;

            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    if isempty(asym0(1).direction)==1
                        continue
                    end
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    val4=asym0(1).direction(origin_group_num,target_group_num).val;                    
                    val5=asym0(3).direction(origin_group_num,target_group_num).val;
                    val333=sorty(val5-val4);
                    [x13,y13]=size(val333);
                    imagesc(val333);title(['asym:',cuename,periodname,apavname,':',origin_region,'->',target_region,':',num2str(y13)]);  colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym1',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
            close
            %asym2           
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(asym2(1).direction)==1
                        continue
                    end
                    val4=asym2(1).direction(origin_group_num,target_group_num).val2;                    
                    val5=asym2(3).direction(origin_group_num,target_group_num).val2;
                    calc_mean1(val4,'-b');hold on;calc_mean1(val5,'-r');hold on;
                    calc_diffmean1(val4,val5,'-b','-r',0.05);
                    calc_mean1(val4,'-b');hold on;calc_mean1(val5,'-r');hold on;
                    if (origin_group_num==4 && target_group_num==3) || (origin_group_num==3 && target_group_num==4)
                        axis([1 26 0 0.5]);
                    elseif (origin_group_num==4 && target_group_num==2) || (origin_group_num==2 && target_group_num==4)
                        axis([1 26 -0.2 0.3]);
                    elseif (origin_group_num==4 && target_group_num==1) || (origin_group_num==1 && target_group_num==4)
                        axis([1 26 -0.1 0.4]);
                    elseif (origin_group_num==3 && target_group_num==1) || (origin_group_num==1 && target_group_num==3)
                        axis([1 26 -0.1 0.2]);
                    elseif (origin_group_num==3 && target_group_num==2) || (origin_group_num==2 && target_group_num==3)
                        axis([1 26 -0.2 0.3]);
                    elseif (origin_group_num==2 && target_group_num==1) || (origin_group_num==1 && target_group_num==2)
                        axis([1 26 -0.2 0.4]);
                    end
                    
                    %calc_hist1(val4,val5);
                    %imagesc(val333);
                    title(['asym2:',cuename,'-',periodname,'-',apavname,':',origin_region,'->',target_region]);  colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
           
            figname=['Diff3-1_asym2',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
            
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(asym2(1).direction)==1
                        continue
                    end
                    val4=asym2(1).direction(origin_group_num,target_group_num).val2;                    
                    val5=asym2(3).direction(origin_group_num,target_group_num).val2;
                    %calc_mean1(val5,'-r');hold on;calc_mean1(val4,'-b');
                    calc_hist1(val4,val5);
                    %imagesc(val333);
                    title(['asym4:',cuename,'-',periodname,'-',apavname,':',origin_region,'->',target_region]);  
                    colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym3',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
%            close
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(granger11(1).direction)
                        continue
                    end
                    val4=granger11(1).direction(origin_group_num,target_group_num).val;
                    val5=granger11(1).direction(target_group_num,origin_group_num).val;       
                    val6=granger11(3).direction(origin_group_num,target_group_num).val;
                    val7=granger11(3).direction(target_group_num,origin_group_num).val;                    
                    [h,p,ci,stats] = calc_hist1_ttest(val4,val5,val6,val7);
                     if h==1
                         title({['asym*:p=',num2str(p)];[cuename,'-',periodname,'-',apavname,':',origin_region,'->',target_region,':',num2str(numel(val7))]});
                     else
                         title(['asym:',cuename,'-',periodname,'-',apavname,':',origin_region,'->',target_region,':',num2str(numel(val7))]);
                     end
                    colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym4',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(asym0(1).direction)
                        continue
                    end
                    val4=asym0(1).direction(origin_group_num,target_group_num).val;                    
                    val5=asym0(3).direction(origin_group_num,target_group_num).val;
                    val333=sorty(val5-val4);
                    [x19,y19]=size(val333);
                    imagesc(val333);title(['asym:',cuename,periodname,apavname,':',origin_region,'>',target_region,':',num2str(x19)]);  colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym5',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
%            close
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(asym0(1).direction)
                        continue
                    end
                    val4=asym0(1).direction(origin_group_num,target_group_num).val;                    
                    val5=asym0(3).direction(origin_group_num,target_group_num).val;
                    %calc_mean2(val5,'-r');hold on;calc_mean2(val4,'-b');
                    %calc_diffmean2(val4,val5,'-b','-r');
                    [val39,y09,y19]=calc_hist1(val4,val5);
                    %imagesc(val333);
                    title(['asym:',cuename,periodname,apavname,':',origin_region,'>',target_region,':',num2str(y09),':',num2str(numel(x19))]);  
                    colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym6',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
            h=figure;
            h.Position(3:4) = [4*280 4*210];
            for origin_group_num=1:4%numel(fsum.freq_channel)
                for target_group_num=1:4%numel(fsum.freq_channel)
                    subplot(4,4,target_group_num+4*(origin_group_num-1)); 
                    origin_region=fsum.freq_channel(origin_group_num).group_name;
                    target_region=fsum.freq_channel(target_group_num).group_name;
                    if isempty(asym0(1).direction)
                        continue
                    end
                    val4=asym0(1).direction(origin_group_num,target_group_num).val;                    
                    val5=asym0(3).direction(origin_group_num,target_group_num).val;
                    %calc_mean2(val5,'-r');hold on;calc_mean2(val4,'-b');
                    %calc_diffmean2(val4,val5,'-b','-r');
                    [val39,y09,y19]=calc_hist2(val4,val5);
                    %imagesc(val333);
                    title(['asym:',cuename,periodname,apavname,':',origin_region,'>',target_region,':',num2str(y09),':',num2str(numel(x19))]);  
                    colormap;   
                end
            end
            colormap redblue;colorbar;
            [figure_path0,foldername0]
            figname=['Diff3-1_asym7',cuename,periodname,apavname];
            print_file10([figure_path0,foldername0],figname,1,1);
        end
        save('mymaxy.mat','maxy');
    end
    % granger beta
    nnum1
    nnum2
    clear val0 val1;
    for cue1precue2all3=cue1precue2all3_loop 
        figure;
        for origin_group_num=1:4%numel(fsum.freq_channel)            
            for target_group_num=1:4%numel(fsum.freq_channel)
                subplot(4,4,target_group_num+4*(origin_group_num-1));
                clear bar0 bar1;
                for ap1av2all3=ap1av2all3_loop
                    for period0=period00
                        try
                            val0=periodmat0(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).matrix00{origin_group_num,target_group_num};
                        catch
                            val0=nan;
                        end
                        try
                            val1=periodmat0(period0).cue1precue2all3(cue1precue2all3).appro1_avoid2_all3(ap1av2all3).matrix{origin_group_num,target_group_num};
                        catch
                            val1=nan;
                        end
                        val0(isnan(val0))=0; 
                        bar0(period0,ap1av2all3)=mean(val0);
                        switch ap1av2all3
                            case 1
                                plot(period0,mean(val1),'bo');hold on;
                            case 2
                                plot(period0,mean(val1),'ro');hold on;
                            case 3
                                plot(period0,mean(val1),'ko');hold on;
                        end
                    end
                    switch ap1av2all3
                        case 1
                            plot([1 3],[bar0(1,1),bar0(3,1)],'-b');hold on;
                        case 2
                            plot([1 3],[bar0(1,2),bar0(3,2)],'-r');hold on;
                        case 3
                            plot([1 2],[bar0(1,3),bar0(2,3)],'-k');hold on;
                            try
                                bar0(3,3)
                            catch
                                continue
                            end
                            plot([2 3],[bar0(2,3),bar0(3,3)],'-k');hold on;
                            plot(period0,mean(val1),'ko');hold on;
                    end
                end
%                bar(bar0);hold on;
            end
        end
        switch cue1precue2all3
            case 1
                cuename='cue';
            case 2
                cuename='precue';
            case 3
                cuename='allperiod';
        end
        title(cuename);
    end
    nnum00=0;
    clear channel1val myLabel1 matrix0 matrix2;
    onum0=0;
%    close all;
    % [figure_path0,foldername0]
    % figname=[cuename,periodname,'fig13'];
    % print_file10([figure_path0,foldername0],figname,1,1); 
end

% periodmat0(1).matrix_final{1,1}=14;
% periodmat0(3).matrix_final{1,1}=25;
for nnum1=1:4
    for nnum2=1:4
        if isnan(periodmat0(1).matrix_final(nnum1,nnum2))
            periodmat0(1).matrix_final(nnum1,nnum2)=0;
        else
            periodmat0(1).matrix_final(nnum1,nnum2)=periodmat0(1).matrix_final(nnum1,nnum2)*1000;
        end
        if isnan(periodmat0(3).matrix_final(nnum1,nnum2))
            periodmat0(3).matrix_final(nnum1,nnum2)=0;
        else
            periodmat0(3).matrix_final(nnum1,nnum2)=periodmat0(3).matrix_final(nnum1,nnum2)*1000;
        end
    end
end

periodmat0(1).matrix_final(isnan(periodmat0(1).matrix3))=0;
periodmat0(3).matrix_final(isnan(periodmat0(3).matrix3))=0;
periodmat0(1).matrix_final(1,1)=19*1000;
periodmat0(3).matrix_final(1,1)=25*1000;

return
