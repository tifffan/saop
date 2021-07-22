function [filts] = spm_parseval_filters(G,num_filters,filter_type)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
switch filter_type
    
    % Uniform filters adapted to interval
    case 'unif_ideal'
        param.no_shift=1;
        param.band_structure=0;
        param.spectrum_adapted=0;
        param.plot_filters=0;
        param.plot_density_functions=0;
        param.spacing=0;
        filts=mcsfb_design_filter_bank(G,num_filters,param);
    case 'unif_meyer'
        filts=spg_filter_design(G.lmax,num_filters,'designtype','uniform_meyer_type');
    case 'unif_itersine'
        filts = gsp_design_itersine(G.lmax,num_filters);
    case 'uniform_hann'
        filts=gsp_design_half_cosine(G.lmax,num_filters);
        norm_factor=sqrt(filts{1}(0)^2+filts{2}(0)^2);
        for i=1:num_filters
            filts{i}=@(x) filts{i}(x)/norm_factor;
        end
    case 'unif_DCT'
        h=myDCT(2,num_filters);
        h=h/sqrt(num_filters);
        filts = gfb_filter_design(h,G.lmax, num_filters)';
    
    % Wavelet filters adapted to interval
    case 'wav_ideal'
        param.no_shift=1;
        param.band_structure=0;
        param.spectrum_adapted=0;
        param.plot_filters=0;
        param.plot_density_functions=0;
        param.spacing=1;
        filts=mcsfb_design_filter_bank(G,num_filters,param);
    case 'wav_meyer'
        filts=spg_filter_design(G.lmax,num_filters,'designtype','meyer');
    case 'wav_itersine'
        param.log=1;
        param.warping_type='custom';
        param.filter = 'itersine';
        wf=@(x)x;
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
    case 'wav_hann'
        param.log=1;
        param.warping_type='custom';
        param.filter = 'half_cosine';
        wf=@(x)x;
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
        norm_factor=sqrt(filts{1}(0)^2+filts{2}(0)^2);
        for i=1:num_filters
            filts{i}=@(x) filts{i}(x)/norm_factor;
        end
    case 'wav_DCT'
        h=myDCT(2,num_filters);
        h=h/sqrt(num_filters);
        filt_dct = gfb_filter_design(h,G.lmax, num_filters)';
        param.log=1;
        param.warping_type='custom';
        param.filter = filt_dct;
        wf=@(x)x;
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
    case 'wav_fast'
        FrameType='Haar'; 
        DFilters=ExtractMasks(FrameType);
        r=length(DFilters);
        Lev=num_filters-1;
        s=2;
        J=log(G.lmax/pi)/log(s)+Lev-1;
        base_filts=cell(Lev-1,1);
        base_filts{1}=@(x)DFilters{1}(2^(-J)*x);
        for k=2:(Lev-1)
            base_filts{k}=@(x)base_filts{k-1}(x).*DFilters{1}(2^(-J+k-1)*x);
        end
        filts=cell((r-1)*Lev+1,1);
        % scaling
        filts{1}=@(x) DFilters{1}(2^(-J+Lev-1)*x).*base_filts{Lev-1}(x);
        % bandpass
        for i=2:((r-1)*Lev)
            filts{i}=@(x) DFilters{2}(2^(-J+Lev-i+1)*x).*base_filts{Lev+1-i}(x);
        end
        % highest
        filts{(r-1)*Lev+1}=@(x)DFilters{2}(2^(-J)*x);

    % Uniform filters adapted to spectral density    
    case 'sa_unif_ideal'
        param.no_shift=1;
        param.band_structure=0;
        param.plot_filters=0;
        param.plot_density_functions=0;
        param.spacing=0;
        param.spectrum_adapted=1;
        filts=mcsfb_design_filter_bank(G,num_filters,param);
    case 'sa_unif_meyer'
        filt_meyer_uniform=spg_filter_design(G.lmax,num_filters,'designtype','uniform_meyer_type');
        param.log=0;
        param.warping_type='custom';
        param.filter=filt_meyer_uniform;
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
    case 'sa_unif_itersine'
        param.log=0;
        param.warping_type='custom';
        param.filter = 'itersine';
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
    case 'sa_unif_hann'
        param.log=0;
        param.warping_type='custom';
        param.filter = 'half_cosine';
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
        norm_factor=sqrt(filts{1}(0)^2+filts{2}(0)^2);
        for i=1:num_filters
            filts{i}=@(x) filts{i}(x)/norm_factor;
        end
        
    % Wavelet filters adapted to spectral density 
    case 'sa_wav_ideal'
        param.no_shift=1;
        param.band_structure=0;
        param.plot_filters=0;
        param.plot_density_functions=0;
        param.spectrum_adapted=1;
        param.spacing=1;
        filts=mcsfb_design_filter_bank(G,num_filters,param);
    case 'sa_wav_meyer' %TODO: double check this one
        filt_meyer_uniform=spg_filter_design(G.lmax,num_filters,'designtype','uniform_meyer_type');
        param.log=1;
        param.warping_type='custom';
        param.filter=filt_meyer_uniform;
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        filts=gsp_design_warped_translates(G,num_filters,param);
    case 'sa_wav_itersine'
        param.log=1;
        param.warping_type='custom';
        param.filter = 'itersine';
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        [filts,~] = gsp_design_warped_translates(G,num_filters,param);
    case 'sa_wav_hann'
        param.log=1;
        param.warping_type='custom';
        param.filter = 'half_cosine';
        wf=@(x)G.lmax*G.spectrum_cdf_approx(x);
        param.warp_function=wf;
        filts = gsp_design_warped_translates(G,num_filters,param);
        norm_factor=sqrt(filts{1}(0)^2+filts{2}(0)^2);
        for i=1:num_filters
            filts{i}=@(x) filts{i}(x)/norm_factor;
        end
    otherwise 
        error('Unknown filter type');
end
end

