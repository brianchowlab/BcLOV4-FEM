function [] = StepwiseInterpolate(pde_C,sol_M,props,param,conversion,I)
    psf_params.size = [64,64,64];%[32,32,21]
    psf_params.NA = param.NA;
    psf_params.lambda = 610e-9;
    psf_params.M = param.mag;
    psf_params.ti0 = 100e-6;
    psf_params.resLateral = param.scale_len * 1e-6;%scale_len * 1e-6;
    psf_params.resAxial = param.scale_len * 1e-6;%scale_len * 1e-6;
    psf_params.pZ = 0;%ceil(params.size(3)/2) * axial_resolution*1e-6;
    psf_params.oversampling = 2;
    [PSF_3D] = GenPSF(psf_params,param);

    for i=1:param.interpolation_interval:size(sol_M,2)
        i
        pdem_C = createpde(1);
        gm_C_f = geometryFromMesh(pdem_C,props.nodes',props.elements');
        if i == size(sol_M,2)
            pde_C_i=createPDEResults(pdem_C,pde_C.NodalSolution(:,[i,1]),1:2,'time-dependent');
        else
            pde_C_i=createPDEResults(pdem_C,pde_C.NodalSolution(:,i:i+1),1:2,'time-dependent');
        end
        
        c_intrp = InterpolateCytoplasm(pde_C_i,1,I,param);

        m_intrp = InterpolateMembrane(I,1,props.pm_TR,sol_M(:,i),param);
        m_intrp(isnan(m_intrp)) = 0;


        c_intrp = c_intrp/conversion * param.conc_ratio;
        m_intrp = m_intrp * 4;
        c_intrp(m_intrp ~= 0) = m_intrp(m_intrp~=0);
        c_intrp_blurred = ConvolvePSF(c_intrp,single(PSF_3D));
        imwrite(uint16(squeeze(c_intrp_blurred(:,:,ceil(size(c_intrp,3)/2)))),['./Out_DM_S10000/',num2str(i),'.png']);
    end
end