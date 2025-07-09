function nonlcon(OPT::OPT_struct, FE::FE_struct,SGs::SGs_struct, dv)

    OPT.dv_old = copy(OPT.dv) # save the previous design
    OPT.dv = copy(dv) # update the design

    # Update or perform the analysis only if design has changed
    OPT, FE = perform_analysis(OPT, FE)
    #-------------------------------
    g, gradg = eval(Symbol(OPT.functions.constraints))(OPT, FE, SGs)

    if OPT.functions.constraints == "volume_fraction" || OPT.functions.constraints == "volume_fraction_multi_scales"
        g = copy(g) - OPT.functions.volume_fraction_limit
    end
    return g, gradg, OPT, FE
end