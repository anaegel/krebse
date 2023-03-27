-- My variables
PARAMS = {
    --diffusion rate meters*day^-1
    d_a = 50.0,
    d_b = 5000.0,
    -- growth rates by species
    r_a = -1,
    r_b = 1,
    -- capacity weight
    k_a = 2, -- to 1/0.8
    k_b = 0.1, -- to 1/12
    -- capacity max,
    k_max=1.0,
    -- death rate
    m_a = 0.0, -- not applicable
    m_b = 0.03,
    -- transmission rate _ fromTo
--    i_ab = 0.0001,
--    i_bb = 0.0001,
    i_ab = 0,
    i_bb = 0,
}

