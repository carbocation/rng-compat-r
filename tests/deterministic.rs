use rng_compat_r::{pnorm_with_mode, MathMode, RRng};

#[test]
fn deterministic_normal_bits_are_pinned() {
    let mut rng = RRng::from_seed(42).with_math_mode(MathMode::Deterministic);
    let actual: Vec<_> = (0..32)
        .map(|iteration| {
            let _ = rng.sample_index(3 + iteration);
            rng.rnorm(-1.25, 2.5).to_bits()
        })
        .collect();
    let expected = [
        0x4004_9d11_4380_13f3,
        0xbfd5_e642_a129_3b16,
        0xc010_0bc6_4a19_22e2,
        0xbff7_c94e_fd5a_33d5,
        0xc010_e05c_ab42_45a6,
        0x4000_18ee_c0d0_a769,
        0x4006_4477_ff33_8ae0,
        0xbfff_26cc_065e_e944,
        0x3fd5_c087_3efe_43a0,
        0xbfa6_c72a_f1aa_1a80,
        0xc01d_6798_6804_03a5,
        0xbfe1_488b_4bd1_b661,
        0xc016_d026_4b98_54ac,
        0x4008_76f1_1280_8da4,
        0xc002_9c00_800b_b913,
        0x3fe7_0d1e_b566_1978,
        0xbfb9_89a6_3e0f_17a0,
        0xbfbc_830d_434c_aaf0,
        0xbfea_9acf_6e6e_6882,
        0xc006_2db3_fa90_ea72,
        0xc008_174d_90fe_94ab,
        0xc009_b06e_1c88_d9bb,
        0x3fc9_e5b3_809a_6ec8,
        0xbff2_8e1a_c0de_0f57,
        0xc00e_0fe5_3391_9381,
        0x3fe4_a72e_dd7c_b79a,
        0xc012_68a3_dded_d98f,
        0xbfc5_7f8c_1d1c_6088,
        0x4003_f9bb_bb62_ecc0,
        0xc002_a101_2aac_2778,
        0xbfdc_7df0_ee5e_de1e,
        0x4005_83b9_9322_4653,
    ];
    assert_eq!(actual, expected);
}

#[test]
fn deterministic_pnorm_bits_are_pinned() {
    let values = [-50.0, -10.0, -8.0, -1.0, 0.0, 1.0, 8.0, 10.0, 50.0];
    let actual: Vec<_> = [false, true]
        .into_iter()
        .flat_map(|log_p| {
            [true, false].into_iter().flat_map(move |lower_tail| {
                values.into_iter().map(move |x| {
                    pnorm_with_mode(x, 0.0, 1.0, lower_tail, log_p, MathMode::Deterministic)
                        .to_bits()
                })
            })
        })
        .collect();
    let expected = [
        0x0000_0000_0000_0000,
        0x3b22_6c75_e84f_b10e,
        0x3cc6_69d2_c90d_55cf,
        0x3fc4_4ed0_bb7c_b20b,
        0x3fe0_0000_0000_0000,
        0x3fea_ec4b_d120_d37d,
        0x3fef_ffff_ffff_fffa,
        0x3ff0_0000_0000_0000,
        0x3ff0_0000_0000_0000,
        0x3ff0_0000_0000_0000,
        0x3ff0_0000_0000_0000,
        0x3fef_ffff_ffff_fffa,
        0x3fea_ec4b_d120_d37d,
        0x3fe0_0000_0000_0000,
        0x3fc4_4ed0_bb7c_b20b,
        0x3cc6_69d2_c90d_55cf,
        0x3b22_6c75_e84f_b10e,
        0x0000_0000_0000_0000,
        0xc093_9b53_5055_a3e5,
        0xc04a_9d9a_c076_c031,
        0xc041_81b8_4f11_312b,
        0xbffd_74d3_1cc8_afc1,
        0xbfe6_2e42_fefa_39ef,
        0xbfc6_1ccb_bb95_43af,
        0xbcc6_69d2_c90d_55d1,
        0xbb22_6c75_e84f_b10e,
        0x8000_0000_0000_0000,
        0x8000_0000_0000_0000,
        0xbb22_6c75_e84f_b10e,
        0xbcc6_69d2_c90d_55d1,
        0xbfc6_1ccb_bb95_43af,
        0xbfe6_2e42_fefa_39ef,
        0xbffd_74d3_1cc8_afc1,
        0xc041_81b8_4f11_312b,
        0xc04a_9d9a_c076_c031,
        0xc093_9b53_5055_a3e5,
    ];
    assert_eq!(actual, expected);
}
