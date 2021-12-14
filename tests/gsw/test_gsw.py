from gsw import *
from hypothesis import given, target, strategies as st

# STRATEGIES

# Salinity
global_salinity_st = st.floats(min_value=0, max_value=42, allow_nan=False)
absolute_salinity_st = global_salinity_st  # Paper [0, 42]
practical_salinity_st = global_salinity_st  # Paper [0, 42]
preformed_salinity_st = global_salinity_st  # Guess
knudsen_salinity_st = global_salinity_st  # Guess
reference_salinity_st = global_salinity_st  # Guess
bulk_absolute_salinity_st = global_salinity_st  # Guess
absolute_salinity_of_sea_ice_st = (
    global_salinity_st  # Guess, probably less than sea water?
)

# Temperature
global_temperature_st = st.floats(min_value=-2, max_value=40, allow_nan=False)
global_temperature_ice_st = st.floats(min_value=-20, max_value=0, allow_nan=False)
potential_temperature_st = global_temperature_st  # Paper? [-2, 40]
in_situ_temperature_st = global_temperature_st  # Paper? [-2, 40]
conservative_temperature_st = global_temperature_st  # Paper? [-2, 40]
potential_temperature_with_reference_to_zero_pressure_st = (
    global_temperature_st  # Guess
)
temperature_t68_st = global_temperature_st  # Guess

in_situ_temperature_of_ice_st = global_temperature_ice_st  # Guess [-20, -2]
in_situ_temperature_of_sea_ice_at_pressure_st = (
    global_temperature_ice_st  # Guess [-20, -2]
)
potential_temperature_of_ice_with_reference_to_zero_pressure_st = (
    global_temperature_ice_st  # Guess [-20, -2]
)

# Pressure
global_pressure_st = st.floats(min_value=0, max_value=11000, allow_nan=False)
sea_pressure_st = global_pressure_st  # Paper [0, 11000]
reference_pressure_st = global_pressure_st  # Guess
upper_sea_pressure_st = global_pressure_st  # Guess
lower_sea_pressure_st = global_pressure_st  # Guess
geopotential_at_zero_sea_pressure_st = global_pressure_st  # Guess

# Fractions
global_fractions_st = st.floats(min_value=0, max_value=1, allow_nan=False)
saturation_fraction_st = global_fractions_st  # Documentation [0, 1]
mass_fraction_of_ice_st = global_fractions_st  # Documentation [0, 1]
mass_fraction_of_sea_ice_st = global_fractions_st  # Documentation [0, 1]

# Others
in_situ_seawater_density_rho_st = st.floats(
    min_value=1000, max_value=1060, allow_nan=False
)  # Paper? In situ density?
longitude_st = st.floats(
    min_value=-360, max_value=360, allow_nan=False
)  # Documentation
latitude_st = st.floats(min_value=-90, max_value=90, allow_nan=False)  # Documentation
conductivity_st = st.floats(
    min_value=0, max_value=60, allow_nan=False
)  # https://github.com/TEOS-10/GSW-Python/issues/83#issuecomment-976903688
interpolation_algorithm_st = st.sampled_from(["pchip", "linear"])  # Documentation
geostrophic_stream_function_st = st.sampled_from(
    [
        "geo_strf_dyn_height",
        # "geo_strf_Montgomery",
        # "geo_strf_Cunninhgam",
        # "geo_strf_isopycnal",
    ]
)  # Documentation (only geo_strf_dyn_height has been implemented)
depth_st = st.floats(min_value=0, max_value=10000, allow_nan=False)  # Documentation

# Energy SKIPPED
# specific_enthalpy_st = None
# bulk_enthalpy_st = None
# bulk_potential_enthalpy_st = None
# potential_enthalpy_of_ice_st = None
# specific_entropy_st = None

# Unknown SKIPPED
# axis_st = None
# Rt_st = None
# max_dp_st = None
# shape_st = None
# order_st = None
# f_st = None
# x_st = None
# y_st = None
# xi_st = None
# dynamic_height_anomoly_st = None

# VALID OUTPUT RANGES

# Used more than once
practical_salinity_al_range = (0, 500)  # Gaet'ale Pond has highest salinity of 433g/kg
absolute_salinity_final_al_range = practical_salinity_al_range
absolute_salinity_al_range = practical_salinity_al_range

conservative_temperature_al_range = (
    -12,
    60,
)  # Lowest is -2 ish, highest is red sea 54 ish
conservative_temperature_final_al_range = conservative_temperature_al_range
potential_temperature_al_range = conservative_temperature_final_al_range
potential_temperature_ice_al_range = (-30, 10)  # Glaciers can be like -20 ish?

# Used more than once SKIPPED
w_ih_final_al_range = None
alpha_al_range = None
beta_al_range = None
specific_entropy_al_range = None
pressure_midpoint_array_range = None

# Used twice SKIPPED
ct_freeze_al_range = None
ct_freezing_al_range = None
ct_freezing_p_al_range = None
ct_freezing_sa_al_range = None
dct_dp_frazil_al_range = None
dsa_dct_frazil_al_range = None
dsa_dp_frazil_al_range = None
h_ct_al_range = None
h_ct_ct_al_range = None
h_sa_al_range = None
h_sa_ct_al_range = None
h_sa_sa_al_range = None
latentheat_evap_al_range = None
melting_ice_sa_ct_ratio_al_range = None
melting_seaice_equilibrium_sa_ct_ratio_al_range = None
melting_seaice_sa_ct_ratio_al_range = None
o2_solar_al_range = None
p_ref_al_range = None
pot_enthalpy_ice_al_range = None
pot_enthalpy_ice_freezing_al_range = None
pot_enthalpy_ice_freezing_p_al_range = None
pot_enthalpy_ice_freezing_sa_al_range = None
rho_al_range = None
rho_sa_al_range = None
rho_sa_sa_al_range = None
sa_freeze_al_range = None
spec_vol_al_range = None
sstar_al_range = None
t_al_range = None
t_freezing_al_range = None
t_freezing_p_al_range = None
t_freezing_sa_al_range = None

# Used once SKIPPED
adiabatic_lapse_rate_al_range = None
adiabatic_lapse_rate_ice_al_range = None
alpha_on_beta_al_range = None
alpha_wrt_t_exact_al_range = None
alpha_wrt_t_ice_al_range = None
beta_const_t_exact_al_range = None
c_al_range = None
c_p_al_range = None
cabbeling_al_range = None
chem_potential_water_d_al_range = None
chem_potential_water_ice_al_range = None
chem_potential_water_t_exact_al_range = None
cp_ice_al_range = None
cp_t_exact_al_range = None
ct_maxdensity_al_range = None
ct_multiple_al_range = None
ct_p_wrt_t_al_range = None
ct_pt_al_range = None
ct_pt_pt_al_range = None
ct_sa_al_range = None
ct_sa_pt_al_range = None
ct_sa_sa_al_range = None
ct_sa_wrt_t_al_range = None
ct_t_wrt_t_al_range = None
delta_sa_al_range = None
delta_sa_atlas_al_range = None
dilution_coefficient_t_exact_al_range = None
distance_between_points_range = None
dynamic_enthalpy_al_range = None
dynamic_height_array_range = None
enthalpy_al_range = None
enthalpy_ct_exact_al_range = None
enthalpy_diff_al_range = None
enthalpy_ice_al_range = None
enthalpy_t_exact_al_range = None
eta_ct_al_range = None
eta_ct_ct_al_range = None
eta_sa_al_range = None
eta_sa_ct_al_range = None
eta_sa_sa_al_range = None
f_delta_al_range = None
gibbs_ice_part_pt0_al_range = None
gibbs_ice_part_t_al_range = None
gibbs_ice_pt0_pt0_al_range = None
grav_al_range = None
helmholtz_energy_ice_al_range = None
hill_ratio_al_range = None
ice_entropy_al_range = None
internal_energy_al_range = None
internal_energy_ice_al_range = None
ipv_vs_fn_squared_range = None
its_90_range = None
kappa_al_range = None
kappa_const_t_ice_al_range = None
kappa_ice_al_range = None
kappa_t_exact_al_range = None
latentheat_melting_al_range = None
melting_ice_equilibrium_sa_ct_ratio_al_range = None
metling_ice_equilibrium_sa_ct_ratio_al_range = None
mid_lat_array_range = None
mid_lon_array_range = None
n2_array_range = None
p_al_range = None
pot_rho_t_exact_al_range = None
pressure_coefficient_ice_al_range = None
pressure_freezing_al_range = None
pt_ct_al_range = None
pt_ct_ct_al_range = None
pt_ice_al_range = None
pt_sa_al_range = None
pt_sa_ct_al_range = None
pt_sa_sa_al_range = None
pt0_al_range = None
rho_ct_al_range = None
rho_ct_ct_al_range = None
rho_ct_p_al_range = None
rho_h_al_range = None
rho_h_h_al_range = None
rho_ice_al_range = None
rho_p_al_range = None
rho_sa_ct_al_range = None
rho_sa_h_al_range = None
rho_sa_p_al_range = None
rho_t_exact_al_range = None
rsubrho_array_range = None
sa_baltic_al_range = None
saar_al_range = None
sigma0_al_range = None
sigma1_al_range = None
sigma2_al_range = None
sigma3_al_range = None
sigma4_al_range = None
sound_speed_al_range = None
sound_speed_ice_al_range = None
sound_speed_t_exact_al_range = None
sp_baltic_al_range = None
spec_vol_anom_al_range = None
spec_vol_ice_al_range = None
spec_vol_t_exact_al_range = None
spiciness0_al_range = None
spiciness1_al_range = None
spiciness2_al_range = None
thermobaric_al_range = None
tu_array_range = None
v_ct_al_range = None
v_ct_ct_al_range = None
v_ct_p_al_range = None
v_h_al_range = None
v_h_h_al_range = None
v_sa_al_range = None
v_sa_ct_al_range = None
v_sa_h_al_range = None
v_sa_p_al_range = None
v_sa_sa_al_range = None
v_sa_sa_wrt_h_al_range = None
v_sa_wrt_h_al_range = None
velocity_array_range = None
w_ih_al_range = None
w_seaice_al_range = None
yi_array_range = None
z_al_range = None


# Utility functions


def check_value_exists(value):
    """Check a single value is not nan/none"""

    if not np.isscalar(value):
        # Should not be checking against list/array/etc
        raise Exception("Value is not scalar")

    if value == None:
        raise Exception("Value is None")
    elif np.isnan(value):
        raise Exception("Value is np.nan")


def check_value_in_range(value, range):
    """Check value is in the specified range"""

    # Skip if value is nan
    if np.isnan(value):
        return

    if value is None:
        print(f"VALUE AHH: {value}")
        raise Exception("You did not include a value...")

    if not range:
        return  # Lots of empty ranges
        # raise Exception("You did not include a range...")

    minimum, maximum = range

    if minimum is not None:
        assert minimum <= value

    if maximum:
        assert value <= maximum


# Tests


class TestGSW:
    class TestConversionFunctions:
        # @given(
        #     absolute_salinity_st,
        #     specific_enthalpy_st,
        #     sea_pressure_st,
        # )
        # def test_CT_from_enthalpy(self, SA, h, p):
        #     conservative_temperature_al = CT_from_enthalpy(SA, h, p)
        #
        #     check_value_exists(conservative_temperature_al)
        #     check_value_in_range(conservative_temperature_al, conservative_temperature_al_range)

        # @given(
        #     absolute_salinity_st,
        #     specific_entropy_st,
        # )
        # def test_CT_from_entropy(self, SA, entropy):
        #     conservative_temperature_al = CT_from_entropy(SA, entropy)
        #
        #     check_value_exists(conservative_temperature_al)
        #     check_value_in_range(conservative_temperature_al, conservative_temperature_al_range)

        @given(
            absolute_salinity_st,
            potential_temperature_st,
        )
        def test_CT_from_pt(self, SA, pt):
            conservative_temperature_al = CT_from_pt(SA, pt)

            check_value_exists(conservative_temperature_al)
            check_value_in_range(
                conservative_temperature_al, conservative_temperature_al_range
            )

        @given(
            in_situ_seawater_density_rho_st,
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_CT_from_rho(self, rho, SA, p):
            conservative_temperature_al, ct_multiple_al = CT_from_rho(rho, SA, p)

            check_value_exists(conservative_temperature_al)
            check_value_in_range(
                conservative_temperature_al, conservative_temperature_al_range
            )
            check_value_exists(ct_multiple_al)
            check_value_in_range(ct_multiple_al, ct_multiple_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_CT_from_t(self, SA, t, p):
            conservative_temperature_al = CT_from_t(SA, t, p)

            check_value_exists(conservative_temperature_al)
            check_value_in_range(
                conservative_temperature_al, conservative_temperature_al_range
            )

        @given(
            practical_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_C_from_SP(self, SP, t, p):
            c_al = C_from_SP(SP, t, p)

            check_value_exists(c_al)
            check_value_in_range(c_al, c_al_range)

        @given(
            practical_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_SA_from_SP(self, SP, p, lon, lat):
            absolute_salinity_al = SA_from_SP(SP, p, lon, lat)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            preformed_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_SA_from_Sstar(self, Sstar, p, lon, lat):
            absolute_salinity_al = SA_from_Sstar(Sstar, p, lon, lat)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            in_situ_seawater_density_rho_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_SA_from_rho(self, rho, CT, p):
            absolute_salinity_al = SA_from_rho(rho, CT, p)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            conductivity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_SP_from_C(self, C, t, p):
            practical_salinity_al = SP_from_C(C, t, p)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_SP_from_SA(self, SA, p, lon, lat):
            practical_salinity_al = SP_from_SA(SA, p, lon, lat)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            knudsen_salinity_st,
        )
        def test_SP_from_SK(self, SK):
            practical_salinity_al = SP_from_SK(SK)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            reference_salinity_st,
        )
        def test_SP_from_SR(self, SR):
            practical_salinity_al = SP_from_SR(SR)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            preformed_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_SP_from_Sstar(self, Sstar, p, lon, lat):
            practical_salinity_al = SP_from_Sstar(Sstar, p, lon, lat)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            practical_salinity_st,
        )
        def test_SR_from_SP(self, SP):
            practical_salinity_al = SR_from_SP(SP)

            check_value_exists(practical_salinity_al)
            check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_Sstar_from_SA(self, SA, p, lon, lat):
            sstar_al = Sstar_from_SA(SA, p, lon, lat)

            check_value_exists(sstar_al)
            check_value_in_range(sstar_al, sstar_al_range)

        @given(
            practical_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_Sstar_from_SP(self, SP, p, lon, lat):
            sstar_al = Sstar_from_SP(SP, p, lon, lat)

            check_value_exists(sstar_al)
            check_value_in_range(sstar_al, sstar_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_adiabatic_lapse_rate_from_CT(self, SA, CT, p):
            adiabatic_lapse_rate_al = adiabatic_lapse_rate_from_CT(SA, CT, p)

            check_value_exists(adiabatic_lapse_rate_al)
            check_value_in_range(adiabatic_lapse_rate_al, adiabatic_lapse_rate_al_range)

        @given(
            practical_salinity_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_deltaSA_from_SP(self, SP, p, lon, lat):
            delta_sa_al = deltaSA_from_SP(SP, p, lon, lat)

            check_value_exists(delta_sa_al)
            check_value_in_range(delta_sa_al, delta_sa_al_range)

        @given(
            absolute_salinity_st,
            potential_temperature_st,
        )
        def test_entropy_from_pt(self, SA, pt):
            specific_entropy_al = entropy_from_pt(SA, pt)

            check_value_exists(specific_entropy_al)
            check_value_in_range(specific_entropy_al, specific_entropy_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_entropy_from_t(self, SA, t, p):
            specific_entropy_al = entropy_from_t(SA, t, p)

            check_value_exists(specific_entropy_al)
            check_value_in_range(specific_entropy_al, specific_entropy_al_range)

        # @given(
        #     depth_st,
        #     latitude_st,
        #     dynamic_height_anomoly_st,
        #     geopotential_at_zero_sea_pressure_st,
        # )
        # def test_p_from_z(self, z, lat, geo_strf_dyn_height, sea_surface_geopotential):
        #     p_al = p_from_z(z, lat, geo_strf_dyn_height, sea_surface_geopotential)
        #
        #     check_value_exists(p_al)
        #     check_value_in_range(p_al, p_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_pt0_from_t(self, SA, t, p):
            pt0_al = pt0_from_t(SA, t, p)

            check_value_exists(pt0_al)
            check_value_in_range(pt0_al, pt0_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_pt_from_CT(self, SA, CT):
            potential_temperature_al = pt_from_CT(SA, CT)

            check_value_exists(potential_temperature_al)
            check_value_in_range(
                potential_temperature_al, potential_temperature_al_range
            )

        # @given(
        #     absolute_salinity_st,
        #     specific_entropy_st,
        # )
        # def test_pt_from_entropy(self, SA, entropy):
        #     potential_temperature_al = pt_from_entropy(SA, entropy)
        #
        #     check_value_exists(potential_temperature_al)
        #     check_value_in_range(potential_temperature_al, potential_temperature_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
            reference_pressure_st,
        )
        def test_pt_from_t(self, SA, t, p, p_ref):
            potential_temperature_al = pt_from_t(SA, t, p, p_ref)

            check_value_exists(potential_temperature_al)
            check_value_in_range(
                potential_temperature_al, potential_temperature_al_range
            )

        @given(
            temperature_t68_st,
        )
        def test_t90_from_t68(self, t68):
            its_90 = t90_from_t68(t68)

            check_value_exists(its_90)
            check_value_in_range(its_90, its_90_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_t_from_CT(self, SA, CT, p):
            t_al = t_from_CT(SA, CT, p)

            check_value_exists(t_al)
            check_value_in_range(t_al, t_al_range)

        # @given(
        #     sea_pressure_st,
        #     latitude_st,
        #     dynamic_height_anomoly_st,
        #     geopotential_at_zero_sea_pressure_st,
        # )
        # def test_z_from_p(self, p, lat, geo_strf_dyn_height, sea_surface_geopotential):
        #     z_al = z_from_p(p, lat, geo_strf_dyn_height, sea_surface_geopotential)
        #
        #     check_value_exists(z_al)
        #     check_value_in_range(z_al, z_al_range)

    class TestDensity:
        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_alpha(self, SA, CT, p):
            alpha_al = alpha(SA, CT, p)

            check_value_exists(alpha_al)
            check_value_in_range(alpha_al, alpha_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_alpha_on_beta(self, SA, CT, p):
            alpha_on_beta_al = alpha_on_beta(SA, CT, p)

            check_value_exists(alpha_on_beta_al)
            check_value_in_range(alpha_on_beta_al, alpha_on_beta_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_beta(self, SA, CT, p):
            beta_al = beta(SA, CT, p)

            check_value_exists(beta_al)
            check_value_in_range(beta_al, beta_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_kappa(self, SA, CT, p):
            kappa_al = kappa(SA, CT, p)

            check_value_exists(kappa_al)
            check_value_in_range(kappa_al, kappa_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho(self, SA, CT, p):
            rho_al = rho(SA, CT, p)

            check_value_exists(rho_al)
            check_value_in_range(rho_al, rho_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho_alpha_beta(self, SA, CT, p):
            rho_al, alpha_al, beta_al = rho_alpha_beta(SA, CT, p)

            check_value_exists(rho_al)
            check_value_in_range(rho_al, rho_al_range)
            check_value_exists(alpha_al)
            check_value_in_range(alpha_al, alpha_al_range)
            check_value_exists(beta_al)
            check_value_in_range(beta_al, beta_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_rho_t_exact(self, SA, t, p):
            rho_t_exact_al = rho_t_exact(SA, t, p)

            check_value_exists(rho_t_exact_al)
            check_value_in_range(rho_t_exact_al, rho_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_sigma0(self, SA, CT):
            sigma0_al = sigma0(SA, CT)

            check_value_exists(sigma0_al)
            check_value_in_range(sigma0_al, sigma0_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_sigma1(self, SA, CT):
            sigma1_al = sigma1(SA, CT)

            check_value_exists(sigma1_al)
            check_value_in_range(sigma1_al, sigma1_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_sigma2(self, SA, CT):
            sigma2_al = sigma2(SA, CT)

            check_value_exists(sigma2_al)
            check_value_in_range(sigma2_al, sigma2_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_sigma3(self, SA, CT):
            sigma3_al = sigma3(SA, CT)

            check_value_exists(sigma3_al)
            check_value_in_range(sigma3_al, sigma3_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_sigma4(self, SA, CT):
            sigma4_al = sigma4(SA, CT)

            check_value_exists(sigma4_al)
            check_value_in_range(sigma4_al, sigma4_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_sound_speed(self, SA, CT, p):
            sound_speed_al = sound_speed(SA, CT, p)

            check_value_exists(sound_speed_al)
            check_value_in_range(sound_speed_al, sound_speed_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol(self, SA, CT, p):
            spec_vol_al = specvol(SA, CT, p)

            check_value_exists(spec_vol_al)
            check_value_in_range(spec_vol_al, spec_vol_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_alpha_beta(self, SA, CT, p):
            spec_vol_al, alpha_al, beta_al = specvol_alpha_beta(SA, CT, p)

            check_value_exists(spec_vol_al)
            check_value_in_range(spec_vol_al, spec_vol_al_range)
            check_value_exists(alpha_al)
            check_value_in_range(alpha_al, alpha_al_range)
            check_value_exists(beta_al)
            check_value_in_range(beta_al, beta_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_anom_standard(self, SA, CT, p):
            spec_vol_anom_al = specvol_anom_standard(SA, CT, p)

            check_value_exists(spec_vol_anom_al)
            check_value_in_range(spec_vol_anom_al, spec_vol_anom_al_range)

    class TestEnergy:
        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy(self, SA, CT, p):
            enthalpy_al = enthalpy(SA, CT, p)

            check_value_exists(enthalpy_al)
            check_value_in_range(enthalpy_al, enthalpy_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            upper_sea_pressure_st,
            lower_sea_pressure_st,
        )
        def test_enthalpy_diff(self, SA, CT, p_shallow, p_deep):
            enthalpy_diff_al = enthalpy_diff(SA, CT, p_shallow, p_deep)

            check_value_exists(enthalpy_diff_al)
            check_value_in_range(enthalpy_diff_al, enthalpy_diff_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_internal_energy(self, SA, CT, p):
            internal_energy_al = internal_energy(SA, CT, p)

            check_value_exists(internal_energy_al)
            check_value_in_range(internal_energy_al, internal_energy_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_latentheat_evap_CT(self, SA, CT):
            latentheat_evap_al = latentheat_evap_CT(SA, CT)

            check_value_exists(latentheat_evap_al)
            check_value_in_range(latentheat_evap_al, latentheat_evap_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
        )
        def test_latentheat_evap_t(self, SA, t):
            latentheat_evap_al = latentheat_evap_t(SA, t)

            check_value_exists(latentheat_evap_al)
            check_value_in_range(latentheat_evap_al, latentheat_evap_al_range)

    # class TestStability:
    #     # @given(
    #     #     absolute_salinity_st,
    #     #     conservative_temperature_st,
    #     #     sea_pressure_st,
    #     #     reference_pressure_st,
    #     #     axis_st,
    #     # )
    #     # def test_IPV_vs_fNsquared_ratio(self, SA, CT, p, p_ref, axis):
    #     #     ipv_vs_fn_squared, pressure_midpoint_array = IPV_vs_fNsquared_ratio(SA, CT, p, p_ref, axis)
    #     #
    #     #     check_value_exists(ipv_vs_fn_squared)
    #     #     check_value_in_range(ipv_vs_fn_squared, ipv_vs_fn_squared_range)
    #     #     check_value_exists(pressure_midpoint_array)
    #     #     check_value_in_range(pressure_midpoint_array, pressure_midpoint_array_range)

    #     # @given(
    #     #     absolute_salinity_st,
    #     #     conservative_temperature_st,
    #     #     sea_pressure_st,
    #     #     latitude_st,
    #     #     axis_st,
    #     # )
    #     # def test_Nsquared(self, SA, CT, p, lat, axis):
    #     #     n2_array, pressure_midpoint_array = Nsquared(SA, CT, p, lat, axis)
    #     #
    #     #     check_value_exists(n2_array)
    #     #     check_value_in_range(n2_array, n2_array_range)
    #     #     check_value_exists(pressure_midpoint_array)
    #     #     check_value_in_range(pressure_midpoint_array, pressure_midpoint_array_range)

    #     # @given(
    #     #     absolute_salinity_st,
    #     #     conservative_temperature_st,
    #     #     sea_pressure_st,
    #     #     axis_st,
    #     # )
    #     # def test_Turner_Rsubrho(self, SA, CT, p, axis):
    #     #     tu_array, rsubrho_array, pressure_midpoint_array = Turner_Rsubrho(SA, CT, p, axis)
    #     #
    #     #     check_value_exists(tu_array)
    #     #     check_value_in_range(tu_array, tu_array_range)
    #     #     check_value_exists(rsubrho_array)
    #     #     check_value_in_range(rsubrho_array, rsubrho_array_range)
    #     #     check_value_exists(pressure_midpoint_array)
    #     #     check_value_in_range(pressure_midpoint_array, pressure_midpoint_array_range)

    class TestGeostrophy:
        # @given(
        #     longitude_st,
        #     latitude_st,
        #     sea_pressure_st,
        #     axis_st,
        # )
        # def test_distance(self, lon, lat, p, axis):
        #     distance_between_points = distance(lon, lat, p, axis)
        #
        #     check_value_exists(distance_between_points)
        #     check_value_in_range(distance_between_points, distance_between_points_range)

        @given(
            latitude_st,
        )
        def test_f(self, lat):
            result = f(lat)

            check_value_exists(result)

        # @given(
        #     absolute_salinity_st,
        #     conservative_temperature_st,
        #     sea_pressure_st,
        #     reference_pressure_st,
        #     axis_st,
        #     max_dp_st,
        #     interpolation_algorithm_st,
        # )
        # def test_geo_strf_dyn_height(self, SA, CT, p, p_ref, axis, max_dp, interp_method):
        #     dynamic_height_array = geo_strf_dyn_height(SA, CT, p, p_ref, axis, max_dp, interp_method)
        #
        #     check_value_exists(dynamic_height_array)
        #     check_value_in_range(dynamic_height_array, dynamic_height_array_range)

        # @given(
        #     geostrophic_stream_function_st,
        #     longitude_st,
        #     latitude_st,
        #     sea_pressure_st,
        #     axis_st,
        # )
        # def test_geostrophic_velocity(self, geo_strf, lon, lat, p, axis):
        #     velocity_array, mid_lon_array, mid_lat_array = geostrophic_velocity(geo_strf, lon, lat, p, axis)
        #
        #     check_value_exists(velocity_array)
        #     check_value_in_range(velocity_array, velocity_array_range)
        #     check_value_exists(mid_lon_array)
        #     check_value_in_range(mid_lon_array, mid_lon_array_range)
        #     check_value_exists(mid_lat_array)
        #     check_value_in_range(mid_lat_array, mid_lat_array_range)

    class TestIce:
        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_Helmholtz_energy_ice(self, t, p):
            helmholtz_energy_ice_al = Helmholtz_energy_ice(t, p)

            check_value_exists(helmholtz_energy_ice_al)
            check_value_in_range(helmholtz_energy_ice_al, helmholtz_energy_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_adiabatic_lapse_rate_ice(self, t, p):
            adiabatic_lapse_rate_ice_al = adiabatic_lapse_rate_ice(t, p)

            check_value_exists(adiabatic_lapse_rate_ice_al)
            check_value_in_range(
                adiabatic_lapse_rate_ice_al, adiabatic_lapse_rate_ice_al_range
            )

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_alpha_wrt_t_ice(self, t, p):
            alpha_wrt_t_ice_al = alpha_wrt_t_ice(t, p)

            check_value_exists(alpha_wrt_t_ice_al)
            check_value_in_range(alpha_wrt_t_ice_al, alpha_wrt_t_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_chem_potential_water_ice(self, t, p):
            chem_potential_water_ice_al = chem_potential_water_ice(t, p)

            check_value_exists(chem_potential_water_ice_al)
            check_value_in_range(
                chem_potential_water_ice_al, chem_potential_water_ice_al_range
            )

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_cp_ice(self, t, p):
            cp_ice_al = cp_ice(t, p)

            check_value_exists(cp_ice_al)
            check_value_in_range(cp_ice_al, cp_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_ice(self, t, p):
            enthalpy_ice_al = enthalpy_ice(t, p)

            check_value_exists(enthalpy_ice_al)
            check_value_in_range(enthalpy_ice_al, enthalpy_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_entropy_ice(self, t, p):
            ice_entropy_al = entropy_ice(t, p)

            check_value_exists(ice_entropy_al)
            check_value_in_range(ice_entropy_al, ice_entropy_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            in_situ_temperature_of_ice_st,
        )
        def test_ice_fraction_to_freeze_seawater(self, SA, CT, p, t_Ih):
            sa_freeze_al, ct_freeze_al, w_ih_al = ice_fraction_to_freeze_seawater(
                SA, CT, p, t_Ih
            )

            check_value_exists(sa_freeze_al)
            check_value_in_range(sa_freeze_al, sa_freeze_al_range)
            check_value_exists(ct_freeze_al)
            check_value_in_range(ct_freeze_al, ct_freeze_al_range)
            check_value_exists(w_ih_al)
            check_value_in_range(w_ih_al, w_ih_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_internal_energy_ice(self, t, p):
            internal_energy_ice_al = internal_energy_ice(t, p)

            check_value_exists(internal_energy_ice_al)
            check_value_in_range(internal_energy_ice_al, internal_energy_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_kappa_const_t_ice(self, t, p):
            kappa_const_t_ice_al = kappa_const_t_ice(t, p)

            check_value_exists(kappa_const_t_ice_al)
            check_value_in_range(kappa_const_t_ice_al, kappa_const_t_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_kappa_ice(self, t, p):
            kappa_ice_al = kappa_ice(t, p)

            check_value_exists(kappa_ice_al)
            check_value_in_range(kappa_ice_al, kappa_ice_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            in_situ_temperature_of_ice_st,
        )
        def test_melting_ice_SA_CT_ratio(self, SA, CT, p, t_Ih):
            melting_ice_sa_ct_ratio_al = melting_ice_SA_CT_ratio(SA, CT, p, t_Ih)

            check_value_exists(melting_ice_sa_ct_ratio_al)
            check_value_in_range(
                melting_ice_sa_ct_ratio_al, melting_ice_sa_ct_ratio_al_range
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            in_situ_temperature_of_ice_st,
        )
        def test_melting_ice_SA_CT_ratio_poly(self, SA, CT, p, t_Ih):
            melting_ice_sa_ct_ratio_al = melting_ice_SA_CT_ratio_poly(SA, CT, p, t_Ih)

            check_value_exists(melting_ice_sa_ct_ratio_al)
            check_value_in_range(
                melting_ice_sa_ct_ratio_al, melting_ice_sa_ct_ratio_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_melting_ice_equilibrium_SA_CT_ratio(self, SA, p):
            metling_ice_equilibrium_sa_ct_ratio_al = (
                melting_ice_equilibrium_SA_CT_ratio(SA, p)
            )

            check_value_exists(metling_ice_equilibrium_sa_ct_ratio_al)
            check_value_in_range(
                metling_ice_equilibrium_sa_ct_ratio_al,
                metling_ice_equilibrium_sa_ct_ratio_al_range,
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_melting_ice_equilibrium_SA_CT_ratio_poly(self, SA, p):
            melting_ice_equilibrium_sa_ct_ratio_al = (
                melting_ice_equilibrium_SA_CT_ratio_poly(SA, p)
            )

            check_value_exists(melting_ice_equilibrium_sa_ct_ratio_al)
            check_value_in_range(
                melting_ice_equilibrium_sa_ct_ratio_al,
                melting_ice_equilibrium_sa_ct_ratio_al_range,
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            mass_fraction_of_ice_st,
            in_situ_temperature_of_ice_st,
        )
        def test_melting_ice_into_seawater(self, SA, CT, p, w_Ih, t_Ih):
            (
                absolute_salinity_final_al,
                conservative_temperature_final_al,
                w_ih_final_al,
            ) = melting_ice_into_seawater(SA, CT, p, w_Ih, t_Ih)

            check_value_exists(absolute_salinity_final_al)
            check_value_in_range(
                absolute_salinity_final_al, absolute_salinity_final_al_range
            )
            check_value_exists(conservative_temperature_final_al)
            check_value_in_range(
                conservative_temperature_final_al,
                conservative_temperature_final_al_range,
            )
            check_value_exists(w_ih_final_al)
            check_value_in_range(w_ih_final_al, w_ih_final_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            absolute_salinity_of_sea_ice_st,
            in_situ_temperature_of_sea_ice_at_pressure_st,
        )
        def test_melting_seaice_SA_CT_ratio(self, SA, CT, p, SA_seaice, t_seaice):
            melting_seaice_sa_ct_ratio_al = melting_seaice_SA_CT_ratio(
                SA, CT, p, SA_seaice, t_seaice
            )

            check_value_exists(melting_seaice_sa_ct_ratio_al)
            check_value_in_range(
                melting_seaice_sa_ct_ratio_al, melting_seaice_sa_ct_ratio_al_range
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            absolute_salinity_of_sea_ice_st,
            in_situ_temperature_of_sea_ice_at_pressure_st,
        )
        def test_melting_seaice_SA_CT_ratio_poly(self, SA, CT, p, SA_seaice, t_seaice):
            melting_seaice_sa_ct_ratio_al = melting_seaice_SA_CT_ratio_poly(
                SA, CT, p, SA_seaice, t_seaice
            )

            check_value_exists(melting_seaice_sa_ct_ratio_al)
            check_value_in_range(
                melting_seaice_sa_ct_ratio_al, melting_seaice_sa_ct_ratio_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_melting_seaice_equilibrium_SA_CT_ratio(self, SA, p):
            melting_seaice_equilibrium_sa_ct_ratio_al = (
                melting_seaice_equilibrium_SA_CT_ratio(SA, p)
            )

            check_value_exists(melting_seaice_equilibrium_sa_ct_ratio_al)
            check_value_in_range(
                melting_seaice_equilibrium_sa_ct_ratio_al,
                melting_seaice_equilibrium_sa_ct_ratio_al_range,
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_melting_seaice_equilibrium_SA_CT_ratio_poly(self, SA, p):
            melting_seaice_equilibrium_sa_ct_ratio_al = (
                melting_seaice_equilibrium_SA_CT_ratio_poly(SA, p)
            )

            check_value_exists(melting_seaice_equilibrium_sa_ct_ratio_al)
            check_value_in_range(
                melting_seaice_equilibrium_sa_ct_ratio_al,
                melting_seaice_equilibrium_sa_ct_ratio_al_range,
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            mass_fraction_of_sea_ice_st,
            absolute_salinity_of_sea_ice_st,
            in_situ_temperature_of_sea_ice_at_pressure_st,
        )
        def test_melting_seaice_into_seawater(
            self, SA, CT, p, w_seaice, SA_seaice, t_seaice
        ):
            (
                absolute_salinity_final_al,
                conservative_temperature_final_al,
            ) = melting_seaice_into_seawater(SA, CT, p, w_seaice, SA_seaice, t_seaice)

            check_value_exists(absolute_salinity_final_al)
            check_value_in_range(
                absolute_salinity_final_al, absolute_salinity_final_al_range
            )
            check_value_exists(conservative_temperature_final_al)
            check_value_in_range(
                conservative_temperature_final_al,
                conservative_temperature_final_al_range,
            )

        @given(
            potential_temperature_of_ice_with_reference_to_zero_pressure_st,
        )
        def test_pot_enthalpy_from_pt_ice(self, pt0_ice):
            pot_enthalpy_ice_al = pot_enthalpy_from_pt_ice(pt0_ice)

            check_value_exists(pot_enthalpy_ice_al)
            check_value_in_range(pot_enthalpy_ice_al, pot_enthalpy_ice_al_range)

        @given(
            potential_temperature_of_ice_with_reference_to_zero_pressure_st,
        )
        def test_pot_enthalpy_from_pt_ice_poly(self, pt0_ice):
            pot_enthalpy_ice_al = pot_enthalpy_from_pt_ice_poly(pt0_ice)

            check_value_exists(pot_enthalpy_ice_al)
            check_value_in_range(pot_enthalpy_ice_al, pot_enthalpy_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_pressure_coefficient_ice(self, t, p):
            pressure_coefficient_ice_al = pressure_coefficient_ice(t, p)

            check_value_exists(pressure_coefficient_ice_al)
            check_value_in_range(
                pressure_coefficient_ice_al, pressure_coefficient_ice_al_range
            )

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_pt0_from_t_ice(self, t, p):
            potential_temperature_ice_al = pt0_from_t_ice(t, p)

            check_value_exists(potential_temperature_ice_al)
            check_value_in_range(
                potential_temperature_ice_al, potential_temperature_ice_al_range
            )

        # @given(
        #     potential_enthalpy_of_ice_st,
        # )
        # def test_pt_from_pot_enthalpy_ice(self, pot_enthalpy_ice):
        #     potential_temperature_ice_al = pt_from_pot_enthalpy_ice(pot_enthalpy_ice)
        #
        #     check_value_exists(potential_temperature_ice_al)
        #     check_value_in_range(potential_temperature_ice_al, potential_temperature_ice_al_range)

        # @given(
        #     potential_enthalpy_of_ice_st,
        # )
        # def test_pt_from_pot_enthalpy_ice_poly(self, pot_enthalpy_ice):
        #     potential_temperature_ice_al = pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice)
        #
        #     check_value_exists(potential_temperature_ice_al)
        #     check_value_in_range(potential_temperature_ice_al, potential_temperature_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
            reference_pressure_st,
        )
        def test_pt_from_t_ice(self, t, p, p_ref):
            pt_ice_al = pt_from_t_ice(t, p, p_ref)

            check_value_exists(pt_ice_al)
            check_value_in_range(pt_ice_al, pt_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_rho_ice(self, t, p):
            rho_ice_al = rho_ice(t, p)

            check_value_exists(rho_ice_al)
            check_value_in_range(rho_ice_al, rho_ice_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            absolute_salinity_of_sea_ice_st,
            in_situ_temperature_of_sea_ice_at_pressure_st,
        )
        def test_seaice_fraction_to_freeze_seawater(
            self, SA, CT, p, SA_seaice, t_seaice
        ):
            (
                sa_freeze_al,
                ct_freeze_al,
                w_seaice_al,
            ) = seaice_fraction_to_freeze_seawater(SA, CT, p, SA_seaice, t_seaice)

            check_value_exists(sa_freeze_al)
            check_value_in_range(sa_freeze_al, sa_freeze_al_range)
            check_value_exists(ct_freeze_al)
            check_value_in_range(ct_freeze_al, ct_freeze_al_range)
            check_value_exists(w_seaice_al)
            check_value_in_range(w_seaice_al, w_seaice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_sound_speed_ice(self, t, p):
            sound_speed_ice_al = sound_speed_ice(t, p)

            check_value_exists(sound_speed_ice_al)
            check_value_in_range(sound_speed_ice_al, sound_speed_ice_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_ice(self, t, p):
            spec_vol_ice_al = specvol_ice(t, p)

            check_value_exists(spec_vol_ice_al)
            check_value_in_range(spec_vol_ice_al, spec_vol_ice_al_range)

        @given(
            potential_temperature_of_ice_with_reference_to_zero_pressure_st,
            sea_pressure_st,
        )
        def test_t_from_pt0_ice(self, pt0_ice, p):
            t_al = t_from_pt0_ice(pt0_ice, p)

            check_value_exists(t_al)
            check_value_in_range(t_al, t_al_range)

    class TestFreezing:
        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_CT_freezing(self, SA, p, saturation_fraction):
            ct_freezing_al = CT_freezing(SA, p, saturation_fraction)

            check_value_exists(ct_freezing_al)
            check_value_in_range(ct_freezing_al, ct_freezing_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_CT_freezing_first_derivatives(self, SA, p, saturation_fraction):
            ct_freezing_sa_al, ct_freezing_p_al = CT_freezing_first_derivatives(
                SA, p, saturation_fraction
            )

            check_value_exists(ct_freezing_sa_al)
            check_value_in_range(ct_freezing_sa_al, ct_freezing_sa_al_range)
            check_value_exists(ct_freezing_p_al)
            check_value_in_range(ct_freezing_p_al, ct_freezing_p_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_CT_freezing_first_derivatives_poly(self, SA, p, saturation_fraction):
            ct_freezing_sa_al, ct_freezing_p_al = CT_freezing_first_derivatives_poly(
                SA, p, saturation_fraction
            )

            check_value_exists(ct_freezing_sa_al)
            check_value_in_range(ct_freezing_sa_al, ct_freezing_sa_al_range)
            check_value_exists(ct_freezing_p_al)
            check_value_in_range(ct_freezing_p_al, ct_freezing_p_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_CT_freezing_poly(self, SA, p, saturation_fraction):
            ct_freezing_al = CT_freezing_poly(SA, p, saturation_fraction)

            check_value_exists(ct_freezing_al)
            check_value_in_range(ct_freezing_al, ct_freezing_al_range)

        @given(
            conservative_temperature_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_SA_freezing_from_CT(self, CT, p, saturation_fraction):
            absolute_salinity_al = SA_freezing_from_CT(CT, p, saturation_fraction)
            target(-absolute_salinity_al)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            conservative_temperature_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_SA_freezing_from_CT_poly(self, CT, p, saturation_fraction):
            absolute_salinity_al = SA_freezing_from_CT_poly(CT, p, saturation_fraction)
            target(-absolute_salinity_al)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_SA_freezing_from_t(self, t, p, saturation_fraction):
            absolute_salinity_al = SA_freezing_from_t(t, p, saturation_fraction)
            target(-absolute_salinity_al)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_SA_freezing_from_t_poly(self, t, p, saturation_fraction):
            absolute_salinity_al = SA_freezing_from_t_poly(t, p, saturation_fraction)

            check_value_exists(absolute_salinity_al)
            check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_pot_enthalpy_ice_freezing(self, SA, p):
            pot_enthalpy_ice_freezing_al = pot_enthalpy_ice_freezing(SA, p)

            check_value_exists(pot_enthalpy_ice_freezing_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_al, pot_enthalpy_ice_freezing_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_pot_enthalpy_ice_freezing_first_derivatives(self, SA, p):
            (
                pot_enthalpy_ice_freezing_sa_al,
                pot_enthalpy_ice_freezing_p_al,
            ) = pot_enthalpy_ice_freezing_first_derivatives(SA, p)

            check_value_exists(pot_enthalpy_ice_freezing_sa_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_sa_al, pot_enthalpy_ice_freezing_sa_al_range
            )
            check_value_exists(pot_enthalpy_ice_freezing_p_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_p_al, pot_enthalpy_ice_freezing_p_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_pot_enthalpy_ice_freezing_first_derivatives_poly(self, SA, p):
            (
                pot_enthalpy_ice_freezing_sa_al,
                pot_enthalpy_ice_freezing_p_al,
            ) = pot_enthalpy_ice_freezing_first_derivatives_poly(SA, p)

            check_value_exists(pot_enthalpy_ice_freezing_sa_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_sa_al, pot_enthalpy_ice_freezing_sa_al_range
            )
            check_value_exists(pot_enthalpy_ice_freezing_p_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_p_al, pot_enthalpy_ice_freezing_p_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_pot_enthalpy_ice_freezing_poly(self, SA, p):
            pot_enthalpy_ice_freezing_al = pot_enthalpy_ice_freezing_poly(SA, p)

            check_value_exists(pot_enthalpy_ice_freezing_al)
            check_value_in_range(
                pot_enthalpy_ice_freezing_al, pot_enthalpy_ice_freezing_al_range
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            saturation_fraction_st,
        )
        def test_pressure_freezing_CT(self, SA, CT, saturation_fraction):
            pressure_freezing_al = pressure_freezing_CT(SA, CT, saturation_fraction)

            check_value_exists(pressure_freezing_al)
            check_value_in_range(pressure_freezing_al, pressure_freezing_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_t_freezing(self, SA, p, saturation_fraction):
            t_freezing_al = t_freezing(SA, p, saturation_fraction)

            check_value_exists(t_freezing_al)
            check_value_in_range(t_freezing_al, t_freezing_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_t_freezing_first_derivatives(self, SA, p, saturation_fraction):
            t_freezing_sa_al, t_freezing_p_al = t_freezing_first_derivatives(
                SA, p, saturation_fraction
            )

            check_value_exists(t_freezing_sa_al)
            check_value_in_range(t_freezing_sa_al, t_freezing_sa_al_range)
            check_value_exists(t_freezing_p_al)
            check_value_in_range(t_freezing_p_al, t_freezing_p_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_t_freezing_first_derivatives_poly(self, SA, p, saturation_fraction):
            t_freezing_sa_al, t_freezing_p_al = t_freezing_first_derivatives_poly(
                SA, p, saturation_fraction
            )

            check_value_exists(t_freezing_sa_al)
            check_value_in_range(t_freezing_sa_al, t_freezing_sa_al_range)
            check_value_exists(t_freezing_p_al)
            check_value_in_range(t_freezing_p_al, t_freezing_p_al_range)

    class TestMisc:
        @given(
            absolute_salinity_st,
            potential_temperature_st,
        )
        def test_CT_first_derivatives(self, SA, pt):
            ct_sa_al, ct_pt_al = CT_first_derivatives(SA, pt)

            check_value_exists(ct_sa_al)
            check_value_in_range(ct_sa_al, ct_sa_al_range)
            check_value_exists(ct_pt_al)
            check_value_in_range(ct_pt_al, ct_pt_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_CT_first_derivatives_wrt_t_exact(self, SA, t, p):
            (
                ct_sa_wrt_t_al,
                ct_t_wrt_t_al,
                ct_p_wrt_t_al,
            ) = CT_first_derivatives_wrt_t_exact(SA, t, p)

            check_value_exists(ct_sa_wrt_t_al)
            check_value_in_range(ct_sa_wrt_t_al, ct_sa_wrt_t_al_range)
            check_value_exists(ct_t_wrt_t_al)
            check_value_in_range(ct_t_wrt_t_al, ct_t_wrt_t_al_range)
            check_value_exists(ct_p_wrt_t_al)
            check_value_in_range(ct_p_wrt_t_al, ct_p_wrt_t_al_range)

        # @given(
        #     absolute_salinity_st,
        #     specific_enthalpy_st,
        #     sea_pressure_st,
        # )
        # def test_CT_from_enthalpy_exact(self, SA, h, p):
        #     conservative_temperature_al = CT_from_enthalpy_exact(SA, h, p)
        #
        #     check_value_exists(conservative_temperature_al)
        #     check_value_in_range(conservative_temperature_al, conservative_temperature_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_CT_maxdensity(self, SA, p):
            ct_maxdensity_al = CT_maxdensity(SA, p)

            check_value_exists(ct_maxdensity_al)
            check_value_in_range(ct_maxdensity_al, ct_maxdensity_al_range)

        # @given(
        #     absolute_salinity_st,
        #     potential_temperature_st,
        # )
        # def test_CT_second_derivatives(self, SA, pt):
        #     (
        #         ct_sa_sa_al,
        #         ct_sa_pt_al,
        #         p_ref_al,
        #         ct_pt_pt_al,
        #         p_ref_al,
        #     ) = CT_second_derivatives(SA, pt)
        #
        #     check_value_exists(ct_sa_sa_al)
        #     check_value_in_range(ct_sa_sa_al, ct_sa_sa_al_range)
        #     check_value_exists(ct_sa_pt_al)
        #     check_value_in_range(ct_sa_pt_al, ct_sa_pt_al_range)
        #     check_value_exists(p_ref_al)
        #     check_value_in_range(p_ref_al, p_ref_al_range)
        #     check_value_exists(ct_pt_pt_al)
        #     check_value_in_range(ct_pt_pt_al, ct_pt_pt_al_range)
        #     check_value_exists(p_ref_al)
        #     check_value_in_range(p_ref_al, p_ref_al_range)
        #     # Argument 3 and 5 are the same from the documentation??

        @given(
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_Fdelta(self, p, lon, lat):
            f_delta_al = Fdelta(p, lon, lat)

            check_value_exists(f_delta_al)
            check_value_in_range(f_delta_al, f_delta_al_range)

        @given(
            in_situ_temperature_st,
        )
        def test_Hill_ratio_at_SP2(self, t):
            hill_ratio_al = Hill_ratio_at_SP2(t)

            check_value_exists(hill_ratio_al)
            check_value_in_range(hill_ratio_al, hill_ratio_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_O2sol(self, SA, CT, p, lon, lat):
            o2_solar_al = O2sol(SA, CT, p, lon, lat)

            check_value_exists(o2_solar_al)
            check_value_in_range(o2_solar_al, o2_solar_al_range)

        @given(
            practical_salinity_st,
            potential_temperature_st,
        )
        def test_O2sol_SP_pt(self, SP, pt):
            o2_solar_al = O2sol_SP_pt(SP, pt)

            check_value_exists(o2_solar_al)
            check_value_in_range(o2_solar_al, o2_solar_al_range)

        @given(
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_SAAR(self, p, lon, lat):
            saar_al = SAAR(p, lon, lat)

            check_value_exists(saar_al)
            check_value_in_range(saar_al, saar_al_range)

        @given(
            practical_salinity_st,
            longitude_st,
            latitude_st,
        )
        def test_SA_from_SP_Baltic(self, SP, lon, lat):
            sa_baltic_al = SA_from_SP_Baltic(SP, lon, lat)

            check_value_exists(sa_baltic_al)
            check_value_in_range(sa_baltic_al, sa_baltic_al_range)

        @given(
            absolute_salinity_st,
            longitude_st,
            latitude_st,
        )
        def test_SP_from_SA_Baltic(self, SA, lon, lat):
            sp_baltic_al = SP_from_SA_Baltic(SA, lon, lat)

            check_value_exists(sp_baltic_al)
            check_value_in_range(sp_baltic_al, sp_baltic_al_range)

        # @given(
        #     Rt_st,
        #     in_situ_temperature_st,
        # )
        # def test_SP_salinometer(self, Rt, t):
        #     practical_salinity_al = SP_salinometer(Rt, t)
        #
        #     check_value_exists(practical_salinity_al)
        #     check_value_in_range(practical_salinity_al, practical_salinity_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_alpha_wrt_t_exact(self, SA, t, p):
            alpha_wrt_t_exact_al = alpha_wrt_t_exact(SA, t, p)

            check_value_exists(alpha_wrt_t_exact_al)
            check_value_in_range(alpha_wrt_t_exact_al, alpha_wrt_t_exact_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_beta_const_t_exact(self, SA, t, p):
            beta_const_t_exact_al = beta_const_t_exact(SA, t, p)

            check_value_exists(beta_const_t_exact_al)
            check_value_in_range(beta_const_t_exact_al, beta_const_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_cabbeling(self, SA, CT, p):
            cabbeling_al = cabbeling(SA, CT, p)

            check_value_exists(cabbeling_al)
            check_value_in_range(cabbeling_al, cabbeling_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_chem_potential_water_t_exact(self, SA, t, p):
            chem_potential_water_t_exact_al = chem_potential_water_t_exact(SA, t, p)

            check_value_exists(chem_potential_water_t_exact_al)
            check_value_in_range(
                chem_potential_water_t_exact_al, chem_potential_water_t_exact_al_range
            )

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_cp_t_exact(self, SA, t, p):
            cp_t_exact_al = cp_t_exact(SA, t, p)

            check_value_exists(cp_t_exact_al)
            check_value_in_range(cp_t_exact_al, cp_t_exact_al_range)

        @given(
            sea_pressure_st,
            longitude_st,
            latitude_st,
        )
        def test_deltaSA_atlas(self, p, lon, lat):
            delta_sa_atlas_al = deltaSA_atlas(p, lon, lat)

            check_value_exists(delta_sa_atlas_al)
            check_value_in_range(delta_sa_atlas_al, delta_sa_atlas_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_dilution_coefficient_t_exact(self, SA, t, p):
            dilution_coefficient_t_exact_al = dilution_coefficient_t_exact(SA, t, p)

            check_value_exists(dilution_coefficient_t_exact_al)
            check_value_in_range(
                dilution_coefficient_t_exact_al, dilution_coefficient_t_exact_al_range
            )

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_dynamic_enthalpy(self, SA, CT, p):
            dynamic_enthalpy_al = dynamic_enthalpy(SA, CT, p)

            check_value_exists(dynamic_enthalpy_al)
            check_value_in_range(dynamic_enthalpy_al, dynamic_enthalpy_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_CT_exact(self, SA, CT, p):
            enthalpy_ct_exact_al = enthalpy_CT_exact(SA, CT, p)

            check_value_exists(enthalpy_ct_exact_al)
            check_value_in_range(enthalpy_ct_exact_al, enthalpy_ct_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_first_derivatives(self, SA, CT, p):
            h_sa_al, h_ct_al = enthalpy_first_derivatives(SA, CT, p)

            check_value_exists(h_sa_al)
            check_value_in_range(h_sa_al, h_sa_al_range)
            check_value_exists(h_ct_al)
            check_value_in_range(h_ct_al, h_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_first_derivatives_CT_exact(self, SA, CT, p):
            h_sa_al, h_ct_al = enthalpy_first_derivatives_CT_exact(SA, CT, p)

            check_value_exists(h_sa_al)
            check_value_in_range(h_sa_al, h_sa_al_range)
            check_value_exists(h_ct_al)
            check_value_in_range(h_ct_al, h_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_second_derivatives(self, SA, CT, p):
            h_sa_sa_al, h_sa_ct_al, h_ct_ct_al = enthalpy_second_derivatives(SA, CT, p)

            check_value_exists(h_sa_sa_al)
            check_value_in_range(h_sa_sa_al, h_sa_sa_al_range)
            check_value_exists(h_sa_ct_al)
            check_value_in_range(h_sa_ct_al, h_sa_ct_al_range)
            check_value_exists(h_ct_ct_al)
            check_value_in_range(h_ct_ct_al, h_ct_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_second_derivatives_CT_exact(self, SA, CT, p):
            h_sa_sa_al, h_sa_ct_al, h_ct_ct_al = enthalpy_second_derivatives_CT_exact(
                SA, CT, p
            )

            check_value_exists(h_sa_sa_al)
            check_value_in_range(h_sa_sa_al, h_sa_sa_al_range)
            check_value_exists(h_sa_ct_al)
            check_value_in_range(h_sa_ct_al, h_sa_ct_al_range)
            check_value_exists(h_ct_ct_al)
            check_value_in_range(h_ct_ct_al, h_ct_ct_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_enthalpy_t_exact(self, SA, t, p):
            enthalpy_t_exact_al = enthalpy_t_exact(SA, t, p)

            check_value_exists(enthalpy_t_exact_al)
            check_value_in_range(enthalpy_t_exact_al, enthalpy_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_entropy_first_derivatives(self, SA, CT):
            eta_sa_al, eta_ct_al = entropy_first_derivatives(SA, CT)

            check_value_exists(eta_sa_al)
            check_value_in_range(eta_sa_al, eta_sa_al_range)
            check_value_exists(eta_ct_al)
            check_value_in_range(eta_ct_al, eta_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_entropy_from_CT(self, SA, CT):
            specific_entropy_al = entropy_from_CT(SA, CT)

            check_value_exists(specific_entropy_al)
            check_value_in_range(specific_entropy_al, specific_entropy_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_entropy_second_derivatives(self, SA, CT):
            eta_sa_sa_al, eta_sa_ct_al, eta_ct_ct_al = entropy_second_derivatives(
                SA, CT
            )

            check_value_exists(eta_sa_sa_al)
            check_value_in_range(eta_sa_sa_al, eta_sa_sa_al_range)
            check_value_exists(eta_sa_ct_al)
            check_value_in_range(eta_sa_ct_al, eta_sa_ct_al_range)
            check_value_exists(eta_ct_ct_al)
            check_value_in_range(eta_ct_ct_al, eta_ct_ct_al_range)

        # @given(
        #     bulk_absolute_salinity_st,
        #     bulk_enthalpy_st,
        #     sea_pressure_st,
        # )
        # def test_frazil_properties(self, SA_bulk, h_bulk, p):
        #     absolute_salinity_final_al, conservative_temperature_final_al, w_ih_final_al = frazil_properties(SA_bulk, h_bulk, p)
        #
        #     check_value_exists(absolute_salinity_final_al)
        #     check_value_in_range(absolute_salinity_final_al, absolute_salinity_final_al_range)
        #     check_value_exists(conservative_temperature_final_al)
        #     check_value_in_range(conservative_temperature_final_al, conservative_temperature_final_al_range)
        #     check_value_exists(w_ih_final_al)
        #     check_value_in_range(w_ih_final_al, w_ih_final_al_range)

        # @given(
        #     bulk_absolute_salinity_st,
        #     bulk_potential_enthalpy_st,
        #     sea_pressure_st,
        # )
        # def test_frazil_properties_potential(self, SA_bulk, h_pot_bulk, p):
        #     absolute_salinity_final_al, conservative_temperature_final_al, w_ih_final_al = frazil_properties_potential(SA_bulk, h_pot_bulk, p)
        #
        #     check_value_exists(absolute_salinity_final_al)
        #     check_value_in_range(absolute_salinity_final_al, absolute_salinity_final_al_range)
        #     check_value_exists(conservative_temperature_final_al)
        #     check_value_in_range(conservative_temperature_final_al, conservative_temperature_final_al_range)
        #     check_value_exists(w_ih_final_al)
        #     check_value_in_range(w_ih_final_al, w_ih_final_al_range)

        # @given(
        #     bulk_absolute_salinity_st,
        #     bulk_potential_enthalpy_st,
        #     sea_pressure_st,
        # )
        # def test_frazil_properties_potential_poly(self, SA_bulk, h_pot_bulk, p):
        #     absolute_salinity_final_al, conservative_temperature_final_al, w_ih_final_al = frazil_properties_potential_poly(SA_bulk, h_pot_bulk, p)
        #
        #     check_value_exists(absolute_salinity_final_al)
        #     check_value_in_range(absolute_salinity_final_al, absolute_salinity_final_al_range)
        #     check_value_exists(conservative_temperature_final_al)
        #     check_value_in_range(conservative_temperature_final_al, conservative_temperature_final_al_range)
        #     check_value_exists(w_ih_final_al)
        #     check_value_in_range(w_ih_final_al, w_ih_final_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            mass_fraction_of_ice_st,
        )
        def test_frazil_ratios_adiabatic(self, SA, p, w_Ih):
            (
                dsa_dct_frazil_al,
                dsa_dp_frazil_al,
                dct_dp_frazil_al,
            ) = frazil_ratios_adiabatic(SA, p, w_Ih)

            check_value_exists(dsa_dct_frazil_al)
            check_value_in_range(dsa_dct_frazil_al, dsa_dct_frazil_al_range)
            check_value_exists(dsa_dp_frazil_al)
            check_value_in_range(dsa_dp_frazil_al, dsa_dp_frazil_al_range)
            check_value_exists(dct_dp_frazil_al)
            check_value_in_range(dct_dp_frazil_al, dct_dp_frazil_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            mass_fraction_of_ice_st,
        )
        def test_frazil_ratios_adiabatic_poly(self, SA, p, w_Ih):
            (
                dsa_dct_frazil_al,
                dsa_dp_frazil_al,
                dct_dp_frazil_al,
            ) = frazil_ratios_adiabatic_poly(SA, p, w_Ih)

            check_value_exists(dsa_dct_frazil_al)
            check_value_in_range(dsa_dct_frazil_al, dsa_dct_frazil_al_range)
            check_value_exists(dsa_dp_frazil_al)
            check_value_in_range(dsa_dp_frazil_al, dsa_dp_frazil_al_range)
            check_value_exists(dct_dp_frazil_al)
            check_value_in_range(dct_dp_frazil_al, dct_dp_frazil_al_range)

        @given(
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_gibbs_ice_part_t(self, t, p):
            gibbs_ice_part_t_al = gibbs_ice_part_t(t, p)

            check_value_exists(gibbs_ice_part_t_al)
            check_value_in_range(gibbs_ice_part_t_al, gibbs_ice_part_t_al_range)

        @given(
            potential_temperature_with_reference_to_zero_pressure_st,
        )
        def test_gibbs_ice_pt0(self, pt0):
            gibbs_ice_part_pt0_al = gibbs_ice_pt0(pt0)

            check_value_exists(gibbs_ice_part_pt0_al)
            check_value_in_range(gibbs_ice_part_pt0_al, gibbs_ice_part_pt0_al_range)

        @given(
            potential_temperature_with_reference_to_zero_pressure_st,
        )
        def test_gibbs_ice_pt0_pt0(self, pt0):
            gibbs_ice_pt0_pt0_al = gibbs_ice_pt0_pt0(pt0)

            check_value_exists(gibbs_ice_pt0_pt0_al)
            check_value_in_range(gibbs_ice_pt0_pt0_al, gibbs_ice_pt0_pt0_al_range)

        @given(
            latitude_st,
            sea_pressure_st,
        )
        def test_grav(self, lat, p):
            grav_al = grav(lat, p)

            check_value_exists(grav_al)
            check_value_in_range(grav_al, grav_al_range)

        # @given(
        #     shape_st,
        #     axis_st,
        #     order_st,
        # )
        # def test_indexer(self, shape, axis, order):
        #     result = indexer(shape, axis, order)
        #
        #     check_value_exists(result)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_kappa_t_exact(self, SA, t, p):
            kappa_t_exact_al = kappa_t_exact(SA, t, p)

            check_value_exists(kappa_t_exact_al)
            check_value_in_range(kappa_t_exact_al, kappa_t_exact_al_range)

        @given(
            absolute_salinity_st,
            sea_pressure_st,
        )
        def test_latentheat_melting(self, SA, p):
            latentheat_melting_al = latentheat_melting(SA, p)

            check_value_exists(latentheat_melting_al)
            check_value_in_range(latentheat_melting_al, latentheat_melting_al_range)

        # @given(
        #     f_st,
        # )
        # def test_match_args_return(self, f):
        #     result = match_args_return(f)
        #
        #     check_value_exists(result)

        # @given(
        #     x_st,
        #     y_st,
        #     xi_st,
        #     axis_st,
        # )
        # def test_pchip_interp(self, x, y, xi, axis):
        #     yi_array = pchip_interp(x, y, xi, axis)
        #
        #     check_value_exists(yi_array)
        #     check_value_in_range(yi_array, yi_array_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
            reference_pressure_st,
        )
        def test_pot_rho_t_exact(self, SA, t, p, p_ref):
            pot_rho_t_exact_al = pot_rho_t_exact(SA, t, p, p_ref)

            check_value_exists(pot_rho_t_exact_al)
            check_value_in_range(pot_rho_t_exact_al, pot_rho_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_pt_first_derivatives(self, SA, CT):
            pt_sa_al, pt_ct_al = pt_first_derivatives(SA, CT)

            check_value_exists(pt_sa_al)
            check_value_in_range(pt_sa_al, pt_sa_al_range)
            check_value_exists(pt_ct_al)
            check_value_in_range(pt_ct_al, pt_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_pt_second_derivatives(self, SA, CT):
            pt_sa_sa_al, pt_sa_ct_al, pt_ct_ct_al = pt_second_derivatives(SA, CT)

            check_value_exists(pt_sa_sa_al)
            check_value_in_range(pt_sa_sa_al, pt_sa_sa_al_range)
            check_value_exists(pt_sa_ct_al)
            check_value_in_range(pt_sa_ct_al, pt_sa_ct_al_range)
            check_value_exists(pt_ct_ct_al)
            check_value_in_range(pt_ct_ct_al, pt_ct_ct_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho_first_derivatives(self, SA, CT, p):
            rho_sa_al, rho_ct_al, rho_p_al = rho_first_derivatives(SA, CT, p)

            check_value_exists(rho_sa_al)
            check_value_in_range(rho_sa_al, rho_sa_al_range)
            check_value_exists(rho_ct_al)
            check_value_in_range(rho_ct_al, rho_ct_al_range)
            check_value_exists(rho_p_al)
            check_value_in_range(rho_p_al, rho_p_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho_first_derivatives_wrt_enthalpy(self, SA, CT, p):
            rho_sa_al, rho_h_al = rho_first_derivatives_wrt_enthalpy(SA, CT, p)

            check_value_exists(rho_sa_al)
            check_value_in_range(rho_sa_al, rho_sa_al_range)
            check_value_exists(rho_h_al)
            check_value_in_range(rho_h_al, rho_h_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho_second_derivatives(self, SA, CT, p):
            (
                rho_sa_sa_al,
                rho_sa_ct_al,
                rho_ct_ct_al,
                rho_sa_p_al,
                rho_ct_p_al,
            ) = rho_second_derivatives(SA, CT, p)

            check_value_exists(rho_sa_sa_al)
            check_value_in_range(rho_sa_sa_al, rho_sa_sa_al_range)
            check_value_exists(rho_sa_ct_al)
            check_value_in_range(rho_sa_ct_al, rho_sa_ct_al_range)
            check_value_exists(rho_ct_ct_al)
            check_value_in_range(rho_ct_ct_al, rho_ct_ct_al_range)
            check_value_exists(rho_sa_p_al)
            check_value_in_range(rho_sa_p_al, rho_sa_p_al_range)
            check_value_exists(rho_ct_p_al)
            check_value_in_range(rho_ct_p_al, rho_ct_p_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_rho_second_derivatives_wrt_enthalpy(self, SA, CT, p):
            rho_sa_sa_al, rho_sa_h_al, rho_h_h_al = rho_second_derivatives_wrt_enthalpy(
                SA, CT, p
            )

            check_value_exists(rho_sa_sa_al)
            check_value_in_range(rho_sa_sa_al, rho_sa_sa_al_range)
            check_value_exists(rho_sa_h_al)
            check_value_in_range(rho_sa_h_al, rho_sa_h_al_range)
            check_value_exists(rho_h_h_al)
            check_value_in_range(rho_h_h_al, rho_h_h_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_sound_speed_t_exact(self, SA, t, p):
            sound_speed_t_exact_al = sound_speed_t_exact(SA, t, p)

            check_value_exists(sound_speed_t_exact_al)
            check_value_in_range(sound_speed_t_exact_al, sound_speed_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_first_derivatives(self, SA, CT, p):
            v_sa_al, v_ct_al, c_p_al = specvol_first_derivatives(SA, CT, p)

            check_value_exists(v_sa_al)
            check_value_in_range(v_sa_al, v_sa_al_range)
            check_value_exists(v_ct_al)
            check_value_in_range(v_ct_al, v_ct_al_range)
            check_value_exists(c_p_al)
            check_value_in_range(c_p_al, c_p_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_first_derivatives_wrt_enthalpy(self, SA, CT, p):
            v_sa_wrt_h_al, v_h_al = specvol_first_derivatives_wrt_enthalpy(SA, CT, p)

            check_value_exists(v_sa_wrt_h_al)
            check_value_in_range(v_sa_wrt_h_al, v_sa_wrt_h_al_range)
            check_value_exists(v_h_al)
            check_value_in_range(v_h_al, v_h_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_second_derivatives(self, SA, CT, p):
            (
                v_sa_sa_al,
                v_sa_ct_al,
                v_ct_ct_al,
                v_sa_p_al,
                v_ct_p_al,
            ) = specvol_second_derivatives(SA, CT, p)

            check_value_exists(v_sa_sa_al)
            check_value_in_range(v_sa_sa_al, v_sa_sa_al_range)
            check_value_exists(v_sa_ct_al)
            check_value_in_range(v_sa_ct_al, v_sa_ct_al_range)
            check_value_exists(v_ct_ct_al)
            check_value_in_range(v_ct_ct_al, v_ct_ct_al_range)
            check_value_exists(v_sa_p_al)
            check_value_in_range(v_sa_p_al, v_sa_p_al_range)
            check_value_exists(v_ct_p_al)
            check_value_in_range(v_ct_p_al, v_ct_p_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_second_derivatives_wrt_enthalpy(self, SA, CT, p):
            (
                v_sa_sa_wrt_h_al,
                v_sa_h_al,
                v_h_h_al,
            ) = specvol_second_derivatives_wrt_enthalpy(SA, CT, p)

            check_value_exists(v_sa_sa_wrt_h_al)
            check_value_in_range(v_sa_sa_wrt_h_al, v_sa_sa_wrt_h_al_range)
            check_value_exists(v_sa_h_al)
            check_value_in_range(v_sa_h_al, v_sa_h_al_range)
            check_value_exists(v_h_h_al)
            check_value_in_range(v_h_h_al, v_h_h_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_specvol_t_exact(self, SA, t, p):
            spec_vol_t_exact_al = specvol_t_exact(SA, t, p)

            check_value_exists(spec_vol_t_exact_al)
            check_value_in_range(spec_vol_t_exact_al, spec_vol_t_exact_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_spiciness0(self, SA, CT):
            spiciness0_al = spiciness0(SA, CT)

            check_value_exists(spiciness0_al)
            check_value_in_range(spiciness0_al, spiciness0_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_spiciness1(self, SA, CT):
            spiciness1_al = spiciness1(SA, CT)

            check_value_exists(spiciness1_al)
            check_value_in_range(spiciness1_al, spiciness1_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
        )
        def test_spiciness2(self, SA, CT):
            spiciness2_al = spiciness2(SA, CT)

            check_value_exists(spiciness2_al)
            check_value_in_range(spiciness2_al, spiciness2_al_range)

        @given(
            absolute_salinity_st,
            in_situ_temperature_st,
            sea_pressure_st,
        )
        def test_t_deriv_chem_potential_water_t_exact(self, SA, t, p):
            chem_potential_water_d_al = t_deriv_chem_potential_water_t_exact(SA, t, p)

            check_value_exists(chem_potential_water_d_al)
            check_value_in_range(
                chem_potential_water_d_al, chem_potential_water_d_al_range
            )

        @given(
            absolute_salinity_st,
            sea_pressure_st,
            saturation_fraction_st,
        )
        def test_t_freezing_poly(self, SA, p, saturation_fraction):
            t_freezing_al = t_freezing_poly(SA, p, saturation_fraction)

            check_value_exists(t_freezing_al)
            check_value_in_range(t_freezing_al, t_freezing_al_range)

        @given(
            absolute_salinity_st,
            conservative_temperature_st,
            sea_pressure_st,
        )
        def test_thermobaric(self, SA, CT, p):
            thermobaric_al = thermobaric(SA, CT, p)

            check_value_exists(thermobaric_al)
            check_value_in_range(thermobaric_al, thermobaric_al_range)
