"""
Enzo frontend tests using moving7



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load
from yt.frontends.enzo.api import EnzoDataset

_fields = ("temperature", "density", "velocity_magnitude",
           "velocity_divergence")

def check_color_conservation(pf):
    species_names = pf.field_info.species_names
    dd = pf.all_data()
    dens_yt = dd["density"].copy()
    # Enumerate our species here
    for s in sorted(species_names):
        if s == "El": continue
        dens_yt -= dd["%s_density" % s]
    dens_yt -= dd["metal_density"]
    delta_yt = np.abs(dens_yt / dd["density"])

    # Now we compare color conservation to Enzo's color conservation
    dd = pf.all_data()
    dens_enzo = dd["Density"].copy()
    for f in sorted(pf.field_list):
        if not f[1].endswith("_Density") or \
               f[1].startswith("Dark_Matter_")  or \
               f[1].startswith("Electron_") or \
               f[1].startswith("SFR_") or \
               f[1].startswith("Forming_Stellar_") or \
               f[1].startswith("Star_Particle_"):
            continue
        dens_enzo -= dd[f]
    delta_enzo = np.abs(dens_enzo / dd["Density"])
    return assert_almost_equal, delta_yt, delta_enzo

m7 = "DD0010/moving7_0010"
@requires_pf(m7)
def test_moving7():
    pf = data_dir_load(m7)
    yield assert_equal, str(pf), "moving7_0010"
    for test in small_patch_amr(m7, _fields):
        test_moving7.__name__ = test.description
        yield test

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
@requires_pf(g30, big_data=True)
def test_galaxy0030():
    pf = data_dir_load(g30)
    yield check_color_conservation(pf)
    yield assert_equal, str(pf), "galaxy0030"
    for test in big_patch_amr(g30, _fields):
        test_galaxy0030.__name__ = test.description
        yield test

hds0 = "rockstar_halos/halos_0.0.bin"
hds1 = "rockstar_halos/halos_0.1.bin"
@requires_pf(hds0)
@requires_pf(hds1)
def test_halo_mass_function():
	hds = data_dir_load(hds0)
	yield assert_equal, str(hds), "halos_0.0.bin"
	for test in hmf_sim_and_analytic(hds0):
		test_halo_mass_function.__name__ = test.description
		yield test

ecp = "enzo_cosmology_plus/DD0046/DD0046"
@requires_pf(ecp, big_data=True)
def test_ecp():
    pf = data_dir_load(ecp)
    # Now we test our species fields
    yield check_color_conservation(pf)
