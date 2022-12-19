from yt.fields.field_info_container import FieldInfoContainer
from yt.frontends.sph.fields import SPHFieldInfo

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class HACCCatalogFieldInfo(FieldInfoContainer):
    known_particle_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    ('gal_idx', ('', ['pariticle_index'], None)),
    ('gal_tag', ('', [], None)),
    ('fof_halo_tag',('', ['host_halo_index'], None)),
    ('gal_radius', ('Mpccm/h', ['aperture_radius'], None)),
    ('gal_dbscan_radius',("Mpccm/h", ['virial_radius'], None)),
    ('gal_2Rhalf_stellar_mass',("Msun/h", [], None)),
    ('gal_half_stellar_rad', ('Mpc/h', ['half_stellar_radius'], None)),
    ('gal_mass', ('Msun/h', [], None)),
    ('gal_dbscan_mstar',('Msun/h', ['particle_mass'], None)),
    ('gal_center_x', ('Mpccm/h', ['particle_position_x'], None)),
    ('gal_center_y',('Mpccm/h', ['particle_position_y'], None)),
    ('gal_center_z', ('Mpccm/h', ['particle_position_z'], None)),
    ('fof_halo_center_x', ('Mpccm/h', ['particle_position_x'], None)),
    ('fof_halo_center_y',('Mpccm/h', ['particle_position_y'], None)),
    ('fof_halo_center_z',('Mpccm/h', ['particle_position_z'], None)),
    ('fof_halo_mass', ('Msun/h', ['particle_mass'], None)),
    ('sod_halo_radius', ('Mpccm/h', ['virial_radius'], None)),
    ('sod_halo_mass', ('Msun/h', ['particle_sod_mass'], None)),
    ('sod_halo_center_x', ('Mpccm/h', ['particle_sod_position_x'], None)),
    ('sod_halo_center_y',('Mpccm/h', ['particle_sod_position_y'], None)),
    ('sod_halo_center_z',('Mpccm/h', ['particle_sod_position_z'], None)),
    ('radius', ('Mpccm/h', [], None)),
    ('central', ('', [], None)),
    ('vel_disp', ('kmcm/h/s', [], None)),
    ('x', ('Mpccm/h', ['particle_position_x'], None)),
    ('y', ('Mpccm/h', ['particle_position_y'], None)),
    ('z', ('Mpccm/h', ['particle_position_z'], None)),
    ('vx', ('kmcm/s/h', ['particle_velocity_x'], None)),
    ('vy', ('kmcm/s/h', ['particle_velocity_y'], None)),
    ('vz', ('kmcm/s/h', ['particle_velocity_z'], None)),
    ('infall_step', ('',[], None)),
    ('infall_fof_halo_tag', ('', [], None)),
    ('infall_fof_halo_mass', ('Msun/h', [], None)),
    ('merged', ('', [], None))
    )

    known_other_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def __init__(self, ds, field_list, slice_info=None):
        # If you want, you can check self.field_list
        super().__init__(ds, field_list, slice_info=slice_info)


    def setup_particle_fields(self, ptype, *args, **kwargs):
        super().setup_particle_fields(ptype)

           


    def setup_gas_particle_fields(self, ptype):
        pass