from yt.fields.field_info_container import FieldInfoContainer
from yt.frontends.sph.fields import SPHFieldInfo

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class HACCFieldInfo(SPHFieldInfo):
    known_particle_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ('id', ('', ['particle_identifier'], None)),
        ('mass', ('code_mass', ['mass','to','particle_mass'], None)),
        ('x', ('code_length', ['particle_position_x'], None)),
        ('y', ('code_length', ['particle_position_y'], None)),
        ('z', ('code_length', ['particle_position_z'], None)),
        ('vx', ('code_velocity', ['particle_velocity_x'], '$\\nu_x$')),
        ('vy', ('code_velocity', ['particle_velocity_y'], '$\\nu_y$')),
        ('vz', ('code_velocity', ['particle_velocity_z'], '$\\nu_z$')),
        ('uu', ('J/Msun', ['specific_thermal_energy'], None)),
        ('rho', ('code_mass/code_length**3', ['density'], 'Density')),
        ('mu', ('', ['mean_molecular_weight'], None)),
        ('zmet', ('code_metallicity', ['metallicity'], None)),
        ('yhe', ('', ['helium_fraction'], None)),
        ('phi', ('code_velocity**2', ['gravitational_potential'], None)),
        ('hh', ('code_length',['smoothing_length'], None))
        

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
        if ptype in ("Gas", "Wind", "SF_Gas"):
            self.setup_gas_particle_fields(ptype)
           


    def setup_gas_particle_fields(self, ptype):
        def _temperature(field, data):
            gamma = 5.0/3.0
            ret = data[ptype, 'uu'] * data[ptype, 'mass'] * (gamma - 1.0) * data['mean_molecular_weight'] * mp / kb
            return ret.in_units(self.ds.unit_system['temperature'])
        self.add_field((ptype, 'Temperature'),
            sampling_type='particle',
            function=_temperature,
            units = self.ds.unit_system['temperature'])

        self.alias((ptype, 'temperature'), (ptype, 'Temperature'))
