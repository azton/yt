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
        ('mass', ('code_mass', ['particle_mass'], None)),
        ('x', ('code_length', ['particle_position_x'], None)),
        ('y', ('code_length', ['particle_position_y'], None)),
        ('z', ('code_length', ['particle_position_z'], None)),
        ('vx', ('code_velocity', ['particle_velocity_x'], '$\\nu_x$')),
        ('vy', ('code_velocity', ['particle_velocity_y'], '$\\nu_y$')),
        ('vz', ('code_velocity', ['particle_velocity_z'], '$\\nu_z$')),
        ('uu', ('erg/g*h', ['thermal_energy'], None)),
        ('rho', ('code_mass/code_length**3', ['density'], 'Density')),
        ('mu', ('', ['mean_molecular_weight'], None)),
        ('zmet', ('code_metallicity', ['metallicity'], "$Z_{ABS}$")),
        ('yhe', ('', ['helium_fraction'], None)),
        ('phi', ('code_velocity**2', ['gravitational_potential'], None)),
        ('hh', ('Mpc/h',['smoothing_length'], None))
        

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
        if (ptype, "Temperature") not in self.ds.field_list:

            def _temperature(field, data):
                gamma = 5.0/3.0
                mp = data.ds.quan(1.6726e-27, 'kg')
                kb = data.ds.quan(1.3806e-16, 'erg/K')
                ret = data[ptype, 'uu'].to('erg/g') \
                        * (gamma - 1.0) * data[ptype,'mu'] \
                        * mp / kb * 1e10
                return ret.in_units(self.ds.unit_system['temperature'])
            self.add_field((ptype, 'Temperature'),
                sampling_type='particle',
                function=_temperature,
                units = self.ds.unit_system['temperature'])

            self.alias((ptype, 'temperature'), (ptype, 'Temperature'))
            self.alias(("gas", "temperature"), (ptype, "Temperature"))
        if ptype=="Gas":
            self.alias(('gas', "density"), (ptype, 'density'))
            new_ptype='gas'
            def _number_density(field, data):
                mp = data.ds.quan(1.6726e-27, 'kg')
                hfrac = 1-data[new_ptype, 'yhe']
                hefrac = data[new_ptype, 'yhe']
                return data[new_ptype, 'density']*(hfrac/mp + hefrac/(4*mp)) 
            def _h_fraction(field, data):
                return 1-data[new_ptype, 'yhe']
            def _he_fraction(field, data):
                return data[new_ptype, 'yhe']
            def _h_density(field, data):
                mp = data.ds.quan(1.6726e-27, 'kg')
                hfrac = 1-data[new_ptype, 'yhe']
                return data[new_ptype, 'density'] * hfrac
            def _he_density(field, data):
                mp = data.ds.quan(1.6726e-27, 'kg')
                hefrac = data[new_ptype, 'yhe']
                return data[new_ptype, 'density'] * hefrac
            self.add_field((new_ptype, 'number_density'),
                sampling_type='particle',
                function=_number_density,
                units = 'cm**-3')
            self.add_field((new_ptype, 'H_fraction'),
                sampling_type='particle',
                function=_h_fraction,
                units = '')
            self.add_field((new_ptype, 'He_fraction'),
                sampling_type='particle',
                function=_he_fraction,
                units = '')
            self.add_field((new_ptype, 'H_density'),
                sampling_type='particle',
                function=_h_density,
                units='g/cm**3')
            self.add_field((new_ptype, 'He_density'),
                sampling_type='particle',
                function=_he_density,
                units='g/cm**3')
            
            # self.alias(('gas', "number_density"), (ptype, 'number_density'))
            # self.alias(('gas', "H_fraction"), (ptype, 'H_fraction'))
            # self.alias(('gas', "He_fraction"), (ptype, 'He_fraction'))
            # self.alias(('gas', "H_density"), (ptype, 'H_density'))
            # self.alias(('gas', "He_density"), (ptype, 'He_density'))