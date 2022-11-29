from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class HACCFieldInfo(FieldInfoContainer):
    known_particle_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ('id', ('', ['particle_identifier'], None)),
        ('mass', ('code_mass', ['mass','to','particle_mass'], None)),
        ('x', ('code_length', ['particle_position_x'], None)),
        ('y', ('code_length', ['particle_position_y'], None)),
        ('z', ('code_length', ['particle_position_z'], None)),
        ('vx', ('code_velocity', ['particle_velocity_x'], None)),
        ('vy', ('code_velocity', ['particle_velocity_y'], None)),
        ('vz', ('code_velocity', ['particle_velocity_z'], None)),
        ('uu', ('code_velocity**2', ['specific_thermal_energy'], None)),
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


    # def setup_particle_fields(self, ptype, *args, **kwargs):

    #     #convert hh to actual smoothing length
    #     def _smoothing_length(data, field):
    #         return data[(ptype, 'hh')] * 4.0 # hh = SL / nperSL, with n = 4 as compile-time parameter in HACC.  Update as required but cant read it in
    #     self.add_field((ptype, 'smoothing_length'), 
    #                     sampling_type='particle',
    #                     function=partial(_smoothing_length, i=i),
    #                     units='code_length*h')
    #     super().setup_particle_fields(ptype, *args, **kwargs)
    #     if ptype in ("PartType0", "Gas"):
    #         self.setup_gas_particle_fields(ptype)


    def setup_gas_particle_fields(self, ptype):
        def _temperature(field, data):
            gamma = 5.0/3.0
            ret = data[ptype, 'uu'] * (gamma - 1.0) * data.ds.mu * mp / kb
            return ret.in_units(self.ds.unit_system['temperature'])
        self.add_field((ptype, 'Temperature'),
            sampling_type='particle',
            function=_temperature,
            units = self.ds.unit_system['temperature'])

        self.alias((ptype, 'temperature'), (ptype, 'Temperature'))
        self.alias(('gas','temperature'), (ptype, 'Temperature'))


    # def setup_fluid_fields(self):
    #     # Here we do anything that might need info about the dataset.
    #     # You can use self.alias, self.add_output_field (for on-disk fields)
    #     # and self.add_field (for derived fields).
    #     pass

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
