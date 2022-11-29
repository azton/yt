import os
import weakref
import numpy as np
import pygio as pg
import glob
from typing import Type
from yt.data_objects.static_output import ParticleFile
from yt.fields.field_info_container import FieldInfoContainer
from yt.frontends.sph.data_structures import SPHDataset, SPHParticleIndex
from yt.geometry.geometry_handler import Index
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex

from .fields import HACCFieldInfo

class HACCBinaryIndex(SPHParticleIndex):
    def __init__(self, ds, dataset_type):
        super().__init__(ds, dataset_type)
        self._initialize_index()

    def _initialize_index(self):
        # Normally this function is called during field detection. We call it
        # here because we need to know which fields exist on-disk so that we can
        # read in the smoothing lengths for SPH data before we construct the
        # Morton bitmaps.
        self._detect_output_fields()
        super()._initialize_index()

    def _initialize_frontend_specific(self):
        super()._initialize_frontend_specific()
        self.io._float_type = self.ds._header.float_type

class HACCGenericIOFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id, range=None):
        super().__init__(ds, io, filename, file_id, range)

class HACCGenericIOHeader():
    def __init__(self, filename):
        self.filename = filename
        f = pg.PyGenericIO('.'.join(filename.split('.')[:-1])) # filename is the parameter file--need to load actual data
        self.float_type = f.read_variable_dtypes()['x']
        del(f)

class HACCDataset(SPHDataset):
    _index_class: Type[Index] = HACCBinaryIndex
    _file_class: Type[ParticleFile] = HACCGenericIOFile
    _field_info_class: Type[FieldInfoContainer] = HACCFieldInfo
    _particle_mass_name = "mass"
    _sph_ptypes = ('Gas','PartType0')

    # _particle_coordinates_name = 'x'

    def __init__(
        self,
        filename,
        dataset_type="hacc_genericio",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
        index_order=None,
        index_filename=None,
        kdtree_filename=None,
    ):
        self.fluid_types += ("HACC",)
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )
        self.storage_filename = None
        if not os.path.exists(filename):
            print("Parameter file %s not found!"%filename)
            raise RuntimeError
        self._header = HACCGenericIOHeader(filename)
        # self._header.float_type = self.header.float_type
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        self.length_unit = self.quan(1.0, 'Mpccm/h')
        # self.mass_unit = self.quan(1.0, "g")
        self.mass_unit = self.quan(1.0, 'Msun/h')
        # self.time_unit = self.quan(1.0, "s")
        self.time_unit = self.quan(1.0, 's')
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        self.velocity_unit = self.quan(1.0, 'kmcm/s/h')
        # self.magnetic_unit = self.quan(1.0, "gauss")
        #
        # If your frontend uses SI EM units, set magnetic units like this
        # instead:
        # self.magnetic_unit = self.quan(1.0, "T")

        # this minimalistic implementation fills the requirements for
        # this frontend to run, change it to make it run _correctly_ !
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))


    def calculate_current_redshift(self):
        z_i = self.parameters['z_in']
        z_f = 0.0
        nsteps = int(self.parameters['n_steps'])
        alpha = 1

        ai = 1.0 / (1.0 + z_i)
        af = 1.0 / (1.0 + z_f)
        steps = np.arange(1, nsteps+1).tolist()
        a = np.linspace(np.power(ai, alpha), np.power(af, alpha), nsteps+1)
        a = np.power(a, 1.0/alpha)

        a = a[1:]
        z = 1.0 / a - 1.0
        indices = [steps.index(t) for t in self.parameters['full_alive_dump']]
        redshifts = [z[t] for t in indices]
        self.parameters['redshifts'] = redshifts
        

    def parse_parameters(self):
        self.current_timestep = int(self.filename.split('.')[-2]) # expecting <parname>.<step number>.params
        with open(self.filename, 'r') as f:
            parameters = {}
            multipar = ['refresh',
                        'full_alive_dump',
                        'pk_dump',
            ]

            for l in f:
                if l.startswith('#') or l.isspace():
                    continue
                # print('parsing _%s_'%l.strip())
                k = l.split()[0].lower()
                print(k)
                if len(l.split()) == 1: # no value corresponding to key
                    v = None
                elif k not in multipar:
                    v = l.split()[1]
                    try:
                        v = float(v)
                    except:
                        v = str(v)   
                else:
                    v = [int(i) for i in l.split()[1:]]
                parameters[k] = v
            self.parameters = parameters
            self.calculate_current_redshift()


    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        self.parse_parameters()

        self.domain_left_edge = np.zeros(3)
        self.domain_right_edge = np.ones(3) * self.parameters['rl']
        self.dimensionality = 3
        self.domain_dimensions = np.ones(3, 'int32') 
        self._periodicity = [True] * 3
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #  Non cosmological isnt an option here.  Will add when/if we start doing that in HACC
        self.cosmological_simulation = 1
        self.current_redshift = self.parameters['redshifts'][self.parameters['full_alive_dump'].index(self.current_timestep)]
        self.hubble_constant = self.parameters['hubble']
        self.omega_matter = self.parameters['deut'] * self.hubble_constant**2 + self.parameters['omega_cdm']
        self.omega_lambda = 1.0 - self.omega_matter
        self.current_time = self.current_redshift

        self.filename_template = os.path.split('.'.join(self.parameter_filename.split('.')[:-1])+"#")[-1]
        print('template = %s'%self.filename_template)
        self.file_count = len(glob.glob('%s*'%self.filename_template)) # exclude .params file
        self.filename_template = f"{self.filename_template}%(num)s"
        # required. Change this if need be.

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # The only valid files here are genericIO; which are unique to HACC so ...
        valid = False
        with open(filename, 'r') as f:
            for l in f:
                if  'HACC_HEADER_VERSION' in l:
                    valid = True
        return valid

