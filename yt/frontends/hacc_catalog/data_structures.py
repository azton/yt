import os
os.environ['GENERICIO_NO_MPI'] = 'True'
import weakref
import numpy as np
from importlib.util import find_spec as imp
# if imp('pygio'):
import pygio as pg
legacy = False
# else:
#     import sys
#     sys.path.append('/home/azton/genericio/legacy_python')
#     import genericio as pg
#     legacy = True
import glob
from typing import Type
from yt.data_objects.static_output import ParticleFile
from yt.fields.field_info_container import FieldInfoContainer
from yt.geometry.geometry_handler import Index
from yt.data_objects.static_output import Dataset
from yt.data_objects.static_output import (
    ParticleDataset,
    ParticleFile,
    validate_index_order,
)
from yt.utilities.cosmology import Cosmology
from yt.frontends.halo_catalog.data_structures import HaloDataset
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.frontends.hacc_catalog.fields import HACCCatalogFieldInfo
from yt.funcs import setdefaultattr

class HACCCatalogIndex(ParticleIndex):
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
    @property
    def chunksize(self):
        return 256**3
    def _initialize_frontend_specific(self):
        super()._initialize_frontend_specific()
        self.io._float_type = self.ds._header.float_type

    def _setup_filenames(self):
        template = self.dataset.filename_template
        ndoms = self.dataset.file_count
        if ndoms==0:
            print('dataset file_count is zero!!')
            raise(IOError(self.dataset.filename_template))
        cls = self.dataset._file_class
        self.data_files = []
        data_files = glob.glob(f'{self.dataset.filename_template}*')

        # print('catalog:_setup_filenames found ')
        # for f in data_files: print(f)
        fi = 0

        for i in range(int(ndoms)):
            start = 0 
            if self.chunksize > 0:
                end = start + self.chunksize
            else:
                end = None
            while True:
                try:
                    _filename = data_files[i]
                    df = cls(self.dataset, self.io, _filename, fi, (start, end))
                except FileNotFoundError:
                    mylog.warning(
                        "Failed to load '%s' (missing file or directory)", _filename
                    )
                    break
                if max(df.total_particles.values()) == 0:
                    break
                fi += 1
                self.data_files.append(df)
                if end is None:
                    break
                start = end
                end += self.chunksize
class HACCCatalogFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id, range=None):
        super().__init__(ds, io, filename, file_id, range)

class HACCGenericIOHeader():
    def __init__(self):
        self.float_type = 'float32'

class HACCCatalogDataset(ParticleDataset):
    _file_class = HACCCatalogFile
    _index_class = HACCCatalogIndex
    _field_info_class = HACCCatalogFieldInfo
    _particle_types_raw = ('particle',)
    
    def __init__(
        self,
        filename,
        index_order=None,
        dataset_type='hacc_catalog',
        # units_override=None,
        unit_system='cgs',
        default_species_fields=None
    ):
        # self.index_order = validate_index_order(index_order)
        super().__init__(
            filename,
            dataset_type,
            # units_override=units_override,
            unit_system=unit_system,
            # default_species_fields=default_species_fields
        )
        # self.filename = filename
        sn = os.path.split(self.filename)[-1].split('-')[1]
        self.current_timestep = int(sn.split('.')[0])
        try:
            prep = os.path.split(self.filename)[0]
            self.parameter_file = glob.glob(f"{prep}/*{self.current_timestep}*.params")[0]
            print('Read catalog parameter file as %s'%self.parameter_file)
        except:
            print(f'No parameter file found.  Use the parameter file from the original simulation,') 
            print(f'placed in the same directory the input catalog as, e.g., {self.current_timestep}.params')
            print(f"We searched{prep}/*{self.current_timestep}*.params")
            exit(OSError)
        sn = os.path.split(self.filename)[-1].split('-')[1]
        self.current_timestep = int(sn.split('.')[0])
        self.parameters = {}
        self._header = HACCGenericIOHeader()

    def _set_code_unit_attributes(self):

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
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def calculate_current_redshift(self):
        
        stepn = self.current_timestep
        alpha = 1

        z_i = 200
        z_f = 0
        total_steps = 625

        ai = 1.0 / (1.0 + z_i)
        af = 1.0 / (1.0 + z_f)
        steps = np.arange(1, total_steps+1).tolist()
        a = np.linspace(np.power(ai, alpha), np.power(af, alpha), total_steps+1)
        a = np.power(a, 1.0/alpha)

        a = a[1:]
        z = 1.0 / a - 1.0
        redshift = z[stepn]
        self.current_redshift = redshift

    def parse_parameters(self):
        with open(self.parameter_file, 'r') as f:
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
                # print(k)
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
        if 'galaxy' in self.filename or 'halo' in self.filename:
            self.filename_template = self.filename + '#'
        else:
            self.filename_template = self.filename
        # print('filename_template = ', self.filename_template)
        sn = os.path.split(self.filename)[-1].split('-')[1]
        self.current_timestep = int(sn.split('.')[0])
        try:
            prep = os.path.split(self.filename)[0]
            self.parameter_file = glob.glob(f"{prep}/*{self.current_timestep}*.params")[0]
            print('Read catalog parameter file as %s'%self.parameter_file)
        except:
            print(f'No parameter file found.  Use the parameter file from the original simulation,') 
            print(f'placed in the same directory the input catalog as, e.g., {self.current_timestep}.params')
            print(f"We searched{prep}/*{self.current_timestep}*.params")
            exit()
        files = sorted(glob.glob('%s*'%self.filename_template))
        # print('parse parameter file:', files)
        self.parse_parameters()
        if not legacy:
            f = pg.PyGenericIO(files[0])
            self.domain_left_edge = np.array(f.read_phys_origin())
            self.domain_right_edge = np.array(f.read_phys_scale()) 
        else:
            self.domain_left_edge = np.array([0.0, 0.0, 0.0])
            box = self.parameters['rl']
            self.domain_right_edge= np.array([box, box, box])
        self.cosmological_simulation = 1
    
        self.dimensionality = 3
        self.domain_dimensions = np.ones(3, 'int32')
        self._periodicity = [True]*3
        self.file_count = len(glob.glob('%s*'%self.filename_template))
        self.hubble_constant = self.parameters['hubble']
        self.omega_matter = self.parameters['deut'] * self.hubble_constant**2 + self.parameters['omega_cdm']
        self.omega_lambda = 1.0 - self.omega_matter
        co = Cosmology(
            hubble_constant = self.hubble_constant,
            omega_matter = self.omega_matter,
            omega_lambda = self.omega_lambda
        )
        self.current_time = co.t_from_z(self.current_redshift)

        # self.current_time = self.parameters['redshifts'][0]
        self.first_out_file = int(os.path.split(files[0])[-1].split('#')[1]) if len(files) > 1 else None

    @classmethod    
    def _is_valid(cls, filename, *args, **kwargs):
        valid = False
        # print(filename)
        if 'properties' in filename:
            try:
                # print('hacc_catalog checking...')
                files = glob.glob('%s'%filename)
                print(files[0])
                if not legacy:
                    f = pg.PyGenericIO(files[0])
                    a= pg.read_num_elems(files[0])
                else:
                    a = pg.get_num_elements(files[0])
                # print(a)
                return True
            except:
                print("Couldn't verify as hacc catalog")
                valid = False
        return valid