import os
os.environ['GENERICIO_NO_MPI'] = 'True'
from importlib.util import find_spec as imp
# if imp('pygio'):
import pygio
legacy = False
# else:
#     import sys
#     sys.path.append('/home/azton/genericio/legacy_python')
#     import genericio as pygio
#     legacy = True
import glob
import numpy as np

import yt.frontends.hacc_catalog.definitions
from yt.utilities.io_handler import BaseIOHandler
from yt.frontends.sph.io import IOHandlerSPH


class HACCCatalogIOHandler(BaseIOHandler):
    # _particle_reader = True
    _dataset_type = "hacc_catalog"
    relevant_core_fnames= [
        'central',
        'vel_disp',
        'radius',
        *['v%s'%s for s in 'xyz'],
        'fof_halo_tag',
        'infall_step',
        'infall_fof_halo_tag',
        'infall_fof_halo_mass',
        'gal_tag',
        'merged',
        *'xyz',

    ]
    relevant_galaxy_fnames = [
    'gal_idx',
    'gal_tag',
    'fof_halo_tag',
    'gal_radius',
    'gal_dbscan_radius',
    'gal_half_stellar_rad',
    'gal_2Rhalf_stellar_mass',
    'gal_mass',
    'gal_center_x',
    'gal_center_y',
    'gal_center_z'
    ]
    relevant_halo_fnames = [
        'fof_halo_tag',
        'fof_halo_center_x',
        'fof_halo_center_y',
        'fof_halo_center_z',
        'fof_halo_mass',
        'sod_halo_radius',
        'sod_halo_mass',
        'sod_halo_center_x',
        'sod_halo_center_y',
        'sod_halo_center_z'
    ]
    # TODO CORE FIELD NAMES
    galaxy_pos = relevant_galaxy_fnames[-3:]
    halo_pos = relevant_halo_fnames[1:4]
    core_pos = ['x','y','z']
    def _count_particles(self, data_file):
        # print('all_io_datafiles',self.ds.data_files)
        # pass
        si, ei = data_file.start, data_file.end
        print(f"catalog _count_particles: {data_file.filename} {si} {ei}, {ei-si} total particles")
        if not legacy:
            f = pygio.PyGenericIO('%s'%data_file.filename )
            nparts = f.read_num_elems()
        else:
            nparts = pygio.get_num_elements('%s'%data_file.filename)
        # mask = f.read(['mask'], print_stats=False)['mask']
        # this effort tries to get the dark matter into its own field
        # nparts = self.dark_matter(mask)
        if None not in (si, ei):
            nparts = np.clip(nparts-si, 0, ei-si)
        print(f"particle count: {nparts}")
        # print('%s read %d particles'%(data_file.filename, nparts))
        pcnt = {}
        for ptype in data_file.ds._particle_types_raw:
            pcnt[ptype] = nparts
        return pcnt

    def _identify_fields(self, data_file):
        # pass
        self.galaxies = False
        self.halos = False
        self.cores = False
        if not legacy:
            f = pygio.PyGenericIO(data_file.filename)
            vnames = f.read_variable_names()
        else:
            vnames = pygio.get_scalars(data_file.filename)[1]
        # print('vnames',vnames)
        fields = []
        ptypes = data_file.ds._particle_types_raw
        if 'gal_dbscan_radius' in vnames:
            self.galaxies = True
            field_names = self.relevant_galaxy_fnames
        elif 'sod_halo_mass_dm' in vnames:
            self.halos = True
            field_names = self.relevant_halo_fnames
        elif 'infall_fof_halo_tag' in vnames:
            self.cores = True
            field_names = self.relevant_core_fnames
        if not self.cores and not self.halos and not self.galaxies:
            raise RuntimeError("No valid fields found")
            print("FIELDS: ", vnames)
        for ptype in ptypes:
            for f in field_names:
                fields.append((ptype, f))
            
        return fields, {}

    def _read_particle_coords(self, chunks, ptf):
        pass
        # # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # # chunks is a list of chunks, and ptf is a dict where the keys are
        # # ptypes and the values are lists of fields.
        # data_files = set()
        # for chunk in chunks:
        #     for obj in chunk.objs:
        #         data_files.update(obj.data_files)
        # for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
        #     si, ei = data_file.start, data_file.end
        #     # print('read_particle_coords: fname = %s'%data_file.filename)
        #     with pygio.PyGenericIO(data_file.file_name) as f:
        #         for ptype in sorted(ptf):
                    
        #             coords = f.read(['x', 'y', 'z','mask'], print_stats=False)
        #             mask = FILTER_DEF[ptype](coords['mask'])[si:ei]
        #             data = coords['x'][si, ei], coords['y'][si, ei], coords['z'][si, ei]
        #             data[0] = data[0][mask]
        #             data[1] = data[1][mask]
        #             data[2] = data[2][mask]

        #             yield ptype, data
                
    def _yield_coordinates(self, data_file, needed_ptype=None):
        # pass
        si, ei = data_file.start, data_file.end
        print(f"catalog _yield_coordinates: {data_file.filename} {si} {ei}, {ei-si} total particles")
        if not legacy:
            f = pygio.PyGenericIO(data_file.filename)
        for ptype, cnt in data_file.total_particles.items():
            if cnt == 0: continue
            if needed_ptype is not None and ptype != needed_ptype:
                continue
            if not legacy:
                if self.galaxies:
                    pp = f.read(self.galaxy_pos, print_stats=False)
                elif self.halos:
                    pp = f.read(self.halo_pos, print_stats=False)
                elif self.cores:
                    pp = f.read(self.core_pos, print_stats=False)
                else:
                    raise(NotImplementedError(data_file.filename))
                pp = [ax[si:ei] for ax in pp.values()]
            else:
                npart = ei-si
                if self.galaxies:
                    pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.galaxy_pos]
                    # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.galaxy_pos]
                elif self.halos:
                    pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.halo_pos]
                    # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.halo_pos]
                elif self.cores:
                    pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.core_pos]
                    # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.core_pos]
                else:
                    raise(NotImplementedError(data_file.filename))
            pp = np.array(pp).T
            # print(pp.shape)
            yield ptype, pp


    def _read_particle_fields(self, chunks, ptf, selector=None):
        # pass
        # # This gets called after the arrays have been allocated.  It needs to
        # # yield ((ptype, field), data) where data is the masked results of
        # # reading ptype, field and applying the selector to the data read in.
        # # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # # you need to do your masking here.
        # # print(type(chunks[0].objs[0].data_files[0]))
        # # print(chunks[0].objs[0].data_files)
        # # print(dir(chunks[0]))
        # data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                for data_file in obj.data_files:
                    si, ei = data_file.start, data_file.end
                    print(f"catalog _read_particle_fields: {data_file.filename} {si} {ei}, {ei-si} total particles")
                    # print(type(si))
                    # print(type(ei))
                    data_return = {}
                    if not legacy:
                        f = pygio.PyGenericIO(data_file.filename)
                    for ptype, field_list in sorted(ptf.items()):
                        if data_file.total_particles[ptype] == 0:
                            continue
                        # if selector is None or getattr(selector, 'is_all_data', False):
                        #     mask = slice(None, None, None)
                        #     mask_sum = data_file.total_particles[ptype]
                        # else:
                        if not legacy:
                            if self.galaxies:
                                coords = f.read(self.galaxy_pos, print_stats=False)
                                mask = selector.select_points(
                                    *[coords[pname][si:ei] for pname in self.galaxy_pos], 0.0
                                )
                            elif self.halos:
                                coords = f.read(self.halo_pos, print_stats=False)
                                mask = selector.select_points(
                                    *[coords[pname][si:ei] for pname in self.halo_pos], 0.0
                                )
                            elif self.cores:
                                coords = f.read(self.core_pos, print_stats=False)
                                mask = selector.select_points(
                                    *[coords[pname][si:ei] for pname in self.core_pos], 0.0
                                )
                            else:
                                raise(NotImplementedError(data_file.filename))
                        else:
                            npart = ei-si
                            if self.galaxies:
                                pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.galaxy_pos]
                                # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.galaxy_pos]
                            elif self.halos:
                                pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.halo_pos]
                                # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.halo_pos]
                            elif self.cores:
                                pp = [pygio.read_globalOffset_scalar(data_file.filename, ax, si, npart) for ax in self.core_pos]
                                # pp = [pygio.read_scalar(data_file.filename, ax)[si:ei] for ax in self.core_pos]
                            else:
                                raise(NotImplementedError(data_file.filename))
                            coords = pp
                            mask = selector.select_points(*coords, 0.0)
                        if mask is not None:
                            mask_sum = mask.sum()
                        if mask is None: continue
                        for field in field_list:
                            # print(si, ei, mask)
                            if not legacy:
                                data = f.read([field], print_stats=False)[field][si:ei]
                            elif legacy:
                                npart = ei-si
                                data=pygio.read_globalOffset_scalar(data_file.filename, field, si, npart)
                                # data = pygio.read_scalar(data_file.filename, field)[si:ei]
                            data = np.array(data)[mask]
                            
                        # data_return[(ptype, field)] = data
                            # print('%s: read_fields read %d'%(data_file.filename, data.size))
                            yield ((ptype, field), data)


    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        pass
    
    def _read_particle_data(self, data_file, ptf, selector=None):
        pass
        # open relevant file given by data_file.file_name
        # probably have to open the subgrid file as well...  
