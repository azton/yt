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
from yt.frontends.hacc.definitions import *
from yt.utilities.io_handler import BaseIOHandler
from yt.frontends.sph.io import IOHandlerSPH


class HACCIOHandler(IOHandlerSPH):
    # _particle_reader = True
    _dataset_type = "hacc_genericio"

    def dark_matter(self, mask):
        nf = fdark_matter(mask).sum()
        return nf



    def _count_particles(self, data_file):
        # print('Counting particles: %s'%data_file.filename)
        # print('si, ei = %d, %d'%(data_file.start, data_file.end))
        si, ei = data_file.start, data_file.end
        print(f"_count_particles {data_file.filename} from {si} to {ei} ({ei-si} particles))")
        if not legacy:
            f = pygio.PyGenericIO('%s'%data_file.filename )
            nparts = f.read_num_elems()
        else:
            nparts = pygio.get_num_elements('%s'%data_file.filename)
        print("nparts = ", nparts)
        # mask = f.read(['mask'], print_stats=False)['mask']
        # this effort tries to get the dark matter into its own field
        # nparts = self.dark_matter(mask)
        if None not in (si, ei):
            nparts = np.clip(nparts-si, 0, ei-si)
        # print('%s read %d particles'%(data_file.filename, nparts))
        pcnt = {}
        for ptype in data_file.ds._sph_ptypes:
            pcnt[ptype] = nparts
        return pcnt

    def _identify_fields(self, data_file):
        if not legacy:
            f = pygio.PyGenericIO(data_file.filename)
            fields = []
            ptypes = data_file.ds._sph_ptypes
            for var in f.get_variables():
                if var.name != 'status':
                    for ptype in ptypes:
                        fields.append((ptype, var.name))
        else:
            fields = []
            ptypes = data_file.ds._sph_ptypes
            for var in pygio.get_scalars(data_file.filename)[1]:
                if var != 'status':
                    for ptype in ptypes:
                        fields.append((ptype, var))
        return fields, {}
    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            print(f"_read_particle_coords {data_file.filename} from {si} to {ei} ({ei-si} particles)")
            # print('read_particle_coords: fname = %s'%data_file.filename)
            for ptype in sorted(ptf):
                    
                if legacy:
                    print('Reading global offset')
                    npart = ei-si
                    # mask = pygio.read_scalar(data_file.filename, 'mask')['mask'][si:ei]
                    mask = pygio.read_globalOffset_scalar(data_file.filename, 'mask', si, npart)
                    mask = FILTER_DEF[ptype](coords['mask'])
                    # data = [pygio.read_scalar(data_file.filename, axis)[axis][si:ei] for axis in 'xyz']
                    data = [pygio.read_globalOffset_scalar(data_file.filename, axis, si, npart) for axis in 'xyz']
                    # data = [ax for ax in pos]
                    data[0] = data[0][mask]
                    data[1] = data[1][mask]
                    data[2] = data[2][mask]
                else:
                    f = pygio.PyGenericIO(data_file.filename)
                    coords = f.read(['x', 'y', 'z','mask'], print_stats=False)
                    mask = FILTER_DEF[ptype](coords['mask'])[si:ei]
                    data = coords['x'][si, ei], coords['y'][si, ei], coords['z'][si, ei]
                    data[0] = data[0][mask]
                    data[1] = data[1][mask]
                    data[2] = data[2][mask]
            
                yield ptype, data
                
    def _yield_coordinates(self, data_file, needed_ptype=None):
        si, ei = data_file.start, data_file.end
        print(f"_yield_coordinates {data_file.filename} from {si} to {ei} ({ei-si} particles)")
        for ptype, cnt in data_file.total_particles.items():
            if cnt == 0: continue
            if needed_ptype is not None and ptype != needed_ptype:
                continue
            if legacy:
                npart = ei-si
                # legacy octree access
                # pos = [pygio.read_scalar(data_file.filename, axis)[axis][si:ei] for axis in 'xyz']
                pos = [pygio.read_globalOffset_scalar(data_file.filename, axis, si, npart) for axis in 'xyz']
                pp = np.array(pos).T
            else:
                f = pygio.PyGenericIO(data_file.filename)
                pp = f.read(['x','y','z'], print_stats=False)
                pp = [ax[si:ei] for ax in pp.values()]
                pp = np.array(pp).T
                # print(pp.shape)
            yield ptype, pp
    def _get_smoothing_length(self, data_file, pdtype, pshape):
        ptype='DarkMatter'
        si, ei = data_file.start, data_file.end 
        print(f"_get_smoothing_length {data_file.filename} from {si} to {ei} ({ei-si} particles)")
        if legacy:
            npart = ei-si
            sl = pygio.read_globalOffset_scalar(data_file.filename, 'hh', si, npart) * 2.0
        else:
            f = pygio.PyGenericIO(data_file.filename)
            sl = f.read(['hh'], print_stats=False)['hh'][si:ei] * 2.0
        return sl

    def _read_particle_fields(self, chunks, ptf, selector=None):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        # print(type(chunks[0].objs[0].data_files[0]))
        # print(chunks[0].objs[0].data_files)
        # print(dir(chunks[0]))
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                for data_file in obj.data_files:
                    si, ei = data_file.start, data_file.end
                    print(f"_read_particle_fields {data_file.filename} from {si} to {ei} ({ei-si} particles)")
                    # print(type(si))
                    # print(type(ei))
                    data_return = {}
                    for ptype, field_list in sorted(ptf.items()):
                        if data_file.total_particles[ptype] == 0:
                            continue
                        # if selector is None or getattr(selector, 'is_all_data', False):
                        #     mask = slice(None, None, None)
                        #     mask_sum = data_file.total_particles[ptype]
                        # else:
                        if legacy:
                            npart = ei-si
                            mask = np.array(pygio.read_globalOffset_scalar(data_file.filename, 'mask', si, npart))
                            # mask = pygio.read_scalar(data_file.filename, 'mask')['mask'][si:ei]
                            pmask = FILTER_DEF[ptype](mask)
                            pos = [np.array(pygio.read_globalOffset_scalar(data_file.filename, axis, si, npart)) for axis in 'xyz']
                            # pos = [pygio.read_scalar(data_file.filename, axis)[axis][si:ei] for axis in 'xyz']
                            hh = pygio.read_globalOffset_scalar(data_file.filename, 'hh', si, npart)
                            # hh = pygio.read_scalar(data_file.filename, 'hh')['hh'][si:ei]
                            sl = hh * 2.0
                            # pos = np.array(pos)
                            print('_read_particle_coords: position shape:', pos[0].shape, mask.shape)
                            data = []
                            data.append(pos[0][mask])
                            data.append(pos[1][mask])
                            data.append(pos[2][mask])
                            data = np.array(data)
                            pmask = mask
                            mask = selector.select_points(
                                data[0], data[1], data[2], sl
                            )
                        else:
                            f = pygio.PyGenericIO(data_file.filename )
                            coords = f.read(['x','y','z', 'mask','hh'], print_stats=False)
                            sl = coords['hh'] * 2.0
                            mask = selector.select_points(
                                coords['x'][si:ei], coords['y'][si:ei], coords['z'][si:ei], sl
                            )
                            pmask = FILTER_DEF[ptype](coords['mask'])[si:ei]

                        if mask is not None:
                            mask_sum = mask.sum()
                        if mask is None: continue
                        for field in field_list:
                            # print(si, ei, mask)
                            if field != 'hh':
                                if legacy:
                                    npart = ei-si
                                    data = pygio.read_globalOffset_scalar(data_file.filename, field, si, npart)
                                    # data = pygio.read_scalar(data_file.filename, field)[field][si:ei]
                                    data = np.array(data)[mask & pmask]
                                else:
                                    f = pygio.PyGenericIO(data_file.filename )
                                    data = f.read([field], print_stats=False)
                                data = np.array(data[field][si:ei])[mask & pmask]
                            else:
                                if legacy:
                                    npart = ei-si
                                    data = pygio.read_globalOffset_scalar(data_file.filename, field, si, npart)
                                    # data = pygio.read_scalar(data_file.filename, field)[field][si:ei] * 2.0
                                    data = np.array(data)[mask & pmask]    
                                else:
                                    f = pygio.PyGenericIO(data_file.filename )
                                    data = f.read(['hh'], print_stats=False)['hh'] * 2.0
                                    data = np.array(data[si:ei])[mask & pmask]                            
                                if field == 'hh' and ptype == 'Star':
                                    data = 1e-3 * data
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
        si, ei = data_file.start, data_file.end

        data_return = {}
        # open relevant file given by data_file.filename
        # probably have to open the subgrid file as well...  

        # iterate particle types and field list

            # iterate field_list