import pygio
import glob
import numpy as np
from yt.utilities.io_handler import BaseIOHandler
from yt.frontends.sph.io import IOHandlerSPH


class HACCIOHandler(IOHandlerSPH):
    # _particle_reader = True
    _dataset_type = "hacc_genericio"


    def _count_particles(self, data_file):
        print('Counting particles: %s'%data_file.filename)
        print('si, ei = %d, %d'%(data_file.start, data_file.end))
        si, ei = data_file.start, data_file.end
        f = pygio.PyGenericIO('%s'%data_file.filename)
        nparts = np.array([f.read_num_elems()])
        if None not in (si, ei):
            np.clip(nparts-si, 0, ei-si, out=nparts)
        nparts = {'PartType0':nparts}
        return nparts

    def _identify_fields(self, data_file):
        f = pygio.PyGenericIO(data_file.filename)
        fields = []
        for var in f.get_variables():
            if var.name != 'status':
                fields.append(('PartType0', var.name))
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
            print('read_particle_coords: fname = %s'%data_file.filename)
            with pygio.PyGenericIO(data_file.file_name) as f:
                for ptype in sorted(ptf):
                    
                    coords = f.read(['x', 'y', 'z'])
                    data = coords['x'][si, ei], coords['y'][si, ei], coords['z'][si, ei]


                    yield ptype, data
                
    def _yield_coordinates(self, data_file, needed_ptype=None):
        f = pygio.PyGenericIO(data_file.filename)
        for ptype, cnt in data_file.total_particles.items():
            if cnt == 0: continue
            if needed_ptype is not None and ptype != needed_ptype:
                continue
            pp = f.read(['x','y','z'])
            pp = [ax for ax in pp.values()]
            pp = np.array(pp).T
            # print(pp.shape)
            yield ptype, pp

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
                    si, ei = data_file.start, data_file.end[0]
                    # print(type(si))
                    # print(type(ei))
                    data_return = {}
                    f = pygio.PyGenericIO(data_file.filename)
                    for ptype, field_list in sorted(ptf.items()):
                        if data_file.total_particles[ptype] == 0:
                            continue
                        if selector is None or getattr(selector, 'is_all_data', False):
                            mask = slice(None, None, None)
                            mask_sum = data_file.total_particles[ptype]
                        else:
                            coords = f.read(['x','y','z'])
                            sl = f.read(['hh'])['hh'] * 2.0
                            mask = selector.select_points(
                                coords['x'][si:ei], coords['y'][si:ei], coords['z'][si:ei], sl
                            )
                            if mask is not None:
                                mask_sum = mask.sum()
                        if mask is None: continue
                        for field in field_list:
                            # print(si, ei, mask)
                            if field != 'hh':
                                data = f.read([field])
                                data = np.array(data[field][si:ei])[mask]
                            else:
                                data = f.read(['hh'])['hh'] * 2.0
                                data = np.array(data[si:ei])[mask]                            
                        # data_return[(ptype, field)] = data
                    
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
        # open relevant file given by data_file.file_name
        # probably have to open the subgrid file as well...  

        # iterate particle types and field list

            # iterate field_list