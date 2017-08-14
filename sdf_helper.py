import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.offsetbox import HPacker, VPacker, TextArea, AnchoredOffsetbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
try:
    import builtins
except:
    import __builtin__ as builtins

try:
    import sdf
    got_sdf = True
except ImportError:
    got_sdf = False
try:
    from matplotlib.pyplot import *  # NOQA
    got_mpl = True
except ImportError:
    got_mpl = False

# mpl.rcParams['interactive'] = True
# hold()

global data, t, step, p, ppi, ppe, rho, ei, ee, vx, vy, vz, bx, by, bz
global x, y, z, xc, yc, zc, grid, grid_mid, mult_x, mult_y
global old_mtime, wkdir, verbose, fig, im, cbar

verbose = True

wkdir = 'Data'
old_mtime = 0
old_size = 0
old_filename = ''
cached = False
fig = None
im = None
cbar = None
mult_x = 1
mult_y = 1

__version__ = "2.2.0"
_module_name = "sdf_helper"
_sdf_version = __version__


def _error_message():
    if not got_sdf:
        raise ImportError(
            r"This module relies on the sdf python module \n"
            "which we were unable to load\n")
    raise ImportError(
        r"Your sdf python module is too old for this version "
        "of " + _module_name + ".\n"
        "Either upgrade to sdf python " + _sdf_version + " or newer, or "
        "downgrade " + _module_name)


def _check_validity():
    if not got_sdf or not hasattr(sdf, "__version__"):
        return _error_message()
    our_version = list(map(int, _sdf_version.split(".")))
    lib_version = list(map(int, sdf.__version__.split(".")))
    # Check that major version number matches, and minor version is at least
    # as big as that specified
    if our_version[0] != lib_version[0]:
        return _error_message()
    if our_version[1] > lib_version[1]:
        return _error_message()

_check_validity()


class ic_type():
    NEW = 1
    RESTART = 2
    SOD = 3
    SALTZMAN = 4
    NOH = 5
    SEDOV = 6
    BRIO = 7
    WAVE = 8
    ADVECT = 9


def get_si_prefix(scale):
    scale = abs(scale)
    mult = 1
    sym = ''
    if scale > 1e-24:
        if scale < 1e-13:
            mult = 1e15
            sym = 'f'
        elif scale < 1e-10:
            mult = 1e12
            sym = 'p'
        elif scale < 1e-7:
            mult = 1e9
            sym = 'n'
        elif scale < 1e-4:
            mult = 1e6
            sym = '{\mu}'
        elif scale < 1e-1:
            mult = 1e3
            sym = 'm'
        elif scale >= 1e12:
            mult = 1e-12
            sym = 'T'
        elif scale >= 1e9:
            mult = 1e-9
            sym = 'G'
        elif scale >= 1e6:
            mult = 1e-6
            sym = 'M'
        elif scale >= 1e3:
            mult = 1e-3
            sym = 'k'
    return mult, sym


def get_title(geom=False):
    global data

    t = data.Header['time']
    mult, sym = get_si_prefix(t)

    stitle = r'$t = {:.3}{}s$'.format(mult * t, sym)

    if hasattr(data, 'Logical_flags'):
        if hasattr(data.Logical_flags, 'use_szp') \
                and data.Logical_flags.use_szp:
            stitle += r', $m_f = {:.1}$'.format(data.Real_flags.m_f)
        if hasattr(data.Logical_flags, 'use_tts') \
                and data.Logical_flags.use_tts:
            stitle += r', TTS'
        if hasattr(data.Logical_flags, 'use_tav') \
                and data.Logical_flags.use_tav:
            stitle += r', Tensor viscosity'
        else:
            if hasattr(data.Logical_flags, 'use_qmono') \
                    and data.Logical_flags.use_qmono:
                stitle += r', Edge viscosity'
            if hasattr(data.Logical_flags, 'use_edge') \
                    and data.Logical_flags.use_edge:
                stitle += r', Edge viscosity'

    if hasattr(data, 'Real_flags'):
        if hasattr(data.Real_flags, 'visc1'):
            stitle += r', $c_1 = ' + str(data.Real_flags.visc1) + '$'
        if hasattr(data.Real_flags, 'visc2'):
            stitle += r'$,\,c_2 = ' + str(data.Real_flags.visc2) + '$'

    if geom:
        if hasattr(data, 'Logical_flags'):
            if hasattr(data.Logical_flags, 'use_rz') \
                    and data.Logical_flags.use_rz:
                stitle += r', R-Z'
            else:
                stitle += r', Cartesian'

            if hasattr(data.Logical_flags, 'polar_grid') \
                    and data.Logical_flags.polar_grid:
                stitle += r' Polar'

    return stitle


def get_default_iso(data):
    iso = True
    if hasattr(data, 'Integer_flags') \
            and hasattr(data.Integer_flags, 'ic_type'):
        ic = data.Integer_flags.ic_type
        if ic == ic_type.NOH or ic == ic_type.SEDOV:
            iso = True
    return iso


def get_time(time=0, wkd=None):
    global data, wkdir

    if wkd is not None:
        wkdir = wkd

    flist = glob.glob(wkdir + "/[0-9]*.sdf")
    if len(flist) == 0:
        flist = glob.glob("./[0-9]*.sdf")
        if len(flist) == 0:
            print("No SDF files found")
            return
        wkdir = '.'

    t = None
    t_old = -1
    jobid = None
    fname = None
    for f in sorted(flist):
        dat_tmp = sdf.read(f)
        jobid_tmp = dat_tmp.Header['jobid1']
        if jobid is None:
            jobid = jobid_tmp
        elif jobid != jobid_tmp:
            continue

        t = dat_tmp.Header['time']
        if time is None:
            if t > t_old:
                fname = f
                t_old = t
        else:
            if t >= time - 1e-30:
                fname = f
                break

    if fname is None:
        raise Exception("No valid file found in directory: " + wkdir)

    data = getdata(fname, verbose=False)
    return data


def get_latest(wkd=None):
    return get_time(time=None, wkd=wkd)


def set_wkdir(wkd):
    global wkdir
    wkdir = wkd


def get_wkdir():
    global wkdir
    return wkdir


def sdfr(filename):
    return sdf.read(filename)


def oplot1d(*args, **kwargs):
    kwargs['set_ylabel'] = False
    kwargs['hold'] = True
    plot1d(*args, **kwargs)


def plot1d(var, fmt=None, xdir=None, idx=-1, xscale=0, yscale=0, cgs=False,
           title=True, sym=True, set_ylabel=True, hold=False, **kwargs):
    global data
    global x, y, mult_x, mult_y

    if len(var.dims) != 1 and len(var.dims) != 2:
        print("error: Not a 1d dataset")
        return

    if not hold:
        try:
            clf()
        except:
            pass

    if var.dims[0] == var.grid.dims[0]:
        grid = var.grid
    else:
        grid = var.grid_mid

    if len(var.dims) > 1:
        if xdir is None:
            if var.dims[1] < var.dims[0]:
                xdir = 0
            else:
                xdir = 1
        if xdir == 0:
            if idx == -1:
                idx = var.dims[1] / 2
            s = [slice(None), idx]
        else:
            if idx == -1:
                idx = var.dims[0] / 2
            s = [idx, slice(None)]
        Y = var.data[s]
    else:
        Y = var.data

    if np.ndim(var.grid.data[0]) == 1:
        X = grid.data[0]
    else:
        X = grid.data[xdir][s]

    if xdir is None:
        xdir = 0

    if xscale == 0:
        length = max(abs(X[0]),abs(X[-1]))
        mult_x, sym_x = get_si_prefix(length)
    else:
        mult_x, sym_x = get_si_prefix(xscale)

    if yscale == 0:
        length = max(abs(Y[0]),abs(Y[-1]))
        mult_y, sym_y = get_si_prefix(length)
    else:
        mult_y, sym_y = get_si_prefix(yscale)

    X = mult_x * X
    Y = mult_y * Y

    if fmt:
        plot(X, Y, fmt, **kwargs)
    else:
        plot(X, Y, **kwargs)

    xlabel(grid.labels[xdir] + ' $(' + sym_x + grid.units[xdir] + ')$')
    if set_ylabel:
        ylabel(var.name + ' $(' + sym_y + var.units + ')$')

    if title:
        plt.title(get_title(), fontsize='large', y=1.03)

    f = gcf()
    f.set_tight_layout(True)
    draw()


def oplot2d(*args, **kwargs):
    kwargs['hold'] = True
    plot2d(*args, **kwargs)


def plot2d(var, iso=None, fast=None, title=False, full=True, vrange=None,
           ix=None, iy=None, iz=None, reflect=0, norm=None, irange=None,
           jrange=None, hold=False, xscale=0, yscale=0):
    global data, fig, im, cbar
    global x, y, mult_x, mult_y

    si = slice(None, irange)
    sj = slice(None, jrange)

    i0 = 0
    i1 = 1
    if len(var.dims) == 3:
        if ix is not None:
            if iz is None:
                if ix < 0:
                    ix = var.dims[0] / 2
                i0 = 1
                i1 = 2
                ss = [ix,si,sj]
            else:
                if ix < 0:
                    ix = var.dims[2] / 2
                i0 = 0
                i1 = 2
                ss = [si,sj,ix]
        elif iy is not None:
            if iy < 0:
                iy = var.dims[1] / 2
            i0 = 0
            i1 = 2
            ss = [si,iy,sj]
            if iz is not None:
                i0 = 0
                i1 = 1
                ss = [si,iy,sj]
        elif iz is not None:
            if iz < 0:
                iz = var.dims[2] / 2
            i0 = 0
            i1 = 1
            ss = [si,sj,iz]
        else:
            print("error: Not a 2d dataset")
            return
        i2 = i0 + 3
        i3 = i1 + 3
    elif len(var.dims) != 2:
        print("error: Not a 2d dataset")
        return
    else:
        i2 = i0 + 2
        i3 = i1 + 2
        ss = [si,sj]

    var_data = var.data[ss]
    if np.ndim(x) == 1:
        x = var.grid.data[i0][si]
        y = var.grid.data[i1][sj]
    else:
        x = var.grid.data[i0][ss]
        y = var.grid.data[i1][ss]

    cmap = get_cmap()
    if norm is not None:
        v0 = np.min(var_data) - norm
        v1 = np.max(var_data) - norm
        if abs(v0/v1) > 1:
            low = 0
            high = 0.5 * (1 - v1/v0)
        else:
            low = 0.5 * (1 + v0/v1)
            high = 1.0

        cmap = cm.colors.LinearSegmentedColormap.from_list('tr',
                cmap(np.linspace(low,high,256)))

    if not hold:
        try:
            clf()
        except:
            pass

    if iso is None:
        iso = get_default_iso(data)

    ext = list(var.grid.extents)
    if xscale == 0:
        length = max(abs(ext[i2]),abs(ext[i0]))
        mult_x, sym_x = get_si_prefix(length)
    else:
        mult_x, sym_x = get_si_prefix(xscale)

    if yscale == 0:
        length = max(abs(ext[i3]),abs(ext[i1]))
        mult_y, sym_y = get_si_prefix(length)
    else:
        mult_y, sym_y = get_si_prefix(yscale)

    if vrange == 1:
        v = np.max(abs(var_data))
        vrange = [-v,v]

    if fast:
        if reflect == 1:
            # about x=0
            ext[i0] = 2 * ext[i0] - ext[i2]
            var_data = np.vstack((np.flipud(var_data), var_data))
        elif reflect == 2:
            # about y=0
            ext[i1] = 2 * ext[i1] - ext[i3]
            var_data = np.hstack((np.fliplr(var_data), var_data))
        elif reflect == 3:
            # about x=0, y=0
            ext[i0] = 2 * ext[i0] - ext[i2]
            ext[i1] = 2 * ext[i1] - ext[i3]
            var_data = np.vstack((np.flipud(var_data), var_data))
            var_data = np.hstack((np.fliplr(var_data), var_data))

    if np.ndim(x) == 1:
        if fast is None:
            fast = True

        if not fast:
            Y, X = np.meshgrid(y, x)
    else:
        if fast is None:
            fast = False

        if not fast:
            X = x
            Y = y

    if fast:
        ext[i0] = mult_x * ext[i0]
        ext[i1] = mult_y * ext[i1]
        ext[i2] = mult_x * ext[i2]
        ext[i3] = mult_y * ext[i3]
        e = ext[i1]
        ext[i1] = ext[i2]
        ext[i2] = e
        if vrange is None:
            im = imshow(var_data.T, interpolation='none', origin='lower',
                        extent=ext, cmap=cmap)
        else:
            im = imshow(var_data, interpolation='none', origin='lower',
                        extent=ext, cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    else:
        X = np.multiply(mult_x, X)
        Y = np.multiply(mult_y, Y)
        if vrange is None:
            im = pcolormesh(X, Y, var_data, cmap=cmap)
        else:
            im = pcolormesh(X, Y, var_data, cmap=cmap,
                            vmin=vrange[0], vmax=vrange[1])

    xlabel(var.grid.labels[i0] + ' $(' + sym_x + var.grid.units[i0] + ')$')
    ylabel(var.grid.labels[i1] + ' $(' + sym_y + var.grid.units[i1] + ')$')
    if full:
        plt.title(var.name + ' $(' + var.units + ')$, ' + get_title(),
                  fontsize='large', y=1.03)
    elif title:
        plt.title(var.name + ' $(' + var.units + ')$',
                  fontsize='large', y=1.03)

    axis('tight')
    if iso:
        axis('image')

    fig = plt.gcf()
    if not hold:
        ca = plt.gca()
        divider = make_axes_locatable(ca)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = colorbar(im, cax=cax)
        draw()
        plt.sca(ca)

    fig.set_tight_layout(True)
    draw()


def plot2d_update(var):
    global fig, im

    im.set_array(var.ravel())
    im.autoscale()
    fig.canvas.draw()


def plot_levels(var, r0=None, r1=None, nl=10, iso=None, out=False,
                title=True, levels=True):
    global data

    try:
        clf()
    except:
        pass

    if iso is None:
        iso = get_default_iso(data)

    if np.ndim(var.grid.data[0]) == 1:
        X, Y = np.meshgrid(var.grid.data[1], var.grid.data[0])
    else:
        if var.grid.dims == var.dims:
            X = var.grid.data[0]
            Y = var.grid.data[1]
        else:
            X = var.grid_mid.data[0]
            Y = var.grid_mid.data[1]

    if r0 is None:
        r0 = np.min(var.data)
        r1 = np.max(var.data)
        dr = (r1 - r0) / (nl + 1)
        r0 += dr
        r1 -= dr

    rl = r0 + 1.0 * (r1 - r0) * np.array(range(nl)) / (nl - 1)

    fig = gcf()
    if out:
        gs = GridSpec(1, 1)  # , width_ratios=[8, 1])
        ax = subplot(gs[0])
    else:
        ax = gca()

    cs = ax.contour(X, Y, var.data, levels=rl, colors='k', linewidths=0.5)

    if levels:
        fmt = {}
        for l, i in zip(cs.levels, range(1, len(cs.levels)+1)):
            fmt[l] = str(i)

        sidx = ""
        slvl = ""
        for l, i in reversed(zip(cs.levels, range(1, len(cs.levels)+1))):
            # sidx += rtn + "%i" % i
            # slvl += rtn + "%-6.4g" % l
            # rtn = "\n"
            sidx += "%i\n" % i
            slvl += "%-6.4g\n" % l

        t1 = TextArea('Level', textprops=dict(color='k', fontsize='small'))
        t2 = TextArea(sidx, textprops=dict(color='k', fontsize='small'))
        tl = VPacker(children=[t1, t2], align="center", pad=0, sep=2)

        lname = var.name.replace("_node", "")
        t3 = TextArea(lname, textprops=dict(color='k', fontsize='small'))
        t4 = TextArea(slvl, textprops=dict(color='k', fontsize='small'))
        tr = VPacker(children=[t3, t4], align="center", pad=0, sep=2)

        t = HPacker(children=[tl, tr], align="center", pad=0, sep=8)

        if out:
            t = AnchoredOffsetbox(loc=2, child=t, pad=.4,
                                  bbox_to_anchor=(1.01, 1), frameon=True,
                                  bbox_transform=ax.transAxes, borderpad=0)
        else:
            t = AnchoredOffsetbox(loc=1, child=t, pad=.4,
                                  bbox_to_anchor=(1, 1), frameon=True,
                                  bbox_transform=ax.transAxes, borderpad=.4)
        t.set_clip_on(False)
        ax.add_artist(t)

        clabel(cs, cs.levels, fmt=fmt, inline_spacing=2, fontsize=8)

    ax.set_xlabel(var.grid.labels[0] + ' $(' + var.grid.units[0] + ')$')
    ax.set_ylabel(var.grid.labels[1] + ' $(' + var.grid.units[1] + ')$')

    if title:
        if out:
            # suptitle(get_title(), fontsize='large')
            suptitle(get_title(), fontsize='large', y=0.92)
        else:
            plt.title(get_title(), fontsize='large', y=1.03)

    axis('tight')
    if iso:
        axis('image')

    draw()
    if out:
        gs.tight_layout(fig, rect=[0, -0.01, 0.95, 0.92])
        # fw = fig.get_window_extent().width
        # tw = fw*.06+t1.get_children()[0].get_window_extent().width \
        #     + t3.get_children()[0].get_window_extent().width
        # print(1-tw/fw)
        # gs.update(right=1-tw/fw)
        # ax.set_position([box.x0, box.y0, box.width + bw, box.height])
    else:
        fig.set_tight_layout(True)
    draw()


def plot_contour(var, r0=None, r1=None, nl=10, iso=None, title=True):
    return plot_levels(var, r0=r0, r1=r1, nl=nl, iso=iso, out=False,
                       title=title, levels=False)


def getdata(fname, wkd=None, verbose=True):
    global data, t, step, p, ppi, ppe, rho, ei, ee, vx, vy, vz, bx, by, bz
    global x, y, z, xc, yc, zc, grid, grid_mid
    global old_mtime, old_filename, old_size, cached
    global wkdir

    if wkd is not None:
        wkdir = wkd

    if isinstance(fname, int):
        filename = wkdir + "/%0.4i.sdf" % fname
    else:
        filename = fname

    try:
        st = os.stat(filename)
    except OSError as e:
        filename = "./%0.4i.sdf" % fname
        try:
            st = os.stat(filename)
            wkdir = '.'
        except OSError as e:
            print("ERROR opening file {0}: {1}".format(filename, e.strerror))
            raise

    if st.st_mtime != old_mtime or st.st_size != old_size \
            or filename != old_filename:
        if verbose:
            print("Reading file " + filename)
        data = sdf.read(filename)
        old_mtime = st.st_mtime
        old_size = st.st_size
        old_filename = filename
    else:
        cached = True
        return data

    cached = False
    fdict = {}

    sdfdict = {}
    for key, value in data.__dict__.items():
        if hasattr(value, "id"):
            sdfdict[value.id] = value
        else:
            sdfdict[key] = value

    table = {'time': 't'}
    k = 'Header'
    if k in sdfdict:
        h = sdfdict[k]
        k = list(table.keys())[0]
        if k in h:
            key = table[k]
            var = h[k]
            if verbose:
                print(key + str(np.shape(var)) + ' = ' + k)
            fdict[key] = var
            globals()[key] = var
            builtins.__dict__[key] = var

    table = {'Pressure': 'p',
             'Pressure_ion': 'ppi',
             'Pressure_electron': 'ppe',
             'Rho': 'rho',
             'Energy_ion': 'ei',
             'Energy_electron': 'ee',
             'Vx': 'vx',
             'Vy': 'vy',
             'Vz': 'vz',
             'Vr': 'vx',
             'VTheta': 'vz',
             'Bx': 'bx',
             'By': 'by',
             'Bz': 'bz',
             'Br': 'bx',
             'Bt': 'bz',
             'bx': 'bx',
             'by': 'by',
             'bz': 'bz',
             'ex': 'ex',
             'ey': 'ey',
             'ez': 'ez',
             'jx': 'jx',
             'jy': 'jy',
             'jz': 'jz'}

    rz = False
    if 'Vr' in sdfdict:
        rz = True

    if rz:
        table['Vz'] = 'vy'
        table['Bz'] = 'by'

    inv_table = {}
    for k, v in table.items():
        inv_table[v] = inv_table.get(v, [])
        inv_table[v].append(k)

    for k in table:
        if k in sdfdict:
            key = table[k]
            if hasattr(sdfdict[k], "data"):
                var = sdfdict[k].data
            else:
                var = sdfdict[k]
            dims = str(tuple(int(i) for i in sdfdict[k].dims))
            if verbose:
                print(key + dims + ' = ' + k)
            fdict[key] = var
            globals()[key] = var
            builtins.__dict__[key] = var

    k = 'grid'
    if k in sdfdict:
        vargrid = sdfdict[k]
        grid = vargrid
        keys = 'x', 'y', 'z'
        for n in range(np.size(vargrid.dims)):
            key = keys[n]
            var = vargrid.data[n]
            dims = str(tuple(int(i) for i in sdfdict[k].dims))
            if verbose:
                print(key + dims + ' = ' + k)
            fdict[key] = var
            globals()[key] = var
            builtins.__dict__[key] = var

    k = 'grid_mid'
    if k in sdfdict:
        vargrid = sdfdict[k]
        grid_mid = vargrid
        keys = 'xc', 'yc', 'zc'
        for n in range(np.size(vargrid.dims)):
            key = keys[n]
            var = vargrid.data[n]
            dims = str(tuple(int(i) for i in sdfdict[k].dims))
            if verbose:
                print(key + dims + ' = ' + k)
            fdict[key] = var
            globals()[key] = var
            builtins.__dict__[key] = var

    # Export particle arrays
    for k, value in data.__dict__.items():
        if type(value) != sdf.BlockPointVariable \
                and type(value) != sdf.BlockPointMesh:
            continue
        key = re.sub(r'[^a-z0-9]', '_', value.id.lower())
        if hasattr(value, "data"):
            var = value.data
        else:
            var = value
        dims = str(tuple(int(i) for i in value.dims))
        if type(value) == sdf.BlockPointVariable:
            if verbose:
                print(key + dims + ' = ' + value.name)
            fdict[key] = var
            globals()[key] = var
            builtins.__dict__[key] = var
        else:
            vargrid = value
            grid = vargrid
            keys = 'x', 'y', 'z'
            for n in range(np.size(value.dims)):
                gkey = keys[n] + '_' + key
                var = value.data[n]
                dims = str(tuple(int(i) for i in value.dims))
                if verbose:
                    print(gkey + dims + ' = ' + k + ' ' + keys[n])
                fdict[gkey] = var
                globals()[gkey] = var
                builtins.__dict__[gkey] = var

    # X, Y = np.meshgrid(x, y)
    return data


def ogrid(skip=None):
    global x, y, mult_x, mult_y
    if np.ndim(x) == 1:
        X, Y = np.meshgrid(x, y)
    else:
        s = slice(None, None, skip)
        X = x[s,s]
        Y = y[s,s]
    X = np.multiply(mult_x, X)
    Y = np.multiply(mult_y, Y)
    plot(X, Y, color='k', lw=0.5)
    plot(X.transpose(), Y.transpose(), color='k', lw=0.5, hold=True)


def plotgrid(fname=None, iso=None, title=True):
    if type(fname) is sdf.BlockList:
        dat = fname
    elif fname is not None:
        dat = getdata(fname, verbose=verbose)

    if iso is None:
        iso = get_default_iso(dat)

    ogrid()

    ax = gca()

    ax.set_xlabel(grid.labels[0] + ' $(' + grid.units[0] + ')$')
    ax.set_ylabel(grid.labels[1] + ' $(' + grid.units[1] + ')$')

    if title:
        plt.title(get_title(), fontsize='large', y=1.03)

    axis('tight')
    if iso:
        axis('image')

    draw()

    fig = gcf()
    fig.set_tight_layout(True)
    draw()


def axis_offset(boxed=False):
    ax = gca()
    xlab = ax.get_xlabel()
    ylab = ax.get_ylabel()

    f = 1e-3

    # for o in ax.findobj():
    for l in ax.get_lines():
        bb = l.get_clip_box()
        bb._bbox = Bbox([[-f, -f], [1+2*f, 1+2*f]])
        l.set_clip_box(bb)
        # l.set_clip_on(False)

    if boxed:
        r = matplotlib.patches.Rectangle((-f, -f), 1+2*f, 1+2*f,
                                         transform=ax.transAxes)
        r.set_color((0, 0, 0, 0))
        r.set_edgecolor('k')
        r.set_clip_on(False)
        ax.add_patch(r)

    w = 1.1
    gap = 8
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', gap))
    ax.spines['left'].set_linewidth(w)
    ax.spines['bottom'].set_position(('outward', gap))
    ax.spines['bottom'].set_linewidth(w)
    ax.tick_params(direction='outwards', width=w, length=4.5, top='off',
                   right='off')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    draw()


pi = 3.141592653589793238462643383279503
q0 = 1.602176565e-19
m0 = 9.10938291e-31
c = 2.99792458e8
kb = 1.3806488e-23
mu0 = pi * 4e-7
epsilon0 = 1.0 / mu0 / c**2
h_planck = 6.62606957e-34
h_bar = h_planck / 2.0 / pi
