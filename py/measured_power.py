import os
import pickle

def get_repo_dir():
    """ Return directory of this repository, using environment variable """

    assert ('SHERWOOD' in os.environ),'Define SHERWOOD environ'
    return os.environ['SHERWOOD']


def get_power_fname(grid,power_type):
    """ Find out (full path) filename for a given measured power spectrum """

    # get full path to repo data
    data_dir=get_repo_dir()+'/data/'
    print('data dir',data_dir)

    # get nametag from grid
    nametag=grid.get_nametag()
    print('name tag',nametag)

    if power_type == "p1d":
        return '{}/flux_p1d/p1d_{}.p'.format(data_dir,nametag)
    elif power_type == "p3d":
        return '{}/flux_p3d/p3d_{}_20_16_20.p'.format(data_dir,nametag)
    elif power_type == "cross":
        return '{}/cross_p3d/cross_{}_20_16_20.p'.format(data_dir,nametag)
    elif power_type == "halo":
        # only one mu bin when not adding RSDs
        n_mu = (16 if grid.add_rsd else 1)
        return '{}/halo_p3d/halo_{}_20_{}_20.p'.format(data_dir,nametag,n_mu)
    else:
        raise ValueError("unknown power spectrum type",power_type)


def get_measured_power(grid=None,power_type=None,fname=None):
    """Return measured power spectrum from pickled file"""

    if fname is None:
        fname = get_power_fname(grid,power_type)
    print('will unpickle from file',fname)
    return pickle.load(open(fname,"rb"))

