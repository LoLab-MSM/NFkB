# import pandas as pd
# import glob
#
# path =r'/Users/geenaildefonso/Projects/NFkB/Hoffman/PSO_results' # use your path
# allFiles = glob.glob(path + "/*.csv")
# frame = pd.DataFrame()
# list_ = []
# for file_ in allFiles:
#     df = pd.read_csv(file_,index_col=None, header=0)
#     list_.append(df)
# frame = pd.concat(list_)
# print(len(frame))
#
# import pandas as pd
# import glob, os
# df = pd.concat(map(pd.read_csv, glob.glob(os.path.join('', "best_fit_*.csv"))))
# df.to_csv('happy.csv')

from old_hoffman_nfkb import model

def read_all_pars(pars_path, new_path=None):
    """
    Reads all pars in file or directory
    :param new_path:
    :param pars_path: Parameter file or directory path
    :return: DataFrame with all the parameters
    """
    if type(pars_path) is list:
        par_sets = pars_path
    elif os.path.isfile(pars_path):
        par_sets = list_pars_infile(pars_path, new_path)
    elif os.path.isdir(pars_path):
        par_sets = listdir_fullpath(pars_path)
    else:
        raise Exception("Not valid parameter file or path")
    pars_0 = read_pars(par_sets[0])
    all_pars = np.zeros((len(par_sets), len(pars_0)))
    all_pars[0] = pars_0
    for idx in range(1, len(par_sets)):
        all_pars[idx] = read_pars(par_sets[idx])
    return all_pars

read_all_pars()