from libtbx import easy_pickle
import iotbx.pdb
import sys, os
from libtbx.utils import null_out
import iotbx.ccp4_map
from scitbx.array_family import flex
from libtbx import easy_mp
from iotbx import map_and_model
from phenix.programs import validation_cryoem
import mmtbx.model

def run_one(arg):
  try:
    bug_log, pkl_result, pdb_file, map_file, map_file_1, map_file_2, \
      info_file = arg
    #
    map_inp   = iotbx.ccp4_map.map_reader(file_name = map_file)
    map_inp_1, map_inp_2 = None,None
    if(map_file_1 is not None):
      map_inp_1 = iotbx.ccp4_map.map_reader(file_name = map_file_1)
      map_inp_2 = iotbx.ccp4_map.map_reader(file_name = map_file_2)
    pdb_inp = iotbx.pdb.input(file_name = pdb_file)
    model = mmtbx.model.manager(
      model_input       = pdb_inp,
      stop_for_unknowns = False,
      log               = null_out())
    cs_consensus = mmtbx.utils.check_and_set_crystal_symmetry(
      models = [model], map_inps=[map_inp, map_inp_1, map_inp_2])
    map_data_1, map_data_2 = None,None
    if(map_file_1 is not None):
      map_data_1 = map_inp_1.map_data()
      map_data_2 = map_inp_2.map_data()
    base = map_and_model.input(
      map_data   = map_inp.map_data(),
      map_data_1 = map_data_1,
      map_data_2 = map_data_2,
      model      = model,
      box        = True)
    #
    inf = easy_pickle.load(info_file)
    resolutions = [inf.d_emdb, inf.d_pdb, inf.d_cif]
    while None in resolutions:
      resolutions.remove(None)
    resolution = min(resolutions)
    #
    params = validation_cryoem.master_params().extract()
    params.resolution = resolution
    params.scattering_table="n_gaussian"
    params.mtriage.include_curves = False
    params.mtriage.include_mask = False
    o = validation_cryoem.validation(
      model      = base.model(),
      map_data   = base.map_data(),
      map_data_1 = base.map_data_1(),
      map_data_2 = base.map_data_2(),
      params     = params).get_results()
    o.source_info = inf
    #
    easy_pickle.dump(file_name = pkl_result, obj = o)
  except Exception, e:
    of = open(bug_log, "w")
    for a in arg:
      print >> of, a
    print >> of, str(e)
    of.close()

def run(NPROC=10):
  path = "/net/cci/share/cryoem/maps_and_models/"
  bug_path = "/net/cci/share/cryoem/bugs/"
  args = []
  size = flex.double()
  for folder in os.listdir(path):
    prefix = folder
    folder = path+folder+"/"
    pkl_result = "%s%s.pkl"%(folder, prefix)
    bug_log = "%s%s.log"%(bug_path, prefix)
    if(not os.path.isdir(folder)): continue
    #
    if(os.path.isfile(pkl_result)):
      can_load=True
      try:
        easy_pickle.load(pkl_result)
      except:
        can_load=False
      if(can_load): continue
    #
    pdb_file   = folder+prefix+".pdb"
    map_file   = folder+prefix+".map"
    map_file_1 = folder+prefix+"_1.map"
    map_file_2 = folder+prefix+"_2.map"
    info_file  = folder+"source_info.pkl"
    assert os.path.isfile(pdb_file)
    assert os.path.isfile(map_file)
    assert os.path.isfile(info_file)
    if(not os.path.isfile(map_file_1)): map_file_1 = None
    if(not os.path.isfile(map_file_2)): map_file_2 = None
    arg = [bug_log, pkl_result, pdb_file, map_file, map_file_1, map_file_2,
           info_file]
    args.append(arg)
    size.append(easy_pickle.load(info_file).n_atoms)
  tmp = []
  for i in flex.sort_permutation(size):
    tmp.append(args[i])
  args = tmp[:]
  print "Total jobs:", len(args)
  sys.stdout.flush()
  #
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = args,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for arg in args:
      run_one(arg)
  return True

if (__name__ == "__main__"):
  run()
