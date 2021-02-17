import iotbx.pdb
import os, sys
from libtbx import group_args
from libtbx import easy_run
import iotbx.ccp4_map
import mmtbx.model
import mmtbx.model
import mmtbx.utils
from libtbx import easy_pickle
from scitbx.array_family import flex
from libtbx import slots_getstate_setstate
from libtbx import adopt_init_args

pdb_files   = "/net/cci/pdb_mirror/pdb/"
mmcif_files = "/net/cci/pdb_mirror/mmcif/"
emdb        = "/net/cci-filer2/raid1/emdb_mirror/structures/"

work_root_folder = "/net/cci-filer2/raid1/share/cryoem/maps_and_models/"

def extract_from_header(header_path):
  pdb_codes = []
  resolutions = []
  for f in os.listdir(header_path):
    fo = open(header_path+f,"r")
    for l in fo.readlines():
      l=l.strip()
      if(l.count("<fittedPDBEntryId>")>0):
        pdb_code = l.replace("<fittedPDBEntryId>","").\
                     replace("</fittedPDBEntryId>","")
        pdb_codes.append(pdb_code)
      if(l.count("<resolutionByAuthor>")>0):
        resolution = l.replace("<resolutionByAuthor>","").\
                       replace("</resolutionByAuthor>","")
        resolutions.append(float(resolution))
      if(l.count("<resolution>")>0):
        resolution = l.replace("<resolution>","").replace("</resolution>","")
        resolutions.append(float(resolution))
    fo.close()
    pdb_codes   = list(set(pdb_codes))
    resolutions = list(set(resolutions))
  return pdb_codes, resolutions

def get_pdb_or_mmcif_file(files, pdb_code):
  pdb_file = None
  fo=open(files+"/INDEX","r")
  for pdb_file_ in fo.readlines():
    pdb_file_ = pdb_file_.strip()
    if(pdb_file_.count(pdb_code)):
      pdb_file = files+pdb_file_
      break
  fo.close()
  return pdb_file

def get_half_maps(other):
  candidates = []
  if(not os.path.isdir(other)): return None, None
  for f in os.listdir(other):
    if(f.count("half")==1):
      candidates.append(f)
  if(len(candidates) != 2): return None, None
  return read_map(map_file = other+candidates[0]), \
         read_map(map_file = other+candidates[1])

def read_map(map_file):
  map_inp = None
  error=None
  try:
    map_file_local = os.path.basename(map_file).strip(".gz")
    get_map(map_in=map_file, map_out=map_file_local)
    map_inp = iotbx.ccp4_map.map_reader(file_name=map_file_local)
  except Exception, e:
    error=str(e)
  return group_args(map_inp=map_inp, error=error, map_file=map_file)

def get_model(f):
  r = group_args(
    model      = None,
    error      = None,
    resolution = None,
    em         = None,
    date       = None)
  try:
    p = iotbx.pdb.input(file_name = f)
    m = mmtbx.model.manager(model_input = p)
    m.expand_with_BIOMT_records()
    resolution = p.resolution()
    em         = p.experiment_type_electron_microscopy()
    date       = p.deposition_date()
  except KeyboardInterrupt: raise
  except Exception, e:
    r.error = str(e)
    return r
  r.model      = m
  r.resolution = resolution
  r.em         = em
  r.date       = date
  return r

def get_map(map_in, map_out):
  files = [map_in, map_out, os.path.basename(map_in),
           os.path.basename(map_in).strip(".gz")]
  if(map_in.endswith(".gz")):
    #easy_run.call("gunzip -dc < %s > %s"%(map_in, map_out))
    # This is a test replacement for line above (faster?)
    easy_run.call("cp %s ."%map_in)
    easy_run.call("gunzip -f %s"%os.path.basename(map_in))
    assert os.path.isfile(map_out), files
  else:
    easy_run.call("cp %s %s"%(map_in, map_out))

def check_work_root_folder(work_root_folder):
  a1 = os.path.abspath(".")
  a2 = os.path.abspath(work_root_folder)
  assert a1 == a2, [a1, a2]

def get_folders_sorted_by_size():
  folders = []
  size = flex.double()
  for d in os.listdir(emdb):
    dm = emdb+d+"/map/"
    if(not os.path.isdir(dm)): continue
    map_file = dm+os.listdir(dm)[0]
    if(os.path.isfile(map_file)):
      folders.append(d)
      size.append(os.path.getsize(map_file))
  tmp = []
  for i in flex.sort_permutation(size):
    tmp.append(folders[i])
  return tmp

def err_and_stop(r, p, m):
  r.error = m
  easy_pickle.dump(
    file_name = work_root_folder + "%s/%s.pkl"%(p,p), obj = r)
  easy_pickle.load(work_root_folder + "%s/%s.pkl"%(p,p))

def map_not_consistent(mapo, mapo1, mapo2):
  m, m1, m2 = mapo.map_data(), mapo1.map_data(), mapo2.map_data()
  if(m.origin() != m1.origin() or
     m1.origin()!= m2.origin() or
     m.all() != m1.all() or
     m1.all()!= m2.all()): return True
  else: return False

def check_and_set_cs(model, mapo, mapo1, mapo2):
  map_inps= []
  cs = None
  error = None
  for _ in [mapo, mapo1, mapo2]:
    if(_ is not None):
      map_inps.append( _.map_inp )
  try:
    cs = mmtbx.utils.check_and_set_crystal_symmetry(
      models = [model], map_inps=map_inps)
  except KeyboardInterrupt: raise
  except Exception, e: error = str(e)
  if(cs is None or cs.is_empty() or cs.is_nonsense()):
    cs = None
    error = "Bad crystal symmetry"
  return group_args(cs = cs, error = error)

def result_template():
  ga = group_args(
    model_file = None,
    map_file   = None,
    map_file_1 = None,
    map_file_2 = None,
    resolution = None,
    date       = None,
    error      = None)
  ga.stop_dynamic_attributes()
  return ga

def run():
  #
  # Folder where this script is supposed to be executed
  check_work_root_folder(work_root_folder)
  #
  # Get list of folders like EMD-xxx
  folders = get_folders_sorted_by_size()
  #
  print "Using EMDB mirror:", emdb
  print "Total folders that contain map:", len(folders)
  #
  work_root_folder_content = os.listdir(work_root_folder)
  for i_emdb, d in enumerate(folders):
    sys.stdout.flush()
    root        = emdb+d
    map_path    = root+"/map/"
    header_path = root+"/header/"
    assert os.path.isdir(map_path)
    assert os.path.isdir(header_path)
    map_code = str(d.replace("EMD-",""))
    #
    # DEBUG
    #if map_code != "0077": continue # conatins half-maps
    #if map_code != "2842": continue
    #if map_code != "30021": continue
    #
    map_file = map_path + "emd_%s.map.gz"%map_code
    pdb_codes, map_resolutions = extract_from_header(header_path = header_path)
    assert os.path.isfile(map_file), [map_path,map_file] # Map file must exist.
    if(len(pdb_codes)==0): continue
    if(len(map_resolutions)==0): continue
    assert len(map_resolutions)>0, [pdb_codes, map_resolutions]
    other = root+"/other/"
    #
    for pdb_code in pdb_codes:
      prefix = pdb_code + "_" + map_code
      if(prefix in work_root_folder_content): continue
      cif_file = get_pdb_or_mmcif_file(files = mmcif_files, pdb_code = pdb_code)
      if(cif_file is None): continue # obsoleted?
      mo = get_model(f = cif_file)
      if(not mo.em): continue
      #
      # Record incidents from now on
      #
      result = result_template()
      os.makedirs(prefix)
      # Model extracted ?
      if(mo.model is None):
        err_and_stop(r=result, p=prefix, m="Failed input model: %s"%str(mo.error))
        continue
      else:
        result.date = mo.date
        result.model_file = cif_file
      # Consensus resolution
      if(mo.resolution in map_resolutions):
        result.resolution = mo.resolution
      else:
        err_and_stop(r=result, p=prefix, m="Resolution mismatch.")
        continue
      # Multi-model files
      if(len(mo.model.get_hierarchy().models())>1):
        err_and_stop(r=result, p=prefix, m="Cannot handle multi-model files")
        continue
      # Ignore single-atom models, threshold is arbitrary
      if(mo.model.percent_of_single_atom_residues()>25):
        err_and_stop(r=result, p=prefix, m="Ignore single-atom models")
        continue
      # Main map
      mapo = read_map(map_file=map_file)
      result.map_file = map_file
      if(mapo.map_inp is None):
        err_and_stop(r=result, p=prefix, m="Failed input map: %s"%str(mapo.error))
        continue
      # Half-maps
      mapo1, mapo2 = get_half_maps(other=other)
      if([mapo1, mapo2].count(None)==0):
        if([mapo1.map_inp, mapo2.map_inp].count(None)>0):
          err_and_stop(r=result, p=prefix, m="Failed input half-maps")
          continue
        else:
          result.map_file_1 = mapo1.map_file
          result.map_file_2 = mapo2.map_file
          # Map consistency (same grididng)
          if(map_not_consistent(
             mapo=mapo.map_inp, mapo1=mapo1.map_inp, mapo2=mapo2.map_inp)):
            err_and_stop(r=result, p=prefix,
                         m="Map and half-maps have different gridding")
            continue
      # Crystal symmetry
      cso = check_and_set_cs(
        model=mo.model, mapo=mapo, mapo1=mapo1, mapo2=mapo2)
      if(cso.cs is None):
        err_and_stop(r=result, p=prefix, m=cso.error)
        continue
      #
      # Write files
      #
      try:
        # Move full map
        t = os.path.basename(map_file).strip(".gz")
        d = work_root_folder + "%s/%s.map"%(prefix,prefix)
        easy_run.call("mv %s %s"%(t, d))
        # Move half-maps
        if([mapo1, mapo2].count(None)==0):
          for i, hmf in enumerate([mapo1.map_file, mapo2.map_file]):
            t = os.path.basename(hmf).strip(".gz")
            d = work_root_folder + "%s/%s_%d.map"%(prefix,prefix, i+1)
            easy_run.call("mv %s %s"%(t, d))
        # Write model
        d = work_root_folder + "%s/%s.pdb"%(prefix,prefix)
        with open(d,"w") as mfo:
          mfo.write(mo.model.model_as_pdb())
        d = work_root_folder + "%s/%s.cif"%(prefix,prefix)
        with open(d,"w") as mfo:
          mfo.write(mo.model.model_as_mmcif())
      except KeyboardInterrupt: raise
      except Exception, e:
        err_and_stop(r=result, p=prefix,
          m="Could not write maps and model: %s"%str(e))
        continue
      # Finally, write pkl with the summary
      easy_pickle.dump(
        file_name = work_root_folder + "%s/%s.pkl"%(prefix,prefix), obj=result)
    #
    easy_run.call("rm -rf %s*.map"%work_root_folder)

if (__name__ == "__main__"):
  run()
