import os,sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp
import numpy as np
import scanpy as sc

x=sc.read_csv("../data/work/Hep_for.scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("../data/work/pySCENIC/hep_50k.loom",x.X.transpose(),row_attrs,col_attrs);

x=sc.read_csv("../data/work/EC_for.scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("../data/work/pySCENIC/LSECs.loom",x.X.transpose(),row_attrs,col_attrs);

x=sc.read_csv("../data/work/HSC_for.scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("../data/work/pySCENIC/HSC.loom",x.X.transpose(),row_attrs,col_attrs);

