"""Real-H5AD regression cases for coordinate coercion."""
import contextlib
import io
import json
import anndata as ad
import numpy as np
import pytest
import validate_submission as vs
@pytest.mark.parametrize('coords,accepted',[
 (np.arange(9).reshape(3,3).astype(str),True),
 (np.full((3,3),'bad'),False),
 (np.array([['1','2','3'],['4','bad','6'],['7','8','9']]),False),
 (np.array([['1','2','3'],['4','nan','6'],['7','8','9']]),False),
 (np.array([['1','2','3'],['4','1e40','6'],['7','8','9']]),False),
 (np.arange(9,dtype=float).reshape(3,3),True),
 (np.column_stack([np.arange(9).reshape(3,3),[np.nan,np.inf,1e40]]),True),
 (np.column_stack([np.arange(9).reshape(3,3).astype(str),['bad']*3]),False),
 (np.zeros((3,2)),False),
 (np.zeros((3,3,1)),False),
])
def test_coordinate_coercion(tmp_path,coords,accepted):
 a=ad.AnnData(X=np.log1p(np.tile([4000.,6000.],(3,1))).astype(np.float32))
 a.var_names=['g1','g2'];a.obsm['spatial_3D']=coords
 p=tmp_path/'input.h5ad';a.write_h5ad(p)
 (tmp_path/'genes.txt').write_text('g1\ng2\n')
 (tmp_path/'panels.json').write_text(json.dumps({'test':{'label':'synthetic','min_cells':1,'max_cells':3,'genes_file':'genes.txt','obsm_required':['spatial_3D'],'needs_coords':True}}))
 with contextlib.redirect_stdout(io.StringIO()):
  rc=vs.main(['--board','test','--input',str(p),'--panels',str(tmp_path/'panels.json')])
 assert rc==(0 if accepted else 1)
