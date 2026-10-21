import contextlib
import io
import json
import warnings
import anndata as ad
import numpy as np
import pytest
import validate_submission as vs

@pytest.mark.parametrize('got,panel,failed',[
 (['A','B','C','A'],['A','B','C'],True),
 (['A','B'],['A','B','A'],True),
 (['A','B','C','D'],['A','B','C'],True),
 (['A','B'],['A','B','C'],True),
 (['A','A','B','C'],['A','B','C'],True),
 (['B','A','C'],['A','B','C'],True),
 (['A','B','C'],['A','B','C'],False),
 (['A','B','B'],['A','B','C'],True),
])
def test_gene_contract_stays_strict(got,panel,failed):
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  a=ad.AnnData(np.ones((2,len(got)),dtype=np.float32)); a.var_names=got
  r=vs.Report(); out=io.StringIO()
  with contextlib.redirect_stdout(out): vs.check_genes(r,a,panel)
 assert bool(r.failed)==failed
 if len(got)!=len(panel):
  assert f'基因数 {len(got)}，应为 {len(panel)}' in out.getvalue()


def test_appended_duplicate_h5ad_returns_failure_exit(tmp_path,capsys):
 panel=['A','B','C']
 a=ad.AnnData(np.ones((2,4),dtype=np.float32)); a.var_names=panel+['A']
 with warnings.catch_warnings():
  warnings.simplefilter('ignore'); a.write_h5ad(tmp_path/'pred.h5ad')
 (tmp_path/'genes.txt').write_text('\n'.join(panel))
 (tmp_path/'panels.json').write_text(json.dumps({'TEST':{'label':'synthetic','genes_file':'genes.txt','min_cells':1,'max_cells':10,'obsm_required':[],'needs_coords':False}}))
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  code=vs.main(['--board','TEST','--input',str(tmp_path/'pred.h5ad'),'--panels',str(tmp_path/'panels.json')])
 assert code==1
 out=capsys.readouterr().out
 assert '基因数 4，应为 3' in out
 assert '不合格' in out
