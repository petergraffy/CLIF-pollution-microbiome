#!/usr/bin/env python3
"""Package an allowlisted workflow; optionally include verified national public data."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ['.gitattributes','.gitignore','config/.gitignore','renv/.gitignore','README.md','LICENSE','requirements-buddy.txt','renv/buddy.lock',
        'docs/buddy_testing.md','code/README.md','config/config_template.json','config/README.md',
        'data/public/README.md','resources/mwas/clif_intermittent_med_categories.csv','utils/config.R','utils/clif_io.R']
SCRIPTS = ['21_prepare_acute_mwas.R','22_cache_mwas_exposures.py','23_run_acute_mwas.R','26_cache_acs_zcta_ses.py',
           '27_run_federated_mwas.R','28_pool_federated_mwas.R','29_report_federated_mwas.R','30_audit_duration_candidates.R',
           '31_shared_exposure_checks.R','32_calibrate_mwas.R','33_site_preflight.R','34_run_buddy_site.R',
           '35_buddy_smoke_test.R','36_package_buddy_source.py','37_site_characteristics.R','38_audit_modifier_support.R']
UTILS = ['mwas.R','mwas_federated.R','mwas_poisson.R','mwas_power.R','mwas_severity.R','mwas_targets.R','mwas_characteristics.R']
TESTS = ['test_mwas.R','test_mwas_power.R','test_mwas_windows.R','test_mwas_federated.R','test_mwas_poisson.R',
         'test_mwas_characteristics.R','test_acs_ses.py','test_public_data.py','smoke_buddy_pipeline.R']
FILES = sorted(BASE + ['code/'+s for s in SCRIPTS] + ['utils/'+s for s in UTILS] + ['tests/'+s for s in TESTS])

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def module(name, filename):
    spec=importlib.util.spec_from_file_location(name,ROOT/'code'/filename)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

def public_files():
    exposures=ROOT/'data/public/exposures'
    manifest=module('exposures','22_cache_mwas_exposures.py').validate(exposures)
    acs=ROOT/'data/public/acs/2017'
    acs_manifest=module('acs','26_cache_acs_zcta_ses.py').validate_bundle(acs)
    files=['data/public/exposures/manifest.json','data/public/acs/2017/manifest.json']
    files += ['data/public/exposures/'+r['path'] for r in manifest['files']]
    files += [f'data/public/exposures/{r["product"]}/{r["year"]}_source.json' for r in manifest['sources']]
    files += ['data/public/acs/2017/'+n for n in acs_manifest['bundled_files_sha256']]
    return files

def build(include_public_data=False):
    files=sorted(set(FILES+(public_files() if include_public_data else [])))
    checksums={}
    for name in files:
        path=ROOT/name
        if path.is_symlink() or not path.is_file() or ROOT.resolve() not in path.resolve().parents:
            raise RuntimeError(f'Missing or unsafe bundle file: {name}')
        checksums[name]=digest(path)
    manifest={'format':'clif_mwas_buddy_v2','includes_national_public_data':include_public_data,
              'files_sha256':checksums,'excludes':['clinical data','local config','outputs','site caches','Git history']}
    out=ROOT/'dist';out.mkdir(exist_ok=True)
    archive=out/('clif_mwas_buddy_with_public_data.zip' if include_public_data else 'clif_mwas_buddy_source.zip')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,allowZip64=True) as z:
        for name in files:
            info=zipfile.ZipInfo(name,date_time=(2026,10,5,0,0,0))
            info.compress_type=zipfile.ZIP_STORED if name.endswith('.parquet') else zipfile.ZIP_DEFLATED
            info.external_attr=0o100644<<16
            with (ROOT/name).open('rb') as f,z.open(info,'w',force_zip64=True) as target:
                for block in iter(lambda:f.read(8*1024*1024),b''):target.write(block)
        info=zipfile.ZipInfo('BUNDLE_MANIFEST.json',date_time=(2026,10,5,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
        z.writestr(info,json.dumps(manifest,indent=2)+'\n')
    checksum=digest(archive)
    archive.with_suffix('.zip.sha256').write_text(f'{checksum}  {archive.name}\n')
    print(f'Buddy bundle: {archive} ({len(files)} files, SHA256 {checksum})')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--include-public-data',action='store_true')
    build(parser.parse_args().include_public_data)
