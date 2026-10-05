#!/usr/bin/env python3
"""Validate bundled national exposures; --download explicitly rebuilds public files.
No clinical ZIP list is read or used to construct the national public bundle.
"""
import argparse, concurrent.futures, datetime, hashlib, json, pathlib, subprocess, tempfile

REPO = 'petergraffy/environment_transplant_survival'
PRODUCTS = {
 'pm25': ('lghap-pm25-zcta-daily-v1', 'lghap_pm25_zcta_daily_', ['zip','date','pm25_ug_m3','value_source','fill_distance_m']),
 'o3': ('o3-zcta-daily-v1', 'o3_zcta_daily_', ['zip','date','o3_ppb','value_source','fill_distance_m']),
 'weather': ('gridmet-zcta-daily-v1', 'gridmet_zcta_daily_', ['zip','date','tmean_c','rhmean_pct']),
}
DEFAULT_ROOT = pathlib.Path('data/public/exposures')

def sha256(path):
 digest = hashlib.sha256()
 with pathlib.Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''): digest.update(block)
 return digest.hexdigest()

def curl(url, path):
 subprocess.run(['curl','-L','--fail','--silent','--show-error','--retry','3',url,'-o',str(path)],check=True)

def expected_months(first_year, last_year):
 # Only December of the preceding year is needed by the maximum 28-day lag.
 return [(year,month) for year in range(first_year-1,last_year+1)
         for month in range(1,13) if year>=first_year or month==12]

def validate(root, first_year=2018, last_year=2024):
 manifest = json.loads((root/'manifest.json').read_text())
 if manifest.get('scope') != 'national_public_no_clinical_filter':
  raise ValueError('Expected national public exposure bundle, not a clinical ZIP subset')
 for r in manifest['files']:
  if r['product'] not in PRODUCTS or not isinstance(r['year'],int) or not isinstance(r['month'],int) or not 1<=r['month']<=12:
   raise ValueError('Invalid public file identity')
  expected=f"{r['product']}/{r['product']}_{r['year']}_{r['month']:02d}.parquet"
  if r['path']!=expected:raise ValueError('Unsafe public file path')
 records = {(r['product'],r['year'],r['month']):r for r in manifest['files']}
 if len(records)!=len(manifest['files']):raise ValueError('Duplicate public month records')
 for source in manifest.get('sources',[]):
  if source['product'] not in PRODUCTS or not isinstance(source['year'],int):
   raise ValueError('Invalid public source identity')
 for product in PRODUCTS:
  for year,month in expected_months(first_year,last_year):
   rec = records.get((product,year,month))
   if rec is None: raise ValueError(f'Missing bundled month: {product} {year}-{month:02d}')
   path = root/rec['path']
   if not path.is_file() or sha256(path)!=rec['sha256']:
    raise ValueError(f'Missing/corrupt public file: {path}; explicitly rebuild with --download')
 return manifest

def prepare_year(product, year, release, root, first_year):
 import pyarrow as pa
 import pyarrow.compute as pc
 import pyarrow.parquet as pq
 tag,prefix,cols=PRODUCTS[product]
 asset=next(a for a in release['assets'] if a['name']==f'{prefix}{year}.parquet')
 metadata = root/product/f'{year}_source.json'
 if metadata.exists():
  old=json.loads(metadata.read_text())
  if old.get('asset_digest')==asset.get('digest') and all(
      (root/r['path']).exists() and sha256(root/r['path'])==r['sha256'] for r in old['files']):
   print(f'Already bundled {product} {year}',flush=True);return old
 (root/product).mkdir(parents=True,exist_ok=True)
 with tempfile.TemporaryDirectory(prefix='national-public-',dir=root) as tmp:
  source=pathlib.Path(tmp)/asset['name']
  print(f'Downloading national {product} {year} ({asset["size"]/1e6:.0f} MB)',flush=True)
  curl(asset['browser_download_url'],source)
  actual='sha256:'+sha256(source)
  if not asset.get('digest') or actual!=asset['digest']:
   raise ValueError(f'Unverified source digest: {asset["name"]}')
  table=pq.read_table(source,columns=cols)
  # Preserve values/types without rounding; only omit unused source columns.
  table=table.set_column(table.schema.get_field_index('zip'),'zip',pc.utf8_lpad(pc.cast(table['zip'],pa.string()),5,'0'))
  table=table.set_column(table.schema.get_field_index('date'),'date',pc.cast(table['date'],pa.date32()))
  files=[]
  for month in (range(1,13) if year>=first_year else [12]):
   start=datetime.date(year,month,1)
   end=datetime.date(year+1,1,1) if month==12 else datetime.date(year,month+1,1)
   part=table.filter(pc.and_(pc.greater_equal(table['date'],pa.scalar(start)),pc.less(table['date'],pa.scalar(end))))
   if not part.num_rows:raise ValueError(f'Empty national month: {product} {start}')
   part=part.sort_by([('zip','ascending'),('date','ascending')])
   distinct=part.select(['zip','date']).group_by(['zip','date']).aggregate([]).num_rows
   if distinct!=part.num_rows:raise ValueError('Duplicate public ZIP/date records')
   relative=f'{product}/{product}_{year}_{month:02d}.parquet';out=root/relative
   pq.write_table(part,str(out)+'.partial',compression='zstd',compression_level=9,row_group_size=100000)
   pathlib.Path(str(out)+'.partial').replace(out)
   if out.stat().st_size>=95_000_000:raise ValueError(f'Shard exceeds safe Git file size: {out}')
   files.append(dict(product=product,year=year,month=month,path=relative,rows=part.num_rows,
      n_zctas=pc.count_distinct(part['zip']).as_py(),size_bytes=out.stat().st_size,sha256=sha256(out)))
  record=dict(product=product,year=year,asset_url=asset['browser_download_url'],asset_digest=asset['digest'],
              verified_source_sha256=actual,columns=cols,files=files)
  metadata.write_text(json.dumps(record,indent=2)+'\n')
  print(f'Bundled national {product} {year}: {sum(r["size_bytes"] for r in files)/1e6:.1f} MB',flush=True)
  return record

def main():
 import os
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--run-dir');p.add_argument('--cache-dir',default=os.environ.get('MWAS_EXPOSURE_CACHE',str(DEFAULT_ROOT)))
 p.add_argument('--first-year',type=int,default=2018);p.add_argument('--last-year',type=int,default=2024)
 p.add_argument('--workers',type=int,default=2);p.add_argument('--download',action='store_true')
 a=p.parse_args();root=pathlib.Path(a.cache_dir)
 if a.download:
  root.mkdir(parents=True,exist_ok=True)
  with tempfile.TemporaryDirectory(prefix='releases-') as tmp:
   path=pathlib.Path(tmp)/'releases.json';curl(f'https://api.github.com/repos/{REPO}/releases?per_page=100',path)
   releases={r['tag_name']:r for r in json.loads(path.read_text())}
  jobs=[(prod,year,releases[info[0]]) for prod,info in PRODUCTS.items() for year in range(a.first_year-1,a.last_year+1)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
   sources=list(pool.map(lambda j:prepare_year(*j,root,a.first_year),jobs))
  manifest=dict(format='clif_mwas_national_exposures_v1',scope='national_public_no_clinical_filter',
   admission_years=[a.first_year,a.last_year],lookback_days=28,
   transformations='Required columns only, five-digit ZIP strings, date32, monthly ZIP/date sorted lossless ZSTD parquet; no value rounding/imputation',
   sources=sources,files=[r for source in sources for r in source['files']])
  (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 manifest=validate(root,a.first_year,a.last_year)
 if a.run_dir:
  (pathlib.Path(a.run_dir)/'exposure_manifest.json').write_text(json.dumps(dict(
   source=str(root),bundle_manifest_sha256=sha256(root/'manifest.json'),**manifest),indent=2)+'\n')
 print(f'Validated national public exposures: {len(manifest["files"])} monthly files in {root}',flush=True)
if __name__=='__main__':main()
