#!/usr/bin/env python3
"""Download public yearly releases, verify GitHub asset SHA256, retain required ZCTAs.
Local ZIP linkage file never leaves this computer. Network requests contain only public asset URLs.
"""
import argparse, concurrent.futures, hashlib, json, pathlib, subprocess, tempfile
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = 'petergraffy/environment_transplant_survival'
PRODUCTS = {
 'pm25': ('lghap-pm25-zcta-daily-v1', 'lghap_pm25_zcta_daily_', ['zip','date','pm25_ug_m3','value_source','fill_distance_m']),
 'o3': ('o3-zcta-daily-v1', 'o3_zcta_daily_', ['zip','date','o3_ppb','value_source','fill_distance_m']),
 'weather': ('gridmet-zcta-daily-v1', 'gridmet_zcta_daily_', ['zip','date','tmean_c','rhmean_pct']),
}

def curl(url, path):
 subprocess.run(['curl','-L','--fail','--silent','--show-error','--retry','3',url,'-o',str(path)],check=True)

def cache_one(product, year, release, zips, root):
 tag,prefix,cols=PRODUCTS[product]
 asset=next(a for a in release['assets'] if a['name']==f'{prefix}{year}.parquet')
 out=root/product/f'{product}_{year}.parquet'
 meta=out.with_suffix('.json')
 ziphash=hashlib.sha256('\n'.join(sorted(zips)).encode()).hexdigest()
 if out.exists() and meta.exists():
  prior=json.loads(meta.read_text())
  if prior.get('zip_set_sha256')==ziphash and prior.get('asset_digest')==asset.get('digest'):
   print(f'Cached {product} {year}',flush=True);return prior
 out.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.TemporaryDirectory(prefix='mwas-public-',dir=root) as tmp:
  source=pathlib.Path(tmp)/asset['name']
  print(f'Downloading {product} {year} ({asset["size"]/1e6:.0f} MB)',flush=True)
  curl(asset['browser_download_url'],source)
  digest=hashlib.sha256()
  with source.open('rb') as f:
   for block in iter(lambda:f.read(8*1024*1024),b''):digest.update(block)
  actual='sha256:'+digest.hexdigest()
  if asset.get('digest') and actual!=asset['digest']:raise RuntimeError(f'Checksum mismatch: {source.name}')
  writer=None;count=0
  try:
   for batch in pq.ParquetFile(source).iter_batches(batch_size=250000,columns=cols):
    tab=pa.Table.from_batches([batch]);tab=tab.filter(pc.is_in(pc.cast(tab['zip'],pa.string()),value_set=pa.array(zips)))
    if not tab.num_rows:continue
    if writer is None:writer=pq.ParquetWriter(str(out)+'.partial',tab.schema,compression='zstd')
    writer.write_table(tab);count+=tab.num_rows
  finally:
   if writer:writer.close()
  if not count:raise RuntimeError(f'No matching ZCTAs: {product} {year}')
  pathlib.Path(str(out)+'.partial').replace(out)
 record={'product':product,'year':year,'asset_url':asset['browser_download_url'],
         'asset_digest':asset.get('digest'),'verified_sha256':actual,'zip_set_sha256':ziphash,
         'n_zips_requested':len(zips),'rows_cached':count}
 meta.write_text(json.dumps(record,indent=2));print(f'Ready {product} {year}: {count} rows',flush=True)
 return record

def main():
 p=argparse.ArgumentParser();p.add_argument('--run-dir',required=True);p.add_argument('--cache-dir',default='data/mwas_cache')
 p.add_argument('--first-year',type=int,default=2017);p.add_argument('--last-year',type=int,default=2024)
 p.add_argument('--workers',type=int,default=2);a=p.parse_args()
 import csv
 zips=[r['zip'] for r in csv.DictReader(open(pathlib.Path(a.run_dir)/'private/required_zips.csv'))]
 root=pathlib.Path(a.cache_dir);root.mkdir(parents=True,exist_ok=True)
 releasefile=root/'releases.json';curl(f'https://api.github.com/repos/{REPO}/releases?per_page=100',releasefile)
 releases={r['tag_name']:r for r in json.loads(releasefile.read_text())}
 jobs=[(prod,year,releases[info[0]]) for prod,info in PRODUCTS.items() for year in range(a.first_year,a.last_year+1)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  results=list(pool.map(lambda j:cache_one(*j,zips,root),jobs))
 (pathlib.Path(a.run_dir)/'exposure_manifest.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
