#!/usr/bin/env python3
"""Cache national public ACS ZCTA data. No clinical ZIP list is transmitted."""
import argparse, concurrent.futures, csv, hashlib, io, json, math, os, pathlib, re, subprocess, zipfile
from urllib.parse import urlencode

TABLES = {
    'B17001': [1, 2],
    'B15003': list(range(1, 17)),
    'B23025': [3, 4, 5],
    'B19013': [1],
    'B25014': [1, 5, 6, 7, 11, 12, 13],
}

def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) and x >= 0 else None
    except (ValueError, TypeError):
        return None

def estimate(row, base):
    # Do not silently treat bounded/annotated medians or suppressed values as exact.
    annotation = row.get(base + 'EA')
    return None if annotation not in (None, '', 'null') else number(row.get(base + 'E'))

def total(row, bases):
    values = [estimate(row, b) for b in bases]
    return sum(values) if all(x is not None for x in values) else None

def percent(num, den):
    return 100 * num / den if num is not None and den is not None and den > 0 and 0 <= num <= den else None

def indicators(row):
    poverty_den = estimate(row, 'B17001_001')
    education_den = estimate(row, 'B15003_001')
    labor_den = estimate(row, 'B23025_003')
    housing_den = estimate(row, 'B25014_001')
    return dict(
        poverty_pct=percent(estimate(row, 'B17001_002'), poverty_den),
        no_high_school_pct=percent(total(row, [f'B15003_{n:03}' for n in range(2, 17)]), education_den),
        unemployment_pct=percent(estimate(row, 'B23025_005'), labor_den),
        median_household_income=estimate(row, 'B19013_001'),
        crowding_pct=percent(total(row, [f'B25014_{n:03}' for n in (5, 6, 7, 11, 12, 13)]), housing_den),
        poverty_denominator=poverty_den, education_denominator=education_den,
        labor_force_denominator=labor_den, occupied_housing_denominator=housing_den)

def fetch(url, path):
    if not path.exists():
        tmp = path.with_suffix('.partial')
        try:
            # Keep an optional API key out of process arguments and exceptions.
            subprocess.run(['curl','--config','-','-L','--fail','--silent','--show-error',
                            '--retry','3','--max-time','120','-o',str(tmp)],
                           input='url = '+json.dumps(url)+'\n',text=True,check=True)
        except subprocess.CalledProcessError:
            raise RuntimeError('Census download failed; check connectivity and API-key configuration') from None
        # Validate before atomically installing; HTML/API errors are not cached.
        text = tmp.read_text()
        if '<title>Missing Key</title>' in text:
            raise RuntimeError('Census API data queries require CENSUS_API_KEY; use --source summary instead')
        data = json.loads(text)
        if not isinstance(data, (dict, list)):
            raise ValueError('Unexpected Census response')
        tmp.replace(path)
    return json.loads(path.read_text())

def bulk_download(url, path):
    if not path.exists():
        tmp = path.with_suffix('.partial')
        subprocess.run(['curl','-L','--fail','--silent','--show-error','--retry','3',
                        '--max-time','120',url,'-o',str(tmp)],check=True)
        tmp.replace(path)

def summary_records(root, year):
    base = f'https://www2.census.gov/programs-surveys/acs/summary_file/{year}'
    folder = base + '/data/5_year_seq_by_state/UnitedStates/All_Geographies_Not_Tracts_Block_Groups'
    lookup = root / 'sequence_lookup.txt'
    lookup_url = base + '/documentation/user_tools/ACS_5yr_Seq_Table_Number_Lookup.txt'
    geo = root / f'g{year}5us.csv'
    geo_url = folder + '/' + geo.name
    bulk_download(lookup_url, lookup); bulk_download(geo_url, geo)
    layouts = {}
    for row in csv.DictReader(lookup.open(encoding='latin1')):
        if row['Table ID'] in TABLES and row['Start Position']:
            layouts[row['Table ID']] = (row['Sequence Number'].zfill(4),int(row['Start Position'])-1)
    if set(layouts) != set(TABLES):
        raise ValueError('Incomplete summary-file layout')
    log_to_zcta = {}
    for row in csv.reader(geo.open(encoding='latin1')):
        if row[2]=='860' and row[3]=='00':
            geoid = [x for x in row if re.fullmatch(r'86000US\d{5}',x)]
            if len(geoid)!=1:
                raise ValueError('Unexpected ZCTA geography row')
            log_to_zcta[row[4].zfill(7)] = geoid[0][-5:]
    if not log_to_zcta:
        raise ValueError('No national ZCTA geographies')
    records = {z:{} for z in log_to_zcta.values()}
    tasks = [(folder+f'/{year}5us{seq}000.zip',root/f'{year}5us{seq}000.zip')
             for seq in sorted({seq for seq,_ in layouts.values()})]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda t:bulk_download(*t),tasks))
    for url,path in tasks:
        with zipfile.ZipFile(path) as archive:
            for prefix,suffix in [('e','E'),('m','M')]:
                names = [n for n in archive.namelist() if pathlib.Path(n).name.startswith(prefix) and n.endswith('.txt')]
                if len(names)!=1:raise ValueError('Unexpected sequence ZIP members')
                with archive.open(names[0]) as f:
                    for row in csv.reader(io.TextIOWrapper(f,encoding='latin1')):
                        zcta = log_to_zcta.get(row[5].zfill(7))
                        if not zcta:continue
                        for table,(table_seq,start) in layouts.items():
                            if table_seq!=row[4].zfill(4):continue
                            for n in TABLES[table]:
                                records[zcta][f'{table}_{n:03}{suffix}'] = row[start+n-1]
    for row in records.values():
        if number(row.get('B19013_001E')) in (2499,250001):
            row['B19013_001EA'] = 'bounded'
    return records,[(lookup_url,lookup),(geo_url,geo)]+tasks

def validate_bundle(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    allowed = {'zcta_ses.csv','acs_variables.csv'} | {table + '_metadata.json' for table in TABLES}
    for name, expected in manifest['bundled_files_sha256'].items():
        if name not in allowed: raise ValueError('Unsafe ACS bundle path')
        path = root / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Missing/corrupt bundled ACS file: {path}')
    return manifest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--year', type=int, default=2017)
    parser.add_argument('--cache', default='data/public/acs')
    parser.add_argument('--source-cache', default='data/mwas_cache/acs')
    parser.add_argument('--source', choices=['summary','api'],default='summary')
    parser.add_argument('--download', action='store_true', help='Explicit maintainer rebuild; ordinary site runs are offline')
    args = parser.parse_args()
    root = pathlib.Path(os.environ.get('MWAS_ACS_DIR', str(pathlib.Path(args.cache) / str(args.year))))
    if not args.download:
        validate_bundle(root)
        print(f'Validated bundled national ACS: {root}', flush=True)
        return
    root.mkdir(parents=True, exist_ok=True)
    source_root = pathlib.Path(args.source_cache) / str(args.year)
    source_root.mkdir(parents=True, exist_ok=True)
    api = f'https://api.census.gov/data/{args.year}/acs/acs5'
    metadata_tasks = [(api+f'/groups/{table}.json',source_root/f'{table}_metadata.json') for table in TABLES]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda t:fetch(*t),metadata_tasks))
    if args.source=='summary':
        records,tasks = summary_records(source_root,args.year)
    else:
        records,tasks = api_records(source_root,api)
    tasks += metadata_tasks
    for _, path in metadata_tasks:
        (root / path.name).write_bytes(path.read_bytes())
    write_indicators(root,args,records,tasks)
    validate_bundle(root)

def api_records(root,api):
    key = os.environ.get('CENSUS_API_KEY')
    if not key:raise RuntimeError('Set CENSUS_API_KEY or use --source summary')
    tasks = []
    for table, numbers in TABLES.items():
        variables = [f'{table}_{n:03}{suffix}' for n in numbers for suffix in ('E', 'M', 'EA', 'MA')]
        # Census API allows at most 50 requested variables; leave room for NAME.
        for i in range(0, len(variables), 44):
            url = api + '?' + urlencode({'get': ','.join(['NAME'] + variables[i:i+44]),
                                          'for': 'zip code tabulation area:*'})
            tasks.append((url, root / f'{table}_{i//44}.json'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda t: fetch(t[0]+'&'+urlencode({'key':key}),t[1]), tasks))
    records = {}
    for result in results:
        if isinstance(result, dict):
            continue
        header, *rows = result
        for values in rows:
            row = dict(zip(header, values))
            zcta = row['zip code tabulation area'].zfill(5)
            records.setdefault(zcta, {}).update(row)
    return records,tasks

def write_indicators(root,args,records,tasks):
    output = root / 'zcta_ses.csv'
    data = [dict(zcta=z, acs_year=args.year, **indicators(records[z])) for z in sorted(records)]
    with output.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader(); writer.writerows(data)
    raw = root / 'acs_variables.csv'
    fields = sorted({key for row in records.values() for key in row})
    with raw.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['zcta'] + fields)
        writer.writeheader()
        writer.writerows(dict(zcta=z, **records[z]) for z in sorted(records))
    bundled = [output, raw] + sorted(root.glob('*_metadata.json'))
    manifest = dict(bundled_files_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in bundled}, acs_year=args.year, period=f'{args.year-4}-{args.year}', geography='ZCTA',source=args.source,
        zcta_count=len(data), urls=[url for url, _ in tasks],
        sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for _, p in tasks},
        indicator_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        definitions=dict(poverty='B17001_002 / B17001_001',
            no_high_school='sum B15003_002..016 / B15003_001; age 25+',
            unemployment='B23025_005 / B23025_003; civilian labor force',
            income='B19013_001; median household income in vintage-year dollars',
            crowding='sum B25014_005,006,007,011,012,013 / B25014_001; >1 occupants/room'),
        missing='Negative sentinels, annotations, invalid ratios and zero denominators -> missing; no imputation',
        uncertainty='Required raw estimates/MOE and available annotations preserved in acs_variables.csv; ACS sampling error not propagated into interactions',
        linkage='Exact five-digit residential ZIP-to-ZCTA proxy; unmatched codes remain missing; no boundary crosswalk')
    (root / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Cached {len(data)} national ZCTAs: {output}', flush=True)

if __name__ == '__main__':
    main()
