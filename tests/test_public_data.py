"""Offline integrity, coverage and national-scope tests; no clinical data/network."""
import hashlib, importlib.util, json, pathlib, tempfile, unittest

ROOT=pathlib.Path(__file__).parents[1]
def load(name, file):
    spec=importlib.util.spec_from_file_location(name,ROOT/'code'/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
exposure=load('exposure','22_cache_mwas_exposures.py')
acs=load('acs','26_cache_acs_zcta_ses.py')

class PublicBundleChecks(unittest.TestCase):
    def fixture(self, root):
        files=[]
        for product in exposure.PRODUCTS:
            for year,month in exposure.expected_months(2018,2018):
                name=f'{product}/{product}_{year}_{month:02d}.parquet';p=root/name;p.parent.mkdir(exist_ok=True)
                p.write_bytes(b'synthetic public file')
                files.append(dict(product=product,year=year,month=month,path=name,sha256=exposure.sha256(p)))
        manifest=dict(scope='national_public_no_clinical_filter',files=files)
        (root/'manifest.json').write_text(json.dumps(manifest));return manifest
    def test_full_year_and_lookback_required(self):
        self.assertEqual(exposure.expected_months(2018,2018)[0],(2017,12))
        self.assertEqual(len(exposure.expected_months(2018,2024)),85)
        with tempfile.TemporaryDirectory() as folder:
            root=pathlib.Path(folder);manifest=self.fixture(root)
            exposure.validate(root,2018,2018)
            manifest['files'].pop(0);(root/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'Missing bundled month'):exposure.validate(root,2018,2018)
    def test_corruption_and_clinical_subset_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=pathlib.Path(folder);manifest=self.fixture(root)
            (root/manifest['files'][0]['path']).write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'Missing/corrupt'):exposure.validate(root,2018,2018)
            manifest['scope']='site_zip_subset';(root/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'national public'):exposure.validate(root,2018,2018)
    def test_manifest_cannot_add_private_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            root=pathlib.Path(folder);manifest=self.fixture(root)
            manifest['files'][0]['path']='../../private/cohort.parquet'
            (root/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'Unsafe public file path'):exposure.validate(root,2018,2018)
    def test_acs_source_variables_integrity(self):
        with tempfile.TemporaryDirectory() as folder:
            root=pathlib.Path(folder);p=root/'acs_variables.csv';p.write_bytes(b'public ACS estimates and MOE')
            (root/'manifest.json').write_text(json.dumps(dict(bundled_files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()})))
            acs.validate_bundle(root);p.unlink()
            with self.assertRaisesRegex(ValueError,'Missing/corrupt'):acs.validate_bundle(root)

if __name__=='__main__':unittest.main()
