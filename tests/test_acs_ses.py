import importlib.util, pathlib, unittest
spec = importlib.util.spec_from_file_location('acs',pathlib.Path(__file__).parents[1]/'code/26_cache_acs_zcta_ses.py')
acs = importlib.util.module_from_spec(spec);spec.loader.exec_module(acs)

class ACSChecks(unittest.TestCase):
    def test_denominators_and_missing(self):
        row = {'B17001_001E':'100','B17001_002E':'20','B23025_003E':'80',
               'B23025_004E':'72','B23025_005E':'8','B19013_001E':'60000',
               'B15003_001E':'100','B25014_001E':'100'}
        row.update({f'B15003_{n:03}E':'1' for n in range(2,17)})
        row.update({f'B25014_{n:03}E':'1' for n in (5,6,7,11,12,13)})
        out = acs.indicators(row)
        self.assertEqual(out['poverty_pct'],20)
        self.assertEqual(out['unemployment_pct'],10)  # civilian labor force, not population
        self.assertEqual(out['no_high_school_pct'],15)
        self.assertEqual(out['crowding_pct'],6)
        row['B17001_001E']='0';row['B19013_001E']='-666666666'
        self.assertIsNone(acs.indicators(row)['poverty_pct'])
        self.assertIsNone(acs.indicators(row)['median_household_income'])
        row['B19013_001E']='250001';row['B19013_001EA']='250,000+'
        self.assertIsNone(acs.indicators(row)['median_household_income'])

if __name__=='__main__':unittest.main()
